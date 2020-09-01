"""Download genome annotations, transform for visualizations, upload to GCS
"""

import json
import os
import subprocess
import urllib.request as request

import genome_annotation_metadata
import utils


def get_ensembl_metadata():
    """Get organism, assembly, and annotation release metadata from Ensembl
    """
    ensembl_metadata = {}

    # API docs: https://rest.ensembl.org/documentation/info/species
    url = 'https://rest.ensembl.org/info/species?content-type=application/json'
    with request.urlopen(url) as response:
        data = response.read().decode('utf-8')
    ensembl_species = json.loads(data)['species']

    for species in ensembl_species:
        taxid = species['taxon_id']
        name = species['name']
        assembly = species['assembly']
        strain = species['strain']

        if taxid == '10090' and strain != 'reference (CL57BL6)':
            # Mouse has multiple annotated assemblies; only use reference assembly
            continue

        ensembl_metadata[taxid] = {
            'organism': name,
            'taxid': taxid,
            'assembly_name': assembly,
            'assembly_accession': species['accession'],
            'release': str(species['release']),
        }

    return ensembl_metadata


class GenomeAnnotations(object):
    def __init__(
        self,
        vault_path='',
        input_dir='./',
        local_output_dir='output/',
        # TODO (SCP-2490): Integrate API endpoint for reference bucket
        gcs_bucket='fc-bcc55e6c-bec3-4b2e-9fb2-5e1526ddfcd2',
        remote_output_dir='reference_data_dev/',
        remote_prod_dir='reference_data/',
        use_cache=True,
        scp_species=None,
    ):
        """Download genome annotations, transform for visualizations

        Args:

        vault_path: Path in Vault for GCS service account credentials
        input_dir: Input directory; where to find organisms.tsv
        local_output_dir: Local directory for output
        gcs_bucket: Name of Google Cloud Storage bucket for upload
        remote_output_dir: Remote directory for output in GCS bucket
        remote_prod_dir: Remote directory from which to copy data into
            remote_output_dir.  Use to ensure test data environment is equivalent
            to production data environment.
        use_cache: Whether to use cache
        scp_species: List of lists of species identifiers
        """

        if scp_species is None:
            self.scp_species = utils.get_species_list(input_dir + 'organisms.tsv')
        else:
            # As used in make_toy_data.py
            self.scp_species = scp_species

        self.output_dir = local_output_dir
        self.remote_output_dir = remote_output_dir

        if os.path.exists(local_output_dir) is False:
            os.makedirs(local_output_dir)

        self.context = {
            'vault_path': vault_path,
            'gcs_bucket': gcs_bucket,
            'output_dir': local_output_dir,
            'remote_prod_dir': remote_prod_dir,
            'remote_output_dir': remote_output_dir,
        }

        self.ensembl_metadata = get_ensembl_metadata()
        self.fetch_gtfs(output_dir=self.output_dir)
        self.make_local_reference_dirs()

        self.transform()

    def load(self):
        """Upload transformed GTFs to Google Cloud Storage
        """
        ensembl_metadata = genome_annotation_metadata.upload_ensembl_gtf_products(
            self.ensembl_metadata, self.scp_species, self.context
        )
        genome_annotation_metadata.record_annotation_metadata(
            self.output_dir, ensembl_metadata, self.scp_species
        )

    def get_ensembl_gtf_urls(self, output_dir):
        """Construct the URL of an Ensembl genome annotation GTF file.

        Returns GTF URLs, and sets local GTF path in Ensembl metadata

        Example URL:
        http://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz
        """

        gtf_urls = []
        for species in self.scp_species:
            taxid = species[2]
            organism_metadata = self.ensembl_metadata[taxid]
            release = organism_metadata['release']
            organism = organism_metadata['organism']
            organism_upper = organism[0].upper() + organism[1:]
            assembly = organism_metadata['assembly_name']

            origin = 'http://ftp.ensembl.org'
            dir = '/pub/release-' + release + '/gtf/' + organism + '/'
            filename = organism_upper + '.' + assembly + '.' + release + '.gtf.gz'

            gtf_url = origin + dir + filename
            gtf_urls.append(gtf_url)
            self.ensembl_metadata[taxid]['gtf_path'] = os.path.join(
                output_dir, filename
            )

        return gtf_urls

    def transform_ensembl_gtf(self, gtf_path, ref_dir):
        """Produce sorted GTF and GTF index from Ensembl GTF; needed for igv.js
        """
        # Example:
        # $ sort -k1,1 -k4,4n gencode.vM17.annotation.gtf > gencode.vM17.annotation.possorted.gtf
        # $ bgzip gencode.vM17.annotation.possorted.gtf
        # $ tabix -p gff gencode.vM17.annotation.possorted.gtf.gz
        sorted_filename = gtf_path.replace('.gtf', '.possorted.gtf')
        sorted_filename = ref_dir + sorted_filename.replace(self.output_dir, '')
        outputs = [sorted_filename + '.gz', sorted_filename + '.gz.tbi']
        if os.path.exists(outputs[1]):
            print('  Using cached GTF transforms')
            return outputs
        else:
            print('  Producing GTF transforms for ' + gtf_path)

        # sort by chromosome name, then genomic start position; needed for index
        sort_command = ('sort -k1,1 -k4,4n ' + gtf_path).split(' ')
        sorted_file = open(sorted_filename, 'w')
        subprocess.call(sort_command, stdout=sorted_file)

        # bgzip enables requesting small indexed chunks of a gzip'd file
        bgzip_command = ('bgzip ' + sorted_filename).split(' ')
        subprocess.call(bgzip_command)

        # tabix creates an index for the GTF file, used for getting small chunks
        tabix_command = ('tabix -p gff ' + sorted_filename + '.gz').split(' ')
        subprocess.call(tabix_command)

        return outputs

    def make_local_reference_dirs(self):
        """Create a folder hierarchy on this machine to mirror that planned for GCS
        """
        print('Making local reference directories')

        for species in self.scp_species:
            taxid = species[2]
            organism_metadata = self.ensembl_metadata[taxid]
            organism = organism_metadata['organism']
            asm_name = organism_metadata['assembly_name']
            asm_acc = organism_metadata['assembly_accession']
            release = 'ensembl_' + organism_metadata['release']
            folder = os.path.join(self.output_dir, self.remote_output_dir)
            folder += organism + '/' + asm_name + '_' + asm_acc + '/' + release + '/'

            os.makedirs(folder, exist_ok=True)

            self.ensembl_metadata[taxid]['reference_dir'] = folder

    def fetch_gtfs(self, output_dir=''):
        """Download GTF files in parallel
        """

        gtf_urls = self.get_ensembl_gtf_urls(output_dir)

        print('Fetching GTFs')
        gtfs = utils.batch_fetch(gtf_urls, output_dir)
        print('Got GTFs!  Number: ' + str(len(gtfs)))

        return gtfs

    def transform(self):
        """Transform each local raw Ensembl GTF to position-sorted GTF and index
        """
        transformed_gtfs = []

        print('Transforming GTFs')
        for species in self.scp_species:
            taxid = species[2]
            organism_metadata = self.ensembl_metadata[taxid]
            gtf_path = organism_metadata['gtf_path'].replace('.gz', '')
            ref_dir = organism_metadata['reference_dir']
            transformed_gtfs = self.transform_ensembl_gtf(gtf_path, ref_dir)

            self.ensembl_metadata[taxid]['transformed_gtfs'] = transformed_gtfs

        return self.ensembl_metadata
