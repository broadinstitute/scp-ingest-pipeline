"""Fetch assembly reports from NCBI and write filtered list to disk

See constructor docstring for more detail.
"""

import os

import utils


class GenomeAssemblies(object):
    def __init__(self, input_dir='./', output_dir='output/'):
        """ Fetch assembly reports from NCBI and write filtered list to disk
        Parse genome assemblies from NCBI by their level of completion.

        This class fetches current and historical Assembly Reports from NCBI, and
        extracts metadata on assemblies that are relevant for Single Cell Portal.
        See get_assembly_report docstring for more detail on historical Assembly
        Reports.

        Only assemblies that are chromosome-level or above are included.

        For definitions of assembly-related terms and other domain background, see:
        https://www.ncbi.nlm.nih.gov/assembly/help/#glossary

        Args:
        * input_dir: Directory where to find organisms.tsv
        * output_dir: Directory to send output data to.  Default: output/
        """

        self.input_dir = input_dir
        self.output_dir = output_dir

        self.species_list = utils.get_species_list(input_dir + 'organisms.tsv')

        if os.path.exists(output_dir) is False:
            os.makedirs(output_dir)

        assembly_report = self.get_assembly_report()
        historical_assembly_report = self.get_assembly_report(use_historical=True)

        assemblies = self.parse_assembly_report(assembly_report)
        historical_assemblies = self.parse_assembly_report(historical_assembly_report)

        assemblies += historical_assemblies

        # Sort assemblies by species priority and assembly release date
        assemblies = sorted(assemblies, key=lambda assembly: assembly[4], reverse=True)

        # One-dimensional list of species names
        species_names = [species[0] for species in self.species_list]

        assemblies = sorted(
            assemblies, key=lambda assembly: species_names.index(assembly[0])
        )

        self.write_assemblies_to_file(assemblies)

    def get_assembly_report(self, use_historical=False):
        """ Fetch NCBI Assembly Report, or read it from disk if present

        :param use_historical {Boolean} Whether to use NCBI's historical Assembly
        Report.  Needed to get major assembly versions when patches are available
        (e.g. human, mouse).
        """
        domain = 'https://ftp.ncbi.nlm.nih.gov'
        path = '/genomes/ASSEMBLY_REPORTS/'
        if use_historical:
            report_filename = 'assembly_summary_genbank_historical.txt'
        else:
            report_filename = 'assembly_summary_genbank.txt'
        assembly_report_url = domain + path + report_filename
        assembly_report_path = os.path.join(self.output_dir, report_filename)

        # To consider: Refactor fetch_content so this is less gross.
        assembly_report = utils.fetch_content(
            [[assembly_report_url, assembly_report_path]]
        )
        assembly_report = list(assembly_report.values())[0]

        return assembly_report

    def parse_columns(self, columns):
        """ Extract relevant values from columns of the NCBI Assembly Report
        """
        accession = columns[0]
        refseq_category = columns[4]
        # taxid = columns[5]
        # species_taxid = columns[6] # Helps with dog, for example
        organism_name = columns[7]
        assembly_level = columns[11]
        release_type = columns[12]
        # genome_rep = columns[13]
        release_date = columns[14].replace('/', '-')
        assembly_name = columns[15]
        refseq_accession = columns[17]

        return [
            release_type,
            assembly_level,
            refseq_category,
            assembly_name,
            release_date,
            organism_name,
            accession,
            refseq_accession,
        ]

    def is_relevant_assembly(self, rel_type, asm_level, refseq_category, refseq_acc):
        """ Determine if assembly is suitable for SCP
        """

        if (rel_type in 'Major' and asm_level == 'Chromosome') is False:
            # Exclude patch assemblies, and genomes that lack assembled chromosomes
            return False

        if refseq_category == 'na':
            # Encountered only in historical assembly report
            if refseq_acc == 'na':
                return False
        else:
            if (
                refseq_category
                in ('representative genome', 'reference genome')
                is False
            ):
                return False

        return True

    def update_assemblies(self, columns, assemblies):
        """ Parse and return all relevant assemblies, by organism
        """
        (
            rel_type,
            asm_level,
            refseq_category,
            asm_name,
            rel_date,
            org_name,
            accession,
            refseq_acc,
        ) = self.parse_columns(columns)

        if self.is_relevant_assembly(rel_type, asm_level, refseq_category, refseq_acc):
            if org_name == 'Homo sapiens' and asm_name[:3] != 'GRC':
                # Omit CHM1_1.1, HuRef, etc.
                return assemblies
            assembly = [asm_name, accession, rel_date]
            if org_name not in assemblies:
                assemblies[org_name] = [assembly]
            else:
                assemblies[org_name].append(assembly)
        return assemblies

    def refine_assemblies(self, assemblies):
        """ Filter assemblies to only those for relevant species; add species data
        """
        refined_assemblies = []

        for species in self.species_list:
            scientific_name, common_name, taxid = species
            if scientific_name not in assemblies:
                continue
            asms = assemblies[scientific_name]
            for asm in asms:
                (assembly_name, accession, release_date) = asm
                row = [
                    scientific_name,
                    common_name,
                    taxid,
                    assembly_name,
                    accession,
                    release_date,
                ]
                refined_assemblies.append(row)

        return refined_assemblies

    def parse_assembly_report(self, assembly_report):
        """ Read Assembly Report, filter to relevant assemblies
        """
        assemblies = {}
        for line in assembly_report:
            if len(line) == 0 or line[0] == '#':
                # Skip header lines
                continue
            columns = line.split('\t')

            assemblies = self.update_assemblies(columns, assemblies)

        refined_assemblies = self.refine_assemblies(assemblies)

        return refined_assemblies

    def write_assemblies_to_file(self, assemblies):
        """ Save assembly metadata to disk
        """
        assemblies_str = []
        header = [
            '# scientific_name',
            'common_name',
            'taxid',
            'assembly_name',
            'assembly_accession',
            'assembly_release_date',
        ]
        assemblies.insert(0, header)
        for assembly in assemblies:
            assemblies_str.append('\t'.join(assembly))
        assemblies_str = '\n'.join(assemblies_str)

        output_path = os.path.join(self.output_dir, 'species_metadata_reference.tsv')
        with open(output_path, 'w') as f:
            f.write(assemblies_str)
        print('Wrote assemblies to ' + output_path)
