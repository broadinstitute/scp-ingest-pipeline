"""Extract, transform, and load data into TSV files for to rank genes

Ranked genes are all, at least, mentioned in a study's related publication.
These published genes are then sought among very differentially expressed
genes in the annotation's various labeled groups.  The groups in which each
published gene is very DE'd (if any) are noted.  Finally, data on each gene's
worldwide popularity, ultimately measured via the Gene Hints pipeline using
either Wikipedia page views or PubMed citation counts, is also included.

The output TSV file contains one gene per row. This is parsed by Image Pipeline
to prioritize which genes to cache.

EXAMPLES

python3 ingest/ingest_pipeline.py --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 rank_genes --study-accession SCP367 --publication="https://doi.org/10.1126/sciadv.adf6251" --rank-genes
"""

import logging
import re
import datetime
import glob
import gzip
import tarfile
import json
import urllib
import urllib.request
import csv
import os
from io import BytesIO
import xml.etree.ElementTree as ET

import requests

# Used when importing internally and in tests
from ingest_files import IngestFiles
from monitor import setup_logger

timestamp = datetime.datetime.now().isoformat(sep="T", timespec="seconds")
url_safe_timestamp = re.sub(':', '', timestamp)
log_name = f"rank_genes_{url_safe_timestamp}_log.txt"

logger = setup_logger(
    __name__, log_name, level=logging.INFO, format="support_configs"
)

def download_gzip(url, output_path, cache=0):
    """Download gzip file, decompress, write to output path; use optional cache
    Cached files can help speed development iterations by > 2x, and some
    development scenarios (e.g. on a train or otherwise without an Internet
    connection) can be impossible without it.
    """
    is_targz = url.endswith("tar.gz")
    headers = {"Accept-Encoding": "gzip"} if is_targz else {}
    logger.info(f"Requesting {url}")
    request = urllib.request.Request(url, headers=headers)
    response = urllib.request.urlopen(request)
    if is_targz:
        tar = tarfile.open(name=None, fileobj=BytesIO(response.read()))
        tar.extractall(output_path)
        tar.close()
    else:
        content = gzip.decompress(response.read()).decode()
        with open(output_path, "w") as f:
            f.write(content)


def fetch_pmcid(doi):
    """Convert Digital Object Identifier (DOI) into PubMed Central ID (PMCID)

    Many PMC articles are fully public-access, and machine-readable.
    """
    idconv_base = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
    params = "?" + "&".join(
        [
            "tool=scp-fetch-pmcid",
            "email=scp-dev@broadinstitute.org",
            "format=json",
            f"ids={doi}",
        ]
    )
    idconv_url = idconv_base + params
    with urllib.request.urlopen(idconv_url) as response:
        data = response.read().decode("utf-8")
    idconv_json = json.loads(data)
    record = idconv_json["records"][0]
    if "pmcid" not in record:
        msg = f"PubMed Central ID (PMCID) not found for DOI \"{doi}\""
        raise ValueError(msg)
    pmcid = record["pmcid"]
    return pmcid


def fetch_pmcid_text(pmcid):
    """Get full text for publication, given its PubMed Central ID"""
    oa_url = f"https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id={pmcid}"
    with urllib.request.urlopen(oa_url) as response:
        data = response.read().decode("utf-8")
    oa_xml = ET.fromstring(data)
    oa_package = oa_xml.find(".//link[@format='tgz']")
    if oa_package is None:
        msg = f"This PMCID is not open access: {pmcid}"
        raise ValueError(msg)
    package_url = oa_package.attrib["href"].replace("ftp://", "https://")
    download_gzip(package_url, ".")
    nxml_path = glob.glob(f"{pmcid}/*.nxml")[0]
    article = ET.parse(nxml_path).getroot()
    text = str(ET.tostring(article.find("body"), method="text"))
    return text


def fetch_publication_text(publication):
    """Get full text for a publicly-accessible article given its URL"""
    # Examples:
    # <Is open access> -- <publication URL> -- <Study accession>
    # 0. Yes -- https://www.biorxiv.org/content/10.1101/2021.11.13.468496v1 -- SCP1671
    # ---
    # 1. Not really, but abstract mentions genes -- https://doi.org/10.1210/en.2018-00750 -- SCP2058
    # 2. Yes -- https://www.biorxiv.org/content/10.1101/2022.07.01.498017v1 -- SCP2053
    # 3. Yes -- https://www.biorxiv.org/content/10.1101/2022.09.27.509354v1 -- SCP2046
    # 4. Yes -- https://www.biorxiv.org/content/10.1101/2022.09.27.509354v1 -- SCP2045
    # 5. Yes -- https://www.biorxiv.org/content/10.1101/2022.09.27.509354v1 -- SCP2011 (summary body)
    # 6. Yes (via PMC) -- https://doi.org/10.1038/s41467-022-28372-y -- SCP1985
    # 7. Yes (via biorxiv) -- https://doi.org/10.1101/2022.01.09.475571 -- SCP1921
    # 8. Yes -- https://www.biorxiv.org/content/10.1101/2022.06.29.496888v1 -- SCP1905
    # 9. Yes -- https://www.biorxiv.org/content/10.1101/2022.06.29.496888v1 -- SCP1903
    # 10. Not really, but abstract mentions genes -- https://doi.org/10.1016/j.immuni.2023.01.002 -- SCP1884
    # 11. Not really, not until 2023-10-12, but abstract mentions genes, and full text available in biorxiv -- https://doi.org/10.1093/toxsci/kfac109 -- SCP1875
    doi_split = publication.split("https://doi.org/")
    if len(doi_split) > 1:
        doi = doi_split[1]
        pmcid = fetch_pmcid(doi)
        logger.info('Converted DOI to PMC ID: ' + pmcid)
        text = fetch_pmcid_text(pmcid)
    elif publication.startswith("https://www.biorxiv.org"):
        full_text_url = f"{publication}.full.txt"
        logger.info('full_text_url', full_text_url)
        with urllib.request.urlopen(full_text_url) as response:
            text = response.read().decode('utf-8')
    else:
        file_path = publication
        with open(file_path) as f:
            text = f.read()
    return text


def fetch_gene_cache(organism):
    genes_filename = f"{organism}-genes.tsv"  # e.g. homo-sapiens-genes.tsv
    base_url = "https://cdn.jsdelivr.net/npm/"
    genes_url = f"{base_url}ideogram@1.41.0/dist/data/cache/genes/{genes_filename}.gz"
    download_gzip(genes_url, genes_filename)

    # Populate containers for various name and significance fields
    interest_rank_by_gene = {}
    counts_by_gene = {}
    loci_by_gene = {}
    full_names_by_gene = {}

    i = 0
    with open(genes_filename) as file:
        reader = csv.reader(file, delimiter="\t")
        for row in reader:
            if row[0] == '#' or len(row) < 2:
                continue
            gene = row[4]
            full_name = row[5]
            counts_by_gene[gene] = 0
            loci_by_gene[gene] = {
                "chromosome": row[0],
                "start": row[1],
                "length": row[2],
            }
            full_names_by_gene[gene] = full_name
            # Genes in upstream file are ordered by global popularity
            interest_rank_by_gene[gene] = i
            i += 1

    return [interest_rank_by_gene, counts_by_gene, loci_by_gene, full_names_by_gene]

def screen_false_positive_mentions(publication_text, mentions_by_gene):
    """Screen out hits that are non-gene abbreviations
    """
    fp_for_tf = bool(re.search("transcription factor(s)? \(TF", publication_text))
    if fp_for_tf and "TF" in mentions_by_gene:
        del mentions_by_gene["TF"]
        logger.info(
            'Deleted false positive gene mention for "TF", ' +
            'defined as "transcription factor" in publication'
        )

    # Detect gene-mentions that use have the grammatical structure:
    #
    #   the <term> (<"gene">)
    #
    # like "the Arcuate (ARC", which mentions a non-gene.
    # Source context: https://doi.org/10.1126/sciadv.adf6251
    articular_fp_genes = []
    for gene in mentions_by_gene:
        fp_match =\
            re.search("the ([A-Za-z ]){0,20} \(" + gene, publication_text)
        is_likely_fp = bool(fp_match)
        if is_likely_fp:
            term = fp_match.group(0)
            if (
                "express" in term or # e.g. "the Oxytocin expressing (OXT"
                gene.lower() not in term.split('(')[0].lower() # e.g. "the astrocyte (GJA1"
            ):
                continue
            logger.info(
                f'Deleted false positive gene mention for "{gene}", ' +
                f'defined as "{term}" in publication'
            )
            articular_fp_genes.append(gene)

    for fp_gene in articular_fp_genes:
        del mentions_by_gene[fp_gene]

    return mentions_by_gene

def extract_mentions_and_interest(organism, publication):
    """Get mentions in publication and interest rank for each gene in organism
    """
    logger.info('Extracting mentions and interest')
    [
        interest_rank_by_gene,
        counts_by_gene,
        loci_by_gene,
        full_names_by_gene,
    ] = fetch_gene_cache(organism)

    publication_text = fetch_publication_text(publication)

    publication_words = re.split(r'[ /]', publication_text)
    for word in publication_words:
        raw_word = word.strip('*,.()-+')
        if raw_word in counts_by_gene:
            counts_by_gene[raw_word] += 1

    mentions_by_gene = {}
    for gene in counts_by_gene:
        counts_by_gene[gene]
        count = counts_by_gene[gene]
        if counts_by_gene[gene] > 0:
            mentions_by_gene[gene] = count

    mentions_by_gene =\
        screen_false_positive_mentions(publication_text, mentions_by_gene)

    return [
        interest_rank_by_gene,
        mentions_by_gene,
        loci_by_gene,
        full_names_by_gene
    ]

def extract(organism, publication, bucket, de_dict):
    """Data mine genes from publication"""
    logger.info(f"Exracting gene names from publication: {publication}")

    [
        interest_rank_by_gene,
        mentions_by_gene,
        loci_by_gene,
        full_names_by_gene
    ] = extract_mentions_and_interest(organism, publication)

    de_by_gene = extract_de(bucket, de_dict)

    return [
        interest_rank_by_gene,
        mentions_by_gene,
        de_by_gene,
        loci_by_gene,
        full_names_by_gene
    ]

def sanitize_string(string):
    return string.replace(' ', '_').replace('.', '_').replace('[', '_').replace(']', '_')


def extract_de(bucket, de_dict):
    """Fetch differential expression (DE) files, return DE fields by gene"""

    de_by_gene = {}

    # directory = "tests/data/rank_genes/"
    directory = "_scp_internal/differential_expression/"

    for [group, de_file] in de_dict["de_groups_and_files"]:
        de_gsurl = f"gs://{bucket}/{directory}{de_file}"

        ingest_file = IngestFiles(de_gsurl)
        local_path = ingest_file.resolve_path(de_gsurl)[1]

        with open(local_path) as file:
            reader = csv.reader(file, delimiter="\t")
            i = 0
            # Headers in upstream, precomputed DE files:
            # names	scores	logfoldchanges	pvals	pvals_adj	pct_nz_group	pct_nz_reference
            for row in reader:
                if row[0] == "":
                    continue  # Skip novel header / empty lines
                if i < 500:
                    # 500 is ~2% of 25k-ish genes in Ensembl genome annotation
                    gene = row[1]
                    de_entry = {
                        "group": group,
                        "log2fc": row[3],  # Scanpy `logfoldchanges`
                        "adjusted_pval": row[5],  # Scanpy `pvals_adj`
                        "scores_rank": str(i),  # per Scanpy `scores`
                    }
                    if gene in de_by_gene:
                        de_by_gene[gene].append(de_entry)
                    else:
                        de_by_gene[gene] = [de_entry]
                i += 1

    sorted_de_by_gene = {}
    for gene in de_by_gene:
        de_entries = de_by_gene[gene]
        # Sort each group by Scanpy `scores` rank
        sorted_de = sorted(de_entries, key=lambda de: int(de["scores_rank"]))
        sorted_de_by_gene[gene] = sorted_de

    return sorted_de_by_gene


def get_de_column(de_entries, de_meta_keys):
    """Join DE entries for a gene by keys like log2fc, etc.
    """
    # Collapse DE props for each group, delimit inner fields with "!"
    de_grouped_props = []
    for de in de_entries:
        props = [de[key] for key in de_meta_keys]
        de_grouped_props.append("!".join(props))

    # Collapse grouped props, all DE data for each gene is in one TSV column
    de_column = ";".join(de_grouped_props)

    return de_column


def get_metainformation(meta, de_meta_keys):
    """Get headers about entire content, and inner column formats

    Metainformation fields are also used in the VCF file specification.
    """
    annot = meta["annotation"]
    meta["annotation"] = f"{annot['name']}--group--{annot['scope']}"
    content_meta = ";".join([f"{k}={v}" for (k, v) in meta.items()])
    de_meta = "differential_expression keys: " + ";".join(de_meta_keys)
    metainformation = (
        "## "
        + "\n## ".join(
            ["Ranked genes - Single Cell Portal", content_meta, de_meta]
        )
        + "\n"
    )

    return metainformation


def sort_genes_by_relevance(gene_dicts):
    """Sort genes by # mentions, then DE score rank, then global interest rank"""
    [
        interest_rank_by_gene,
        mentions_by_gene,
        de_by_gene,
        loci_by_gene,
        full_names_by_gene,
    ] = gene_dicts

    genes = []
    for gene in interest_rank_by_gene:
        de = de_by_gene.get(gene, '')
        top_de_rank = de[0]["scores_rank"] if de != '' else -1
        genes.append({
            "top_de_rank": top_de_rank,
            "de": de,
            "mentions": mentions_by_gene.get(gene, 0),
            "interest_rank": interest_rank_by_gene[gene],
            "locus": loci_by_gene[gene],
            "full_name": full_names_by_gene[gene],
            "symbol": gene,
        })

    genes = sorted(genes, key=lambda gene: int(gene["interest_rank"]))
    genes = sorted(genes, key=lambda gene: int(gene["top_de_rank"]))
    genes = sorted(genes, key=lambda gene: -int(gene["mentions"]))

    return genes


def transform(gene_dicts, meta):
    """Transform extracted dicts into TSV content"""

    num_genes_to_output = 100

    de_meta_keys = ["group", "log2fc", "adjusted_pval", "scores_rank"]

    ranked_genes = sort_genes_by_relevance(gene_dicts)

    rows = []
    i = 0
    for gene in ranked_genes:
        locus = gene["locus"]
        chromosome = locus["chromosome"]
        start = locus["start"]
        length = locus["length"]
        full_name = gene["full_name"]
        de = get_de_column(gene["de"], de_meta_keys)
        mentions = str(gene["mentions"])
        interest_rank = str(gene["interest_rank"])
        row = [
            gene["symbol"],
            chromosome,
            start,
            length,
            full_name,
            de,
            mentions,
            interest_rank,
        ]
        rows.append("\t".join(row))
        i += 1

    rows = "\n".join(rows[:num_genes_to_output])
    metainformation = get_metainformation(meta, de_meta_keys)
    header = (
        "# "
        + "\t".join([
                'name',
                'chromosome',
                'start',
                'length',
                'full_name',
                'differential_expression',
                'publication_mentions',
                'interest_rank',
            ])
        + "\n"
    )

    tsv_content = metainformation + header + rows

    return tsv_content


def load(clustering, annotation, tsv_content, bucket_name):
    """Load TSV content into file, write to disk"""

    # TODO: Consider cluster- and annotation-specific gene ranking
    # output_path = f"ranked_genes_{clustering}--{annotation}.tsv"
    output_path = "ranked_genes.tsv"
    with open(output_path, "w") as f:
        f.write(tsv_content)
    logger.info('Wrote content')
    print(tsv_content)

    bucket_folder = "_scp_internal/ranked_genes"
    gs_path = f"{bucket_folder}/{output_path}"
    output_gs_url = f"gs://{bucket_name}/{gs_path}"

    IngestFiles.delocalize_file(output_gs_url, output_path, gs_path)
    logger.info("Uploaded gene ranks to bucket")


def get_scp_api_origin():
    """Get domain etc. for SCP REST API URLs
    """
    db_name = os.environ['DATABASE_NAME']
    db_env = db_name.split('_')[-1]
    origins_by_environment = {
        'development': 'https://localhost:3000',
        'staging': 'https://singlecell-staging.broadinstitute.org',
        'production': 'https://singlecell.broadinstitute.org'
    }
    return origins_by_environment[db_env]


def fetch_context(accession):
    """Get cluster, annotation, and differential expression data from SCP API
    """
    origin = get_scp_api_origin()
    url = f"{origin}/single_cell/api/v1/studies/{accession}/explore"

    response = requests.get(url, verify=False)
    try:
        explore_json = response.json()
        if response.status_code == 401:
            print(
                '*** Received 401 error from SCP API.  ' +
                'Ensure study is public, not private. ***'
            )
    except json.decoder.JSONDecodeError as e:
        if 'staging' in origin:
            print(
                '*** Error requesting SCP API on staging.  ' +
                'Ensure you are connected to VPN. ***'
            )

        raise e

    bucket_id = explore_json["bucketId"]
    raw_organism = explore_json["taxonNames"][0]
    organism = raw_organism.lower().replace(' ', '-')

    try:
        # DE is available; use first eligible clustering and annotation
        # TODO: Consider outputting multiple gene rank lists, 1 per DE annot
        de = explore_json["differentialExpression"][0]
        clustering = de["cluster_name"]
        annotation = {
            "name": de["annotation_name"],
            "scope": de["annotation_scope"]
        }
        de_groups_and_files = de["select_options"]["one_vs_rest"]
        de_dict = {
            "clustering": clustering,
            "annotation": annotation,
            "de_groups_and_files": de_groups_and_files
        }
    except Exception as e:
        # No DE available; use default clustering and annotation
        # TODO: Add downstream handling for study without DE
        raise ValueError(
            "Rank genes only currently supports studies that have one-vs-rest SCP-computed DE"
        )

    return bucket_id, organism, de_dict


class RankGenes:

    def __init__(
        self,
        study_accession: str,
        publication: str
    ):
        """
        :param study_accession Accession of a public SCP study, e.g. "SCP123"
        :param publication URL of the study's publicly-accessible research article, or GS URL or local path to publication text file
        """

        bucket, organism, de_dict = fetch_context(study_accession)

        clustering = de_dict["clustering"]
        annotation = de_dict["annotation"]

        meta = {
            "accession": study_accession,
            "organism": organism,
            "bucket": bucket,
            "clustering": clustering,
            "annotation": annotation
        }

        gene_dicts = extract(organism, publication, bucket, de_dict)
        tsv_content = transform(gene_dicts, meta)

        load(clustering, annotation, tsv_content, bucket)
