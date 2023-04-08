"""Extract, transform, and load data into TSV files for gene leads ideogram

Gene leads are all, at least, mentioned in a study's related publication.
These published genes are then sought among very differentially expressed
genes in the annotation's various labeled groups.  The groups in which each
published gene is very DE'd (if any) are noted.  Finally, data on each gene's
worldwide popularity, ultimately measured via the Gene Hints pipeline using
either Wikipedia page views or PubMed citation counts, is also included.

The output TSV file contains one gene per row. This is parsed by Ideogram.js
in the SCP UI to display gene leads.  The genome visualization is designed to
engage users with relevant information, and prime the study gene search UX.

EXAMPLES

python3 ingest/gene_leads.py --study-accession SCP138 --bucket-name fc-65379b91-5ded-4d28-8e51-ada209542117 --taxon-name="Homo sapiens" --clustering="All Cells UMAP" --annotation-name="General Celltype" --annotation-groups '["B cells", "CSN1S1 macrophages", "dendritic cells", "eosinophils", "fibroblasts", "GPMNB macrophages", "LC1", "LC2", "neutrophils", "T cells"]' --publication-url https://www.biorxiv.org/content/10.1101/2021.11.13.468496v1
"""

import argparse
import ast
import glob
import gzip
import tarfile
import json
import urllib
import urllib.request
import csv
from io import BytesIO
import xml.etree.ElementTree as ET

# Temporarily duplicated from SCP UI code.
#
# TODO (pre-GA): Expose related genes kit internals via Ideogram.js
# so the end client UI (i.e., SCP UI) can handle color, etc.
color_brewer_list = [
    '#e41a1c',
    '#377eb8',
    '#4daf4a',
    '#984ea3',
    '#ff7f00',
    '#a65628',
    '#f781bf',
    '#999999',
    '#66c2a5',
    '#fc8d62',
    '#8da0cb',
    '#e78ac3',
    '#a6d854',
    '#ffd92f',
    '#e5c494',
    '#b3b3b3',
    '#8dd3c7',
    '#bebada',
    '#fb8072',
    '#80b1d3',
    '#fdb462',
    '#b3de69',
    '#fccde5',
    '#d9d9d9',
    '#bc80bd',
    '#ccebc5',
    '#ffed6f',
]


def download_gzip(url, output_path, cache=0):
    """Download gzip file, decompress, write to output path; use optional cache
    Cached files can help speed development iterations by > 2x, and some
    development scenarios (e.g. on a train or otherwise without an Internet
    connection) can be impossible without it.
    """
    is_targz = url.endswith("tar.gz")
    headers = {"Accept-Encoding": "gzip"} if is_targz else {}
    print(f"Requesting {url}")
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
    """Convert Digital Object Identifier (DOI) into PubMed Central ID (PMCID)"""
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
    """Get full text for publication given its in PMC ID"""
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
    text = ET.tostring(article.find("body"), method="text")
    return text


def fetch_publication_text(publication_url):
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
    doi_split = publication_url.split("https://doi.org/")
    if len(doi_split) > 1:
        doi = doi_split[1]
        pmcid = fetch_pmcid(doi)
        text = fetch_pmcid_text(pmcid)
    elif publication_url.startswith("https://www.biorxiv.org"):
        full_text_url = f"{publication_url}.full.txt"
        with urllib.request.urlopen(full_text_url) as response:
            text = response.read().decode('utf-8')
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


def extract(meta, all_groups, publication_url):
    [accssion, organism, bucket, clustering, annotation] = list(meta.values())
    [
        interest_rank_by_gene,
        counts_by_gene,
        loci_by_gene,
        full_names_by_gene,
    ] = fetch_gene_cache(organism)

    publication_text = fetch_publication_text(publication_url)

    publication_words = publication_text.split(' ')
    for word in publication_words:
        raw_word = word.strip('*,.()—')
        if raw_word in counts_by_gene:
            # print(raw_word)
            counts_by_gene[raw_word] += 1

    mentions_by_gene = {}
    for gene in counts_by_gene:
        counts_by_gene[gene]
        count = counts_by_gene[gene]
        if counts_by_gene[gene] > 0:
            mentions_by_gene[gene] = count

    de_by_gene = extract_de(bucket, clustering, annotation, all_groups)

    return [
        interest_rank_by_gene,
        mentions_by_gene,
        de_by_gene,
        loci_by_gene,
        full_names_by_gene,
    ]


def extract_de(bucket, clustering, annotation, all_groups):
    """Fetch differential expression (DE) files, return DE fields by gene"""
    de_by_gene = {}

    origin = "https://storage.googleapis.com"
    directory = "tests/data/gene_leads/_scp_internal%2Fdifferential_expression%2F"
    de_url_stem = f"{origin}/download/storage/v1/b/{bucket}/o/{directory}"
    leaf = "--study--wilcoxon.tsv"
    params = "?alt=media"

    for group in all_groups:
        safe_group = group.replace(' ', '_')
        tmp_dir = directory.replace('%2F', '_')
        de_filename = f"{tmp_dir}{clustering}--{annotation}--{safe_group}{leaf}"

        # TODO (pre-GA): Fetch these from bucket; requires auth token
        # de_url = de_url_stem + de_filename + params
        # with request.urlopen(de_url) as response:
        #     de_content = response.read().decode('utf-8')
        # with open(de_filename, "w") as f:
        #     f.write(de_content)

        with open(de_filename) as file:
            reader = csv.reader(file, delimiter="\t")
            i = 0
            # Headers in upstream, precomputed DE files:
            # names	scores	logfoldchanges	pvals	pvals_adj	pct_nz_group	pct_nz_reference
            for row in reader:
                if row[0] == "":
                    continue  # Skip novel header / empty lines
                if i < 500:
                    # 500 is ~2% of 25k-ish genes in Ensembl genome annotation.
                    # Like the rest of "gene leads", a big factor in this value
                    # choice is its ability to yield an engaging UI.  In this
                    # case, 500 produced a wide-but-not-overwhelming palette
                    # of annotations in the gene leads ideogram.
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


def get_de_and_color_columns(de_entries, de_meta_keys):
    # Collapse DE props for each group, delimit inner fields with "!"
    de_grouped_props = []
    for de in de_entries:
        props = [de[key] for key in de_meta_keys]
        de_grouped_props.append("!".join(props))

    # Collapse grouped props, all DE data for each gene is in one TSV column
    de_column = ";".join(de_grouped_props)

    # Set color to that of the highest scores-ranked group
    # TODO (SCP-4061): Move color handling to SCP UI
    first_group_index = all_groups.index(de_entries[0]["group"])
    color = color_brewer_list[first_group_index]

    return [de_column, color]


def get_metainformation(meta, de_meta_keys):
    """Get headers about entire content, and inner column formats

    Metainformation fields are also used in the VCF file specification.
    """
    content_meta = ";".join([f"{k}={v}" for (k, v) in meta.items()])
    de_meta = "differential_expression keys: " + ";".join(de_meta_keys)
    metainformation = (
        "## "
        + "\n## ".join(
            ["Gene leads ideogram data - Single Cell Portal", content_meta, de_meta]
        )
        + "\n"
    )

    return metainformation


def sort_genes_by_relevance(gene_dicts):
    """Sort genes by DE score rank, then # mentions, then global interest rank"""
    [
        interest_rank_by_gene,
        mentions_by_gene,
        de_by_gene,
        loci_by_gene,
        full_names_by_gene,
    ] = gene_dicts

    # print('de_by_gene')
    # print(de_by_gene)

    genes = []
    for gene in de_by_gene:
        if gene not in interest_rank_by_gene:
            # TODO: Handle synonyms, e.g. CECR1 for ADA2
            continue
        genes.append(
            {
                "top_de_rank": de_by_gene[gene][0]["scores_rank"],
                "de": de_by_gene[gene],
                "mentions": mentions_by_gene.get(gene, 0),
                "interest_rank": interest_rank_by_gene[gene],
                "locus": loci_by_gene[gene],
                "full_name": full_names_by_gene[gene],
                "symbol": gene,
            }
        )

    genes = sorted(genes, key=lambda gene: int(gene["interest_rank"]))
    genes = sorted(genes, key=lambda gene: -int(gene["mentions"]))
    genes = sorted(genes, key=lambda gene: int(gene["top_de_rank"]))

    return genes


def transform(gene_dicts, meta):
    """Transform extracted dicts into TSV content"""

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
        [de, color] = get_de_and_color_columns(gene["de"], de_meta_keys)
        mentions = str(gene["mentions"])
        interest_rank = str(gene["interest_rank"])
        row = [
            gene["symbol"],
            chromosome,
            start,
            length,
            color,
            full_name,
            de,
            mentions,
            interest_rank,
        ]
        rows.append("\t".join(row))
        i += 1

    rows = "\n".join(rows[:30])
    metainformation = get_metainformation(meta, de_meta_keys)
    header = (
        "# "
        + "\t".join(
            [
                'name',
                'chromosome',
                'start',
                'length',
                'color',
                'full_name',
                'differential_expression',
                'publication_mentions',
                'interest_rank',
            ]
        )
        + "\n"
    )

    tsv_content = metainformation + header + rows

    return tsv_content


def load(tsv_content):
    """Load TSV content into file, write to disk"""
    # TODO (pre-GA): Write output files to whichever bucket we write DE to.
    # The # of gene leads files will be many fewer than DE in Q2 '22,
    # i.e. << A-vs-all DE.
    cache_buster = "_v5"  # TODO (pre-GA): Improve handling if needed beyond dev
    output_path = f"gene_leads_{clustering}--{annotation}{cache_buster}.tsv"
    with open(output_path, "w") as f:
        f.write(tsv_content)

    print('Wrote content')
    print(tsv_content)


def gene_leads(
    accession, bucket, organism, clustering, annotation, all_groups, publication_url
):
    """Extract, transform, and load data into TSV files for gene leads ideogram"""
    organism = organism.lower().replace(' ', '-')
    clustering = clustering.replace(' ', '_')
    annotation = annotation.replace(' ', '_')
    meta = {
        "accession": accession,
        "organism": organism,
        "bucket": bucket,
        "clustering": clustering,
        "annotation": annotation,
    }
    gene_dicts = extract(meta, all_groups, publication_url)
    tsv_content = transform(gene_dicts, meta)
    load(tsv_content)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "--study-accession",
        required=True,
        help="Single study accession associated with provided DE input files.",
    )
    parser.add_argument(
        "--bucket-name",
        required=True,
        help=("Name of GCS bucket, e.g. fc-65379b91-5ded-4d28-8e51-ada209541234"),
    )
    parser.add_argument(
        "--taxon-name",
        required=True,
        help=("Scientific name of organism, e.g. \"Homo sapiens\""),
    )
    parser.add_argument("--clustering", required=True, help=("Name of clustering"))
    parser.add_argument("--annotation-name", required=True, help=("Name of annotation"))
    parser.add_argument(
        "--annotation-groups",
        required=True,
        type=ast.literal_eval,
        help=("List of annotation groups, e.g. ['B cells', 'CSN1S1 macrophages']"),
    )
    parser.add_argument(
        "--publication-url",
        required=True,
        help=("URL of the study's publicly-accessible research article"),
    )

    args = parser.parse_args()

    accession = args.study_accession
    bucket = args.bucket_name
    organism = args.taxon_name
    clustering = args.clustering
    annotation = args.annotation_name
    all_groups = args.annotation_groups
    publication_url = args.publication_url

    gene_leads(
        accession, bucket, organism, clustering, annotation, all_groups, publication_url
    )
