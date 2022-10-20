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
"""

import gzip
import urllib
import urllib.request as request
import csv

# TODO (pre-GA): Extract these to CLI arguments and/or SCP API calls
organism = "homo-sapiens"
bucket = "fc-65379b91-5ded-4d28-8e51-ada209542117"
accession = "SCP138"
clustering = "All_Cells_UMAP"
annotation = "General_Celltype"
all_groups = [
    "B cells",
    "CSN1S1 macrophages",
    "dendritic cells",
    "eosinophils",
    "fibroblasts",
    "GPMNB macrophages",
    "LC1",
    "LC2",
    "neutrophils",
    "T cells"
]

# Temporarily duplicated from SCP UI code.
#
# TODO (pre-GA): Expose related genes kit internals via Ideogram.js
# so the end client UI (i.e., SCP UI) can handle color, etc.
color_brewer_list = [
  '#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#a65628',
  '#f781bf', '#999999', '#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3',
  '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3', '#8dd3c7', '#bebada',
  '#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#d9d9d9',
  '#bc80bd', '#ccebc5', '#ffed6f'
]

def download_gzip(url, output_path, cache=0):
    """Download gzip file, decompress, write to output path; use optional cache
    Cached files can help speed development iterations by > 2x, and some
    development scenarios (e.g. on a train or otherwise without an Internet
    connection) can be impossible without it.
    """
    headers={"Accept-Encoding": "gzip"}
    print(f"Requesting {url}")
    request = urllib.request.Request(url, headers=headers)
    response = urllib.request.urlopen(request)
    content = gzip.decompress(response.read()).decode()

    with open(output_path, "w") as f:
        f.write(content)

def extract(meta, all_groups):
    [accession, organism, bucket, clustering, annotation] = list(meta.values())
    # Example URL:
    # https://www.biorxiv.org/content/10.1101/2021.11.13.468496v1.full.txt
    #
    # TODO (pre-GA, story): Resolve publication IDs to text.   This would take in
    # any URL to a public (not-needing-authorization) publication in any major
    # venue (e.g. Nature, Science, bioRxiv, etc.) and resolve it a URL that will
    # return an easily machine-parseable (e.g. plaintext) file containing the
    # publication's entire prose content (excluding e.g. supplementary files).
    #
    # If it fails, log an error to Sentry and/or Mixpanel and exit PAPI
    url = f'https://www.biorxiv.org/content/10.1101/2021.11.13.468496v1.full.txt'
    with request.urlopen(url) as response:
        publication_text = response.read().decode('utf-8')

    genes_filename = f"{organism}-genes.tsv" # e.g. homo-sapiens-genes.tsv
    base_url = "https://cdn.jsdelivr.net/npm/"
    genes_url = f"{base_url}ideogram@1.37.0/dist/data/cache/{genes_filename}.gz"
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
            if row[0] == '#' or len(row) < 2: continue
            gene = row[4]
            full_name = row[5]
            counts_by_gene[gene] = 0
            loci_by_gene[gene] = {
                "chromosome": row[0],
                "start": row[1],
                "length": row[2]
            }
            full_names_by_gene[gene] = full_name
            # Genes in upstream file are ordered by global popularity
            interest_rank_by_gene[gene] = i
            i += 1

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
        full_names_by_gene
    ]

def extract_de(bucket, clustering, annotation, all_groups):
    """Fetch differential expression (DE) files, return DE fields by gene
    """
    de_by_gene = {}

    origin = "https://storage.googleapis.com"
    directory = "_scp_internal%2Fdifferential_expression%2F"
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
                if row[0] == "": continue # Skip novel header / empty lines
                if i < 500:
                    # 500 is ~2% of 25k-ish genes in Ensembl genome annotation.
                    # Like the rest of "gene leads", a big factor in this value
                    # choice is its ability to yield an engaging UI.  In this
                    # case, 500 produced a wide-but-not-overwhelming palette
                    # of annotations in the gene leads ideogram.
                    gene = row[1]
                    de_entry = {
                        "group": group,
                        "log2fc": row[3], # Scanpy `logfoldchanges`
                        "adjusted_pval": row[5], # Scanpy `pvals_adj`
                        "scores_rank": str(i) # per Scanpy `scores`
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
        de_grouped_props.append(
            "!".join(props)
        )

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
    metainformation = "## " + "\n## ".join([
        "Gene leads ideogram data - Single Cell Portal",
        content_meta,
        de_meta
    ]) + "\n"

    return metainformation

def sort_genes_by_relevance(gene_dicts):
    """ Sort genes by DE score rank, then # mentions, then global interest rank
    """
    [
        interest_rank_by_gene,
        mentions_by_gene,
        de_by_gene,
        loci_by_gene,
        full_names_by_gene
    ] = gene_dicts

    # print('de_by_gene')
    # print(de_by_gene)

    genes = []
    for gene in de_by_gene:
        if gene not in interest_rank_by_gene:
            # TODO: Handle synonyms, e.g. CECR1 for ADA2
            continue
        genes.append({
            "top_de_rank": de_by_gene[gene][0]["scores_rank"],
            "de": de_by_gene[gene],
            "mentions": mentions_by_gene.get(gene, 0),
            "interest_rank": interest_rank_by_gene[gene],
            "locus": loci_by_gene[gene],
            "full_name": full_names_by_gene[gene],
            "symbol": gene
        })

    genes = sorted(genes, key=lambda gene: int(gene["interest_rank"]))
    genes = sorted(genes, key=lambda gene: -int(gene["mentions"]))
    genes = sorted(genes, key=lambda gene: int(gene["top_de_rank"]))

    return genes

def transform(gene_dicts, meta):
    """Transform extracted dicts into TSV content
    """

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
            gene["symbol"], chromosome, start, length, color, full_name,
            de, mentions, interest_rank
        ]
        rows.append("\t".join(row))
        i += 1

    rows = "\n".join(rows[:30])
    metainformation = get_metainformation(meta, de_meta_keys)
    header = "# " + "\t".join([
        'name', 'chromosome', 'start', 'length', 'color', 'full_name',
        'differential_expression', 'publication_mentions', 'interest_rank'
    ]) + "\n"

    tsv_content = metainformation + header + rows

    return tsv_content

def load(tsv_content):
    """Load TSV content into file, write to disk
    """
    # TODO (pre-GA): Write output files to whichever bucket we write DE to.
    # The # of gene leads files will be many fewer than DE in Q2 '22,
    # i.e. << A-vs-all DE.
    cache_buster = "_v5" # TODO (pre-GA): Improve handling if needed beyond dev
    output_path = f"gene_leads_{clustering}--{annotation}{cache_buster}.tsv"
    with open(output_path, "w") as f:
        f.write(tsv_content)

    print('Wrote content')
    print(tsv_content)

def main(organism, bucket, clustering, annotation, all_groups):
    """Extract, transform, and load data into TSV files for gene leads ideogram
    """
    meta = {
        "accession": accession,
        "organism": organism,
        "bucket": bucket,
        "clustering": clustering,
        "annotation": annotation
    }
    gene_dicts = extract(meta, all_groups)
    tsv_content = transform(gene_dicts, meta)
    load(tsv_content)

if __name__ == '__main__':
    main(organism, bucket, clustering, annotation, all_groups)
