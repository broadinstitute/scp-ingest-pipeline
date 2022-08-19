"""Compute "gene leads", interesting genes to search in a single-cell study

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
genes_url = f"https://cdn.jsdelivr.net/npm/ideogram@1.37.0/dist/data/cache/{genes_filename}.gz"
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

mentioned_genes = {}
for gene in counts_by_gene:
    counts_by_gene[gene]
    count = counts_by_gene[gene]
    if counts_by_gene[gene] > 0:
        mentioned_genes[gene] = count

de_by_gene = {}

def process_de_files():
    """Fetch, write, and process precomputed differential expression files
    """
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
            # Headers:
            # names	scores	logfoldchanges	pvals	pvals_adj	pct_nz_group	pct_nz_reference
            for row in reader:
                if i < 500:
                    # 500 is ~2% of 25k-ish genes in Ensembl genome annotation.
                    # Like the rest of "gene leads", a big factor in this value
                    # choice is its ability to yield an engaging UI.  In this case,
                    # 500 produced a wide-but-not-overwhelming palette of colors
                    # in the gene leads ideogram.
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

process_de_files()

def get_de_and_color_columns(gene, de_by_gene):
    # TODO (SCP-4061): Move color handling to SCP UI
    if gene not in de_by_gene:
        return ["", "#4d72aa"] # Empty string, default color
    de_entries = de_by_gene[gene]
    # Sort each group by Scanpy `scores` rank
    de_entries = sorted(de_entries, key=lambda de: int(de["scores_rank"]))

    # Collapse DE props for each group, delimit inner fields with "!"
    de_grouped_props = []
    for de in de_entries:
        de_grouped_props.append(
            "!".join([
                de["group"], de["log2fc"],
                de["adjusted_pval"], de["scores_rank"]
            ])
        )

    # Collapse grouped props, all DE data for each gene is in one TSV column
    de_column = ";".join(de_grouped_props)

    # Set color to that of the highest scores-ranked group
    # TODO (SCP-4061): Move color handling to SCP UI
    first_group_index = all_groups.index(de_entries[0]["group"])
    color = color_brewer_list[first_group_index]

    return [de_column, color]

# Process main content for output
rows = []
i = 0
for gene in mentioned_genes:
    loci = loci_by_gene[gene]
    chromosome = loci["chromosome"]
    start = loci["start"]
    length = loci["length"]
    full_name = full_names_by_gene[gene]

    [de, color] = get_de_and_color_columns(gene, de_by_gene)
    mentions = str(mentioned_genes[gene])
    interest_rank = str(interest_rank_by_gene[gene])
    row = [
        gene, chromosome, start, length, color, full_name,
        de, mentions, interest_rank
    ]
    rows.append("\t".join(row))
    i += 1

rows = "\n".join(rows)
header = "##" + "\t".join([
    'name', 'chromosome', 'start', 'length', 'color', 'full_name',
    'differential_expression', 'publication_mentions', 'interest_rank'
]) + "\n"

content = header + rows

cache_buster = "_v3" # TODO (pre-GA): Improve handling if needed beyond dev
output_path = f"gene_leads_{clustering}--{annotation}{cache_buster}.tsv"
with open(output_path, "w") as f:
    f.write(content)

print('Wrote content')
print(content)

# TODO (pre-GA): Write output files to whichever bucket we write DE to.
# The # of gene leads files will be many fewer than DE in Q2 '22,
# i.e. << A-vs-all DE.
