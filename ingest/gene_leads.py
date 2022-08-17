"""Compute "gene leads", interesting genes to search in a single-cell study

Gene leads are all, at least, mentioned in a study's related publication.
These published genes are then sought among very differentially expressed
genes in the annotation's various labeled groups.  The groups in which each
published gene is very DE'd (if any) are noted.  Finally, data on each gene's
worldwide popularity, ultimately measured via the Gene Hints pipeline using
either Wikipedia page views or PubMed citation counts, is also included.

The output TSV file contains a one gene per row. This is parsed by Ideogram.js
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
    "neutrophils",
    "LC1",
    "LC2",
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
genes_url = f"https://github.com/eweitz/ideogram/blob/master/dist/data/cache/genes/{genes_filename}.gz?raw=true"
download_gzip(genes_url, genes_filename)

# Populate containers for various name and significance fields
view_rank_by_gene = {}
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
        view_rank_by_gene[gene] = i
        i += 1

publication_words = publication_text.split(' ')
for word in publication_words:
    raw_word = word.strip('*,.()—')
    if raw_word in counts_by_gene:
        # print(raw_word)
        counts_by_gene[raw_word] += 1

seen_genes = {}
for gene in counts_by_gene:
    counts_by_gene[gene]
    count = counts_by_gene[gene]
    if counts_by_gene[gene] > 0:
        seen_genes[gene] = count

very_de_by_gene = {}

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
                # TODO (pre-GA):
                # - Account for studies with few assayed genes
                # - Include log2FC, pvals_adj for each very-DE group
                if i < 500:
                    # 500 is ~2% of 25k-ish genes in Ensembl genome annotation.
                    # Like the rest of "gene leads", a big factor in this value
                    # choice is its ability to yield an engaging UI.  In this case,
                    # 500 produced a wide-but-not-overwhelming palette of colors
                    # in the gene leads ideogram.
                    gene = row[1]
                    if gene in very_de_by_gene:
                        very_de_by_gene[gene].append(group)
                    else:
                        very_de_by_gene[gene] = [group]
                i += 1

process_de_files()

br2 = "<br/><br/>"
def contextualize_count(count):
    times = f" {count} times" if count > 1 else ""
    return f"Mentioned{times} in linked publication."

def contextualize_rank(rank):
    """Note gene's global popularity, if in top ~2% for organism's known genes
    """
    text = f"Ranked {rank} in global interest." if rank < 500 else ""
    return text

def contextualize_de(gene, very_de_by_gene):
    # TODO (pre-GA):
    # - Include log2FC, pvals_adj for each high-DE group
    # - Improve formatting for cases where gene is high-DE in many groups
    if gene not in very_de_by_gene:
        return ["", "#4d72aa"] # Empty string, default color
    groups = very_de_by_gene[gene]
    pretty_groups = " and ".join(groups)
    text = f"Very differentially expressed in:<br/>{pretty_groups}{br2}"
    first_group_index = all_groups.index(groups[0])
    color = color_brewer_list[first_group_index]
    return [text, color]

# Process main content for output
rows = []
i = 0
for gene in seen_genes:
    loci = loci_by_gene[gene]
    chromosome = loci["chromosome"]
    start = loci["start"]
    length = loci["length"]
    full_name = full_names_by_gene[gene]

    # TODO (pre-GA): Handling of line breaks, prose wording, color, etc. should
    # be handled in SCP UI.  Expose Ideogram kit internals via its JS API to
    # enable that.
    [de_text, color] = contextualize_de(gene, very_de_by_gene)
    count_text = contextualize_count(seen_genes[gene])
    rank_text = contextualize_rank(view_rank_by_gene[gene])
    if rank_text != "":
        count_text += br2
    significance = de_text + count_text + rank_text
    row = [gene, chromosome, start, length, color, full_name, significance]
    rows.append("\t".join(row))
    i += 1

rows = "\n".join(rows)
header = "##name\tchromosome\tstart\tlength\tcolor\tfull_name\tsignificance\n"

content = header + rows

cache_buster = "_v2" # TODO (pre-GA): Improve handling if needed beyond dev
output_path = f"gene_leads_{clustering}--{annotation}{cache_buster}.tsv"
with open(output_path, "w") as f:
    f.write(content)

print('Wrote content')
print(content)

# TODO (pre-GA): Write output files to whichever bucket we write DE to.
# The # of gene leads files will be many fewer than DE in Q2 '22,
# i.e. << A-vs-all DE.
