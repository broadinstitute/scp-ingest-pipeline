import gzip
import urllib
import urllib.request as request
import csv

organism = "homo-sapiens"
bucket = "fc-65379b91-5ded-4d28-8e51-ada209542117"
clustering = "All_Cells_UMAP"
annotation = "--General_Celltype--"
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
# https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/assembly_summary.txt
url = f'https://www.biorxiv.org/content/10.1101/2021.11.13.468496v1.full.txt'
with request.urlopen(url) as response:
    publication_text = response.read().decode('utf-8')

genes_url = f"https://github.com/eweitz/ideogram/blob/master/dist/data/cache/genes/{organism}-genes.tsv.gz?raw=true"
download_gzip(genes_url, 'genes.txt')

view_rank_by_gene = {}
counts_by_gene = {}
loci_by_gene = {}
full_names_by_gene = {}
i = 0
with open('genes.txt') as file:
    reader = csv.reader(file, delimiter="\t")
    for row in reader:
        if row[0] == '#' or len(row) < 2: continue
        gene = row[4]
        full_name = row[5]
        counts_by_gene[gene] = 0
        view_rank_by_gene[gene] = i
        loci_by_gene[gene] = {
            "chromosome": row[0],
            "start": row[1],
            "length": row[2]
        }
        full_names_by_gene[gene] = full_name
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

most_de_by_gene = {}

def write_de_files():
    origin = "https://storage.googleapis.com"
    directory = "_scp_internal%2Fdifferential_expression%2F"
    de_url_stem = f"{origin}/download/storage/v1/b/{bucket}/o/{directory}"
    leaf = "--study--wilcoxon.tsv"
    params = "?alt=media"

    for group in all_groups:
        safe_group = group.replace(' ', '_')
        de_filename = "_scp_internal_differential_expression_" + clustering + annotation + safe_group + leaf

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
                    gene = row[1]
                    if gene in most_de_by_gene:
                        most_de_by_gene[gene].append(group)
                    else:
                        most_de_by_gene[gene] = [group]
                i += 1

write_de_files()

i = 0
with open('genes.txt') as file:
    reader = csv.reader(file, delimiter="\t")
    for row in reader:
        if row[0] == '#' or len(row) < 2: continue
        gene = row[4]
        full_name = row[5]
        counts_by_gene[gene] = 0
        view_rank_by_gene[gene] = i
        loci_by_gene[gene] = {
            "chromosome": row[0],
            "start": row[1],
            "length": row[2]
        }
        full_names_by_gene[gene] = full_name
        i += 1


# sorted_seen_genes = dict(sorted(seen_genes.items(), key=lambda item: -item[1]))

def contextualize_count(count):
    # return f"Mentioned {count} times in linked publication."
    return f"Mentioned in linked publication."

def contextualize_rank(rank):
    return f"Ranked {rank} in worldwide interest."

def contextualize_de(groups):
    pretty_groups = " and ".join(groups)
    text = f"Highly differentially expressed in:<br/>{pretty_groups}"
    first_group_index = all_groups.index(groups[0])
    color = color_brewer_list[first_group_index]
    return [text, color]

print("most_de_by_gene")
print(most_de_by_gene)
rows = []
i = 0
for gene in seen_genes:
    loci = loci_by_gene[gene]
    chromosome = loci["chromosome"]
    start = loci["start"]
    length = loci["length"]
    full_name = full_names_by_gene[gene]
    count_text = contextualize_count(seen_genes[gene])
    rank_text = contextualize_rank(view_rank_by_gene[gene])
    de_text = ""
    color = "#4d72aa"
    if gene in most_de_by_gene:
        de_context = contextualize_de(most_de_by_gene[gene])
        de_text = "<br/><br/>" + de_context[0]
        color = de_context[1]
    significance = f"{count_text}<br/>{rank_text}{de_text}"
    row = [gene, chromosome, start, length, color, full_name, significance]
    rows.append("\t".join(row))
    i += 1

rows = "\n".join(rows)
header = "##name\tchromosome\tstart\tlength\tcolor\tfull_name\tsignificance\n"

content = header + rows
with open("gene_highlights_v4.tsv", "w") as f:
    f.write(content)

print('Wrote content')
print(content)

# print(sorted_seen_genes)

