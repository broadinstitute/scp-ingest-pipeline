import gzip
import urllib
import urllib.request as request
import csv

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

genes_url = 'https://github.com/eweitz/ideogram/blob/master/dist/data/cache/genes/homo-sapiens-genes.tsv.gz?raw=true'
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
    de_url_stem = "https://storage.googleapis.com/download/storage/v1/b/fc-65379b91-5ded-4d28-8e51-ada209542117/o/_scp_internal%2Fdifferential_expression%2F"
    clustering = "All_Cells_UMAP"
    annotation = "--General_Celltype--"
    groups = [
        "B cells",
        "LC2",
        "GPMNB macrophages",
        "neutrophils",
        "T cells",
        "CSN1S1 macrophages",
        "dendritic cells",
        "LC1",
        "eosinophils",
        "fibroblasts"
    ]
    leaf = "--study--wilcoxon.tsv"
    params = "?alt=media"

    for group in groups:
        de_filename = "_scp_internal_differential_expression_" + clustering + annotation + group.replace(' ', '_') + leaf

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
                if i < 10:
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
    return f"Highly differentially expressed in:<br/>{pretty_groups}"

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
    count_context = contextualize_count(seen_genes[gene])
    rank_context = contextualize_rank(view_rank_by_gene[gene])
    de_context = ""
    if gene in most_de_by_gene:
        de_context = "<br/><br/>" + contextualize_de(most_de_by_gene[gene])
    significance = f"{count_context}<br/>{rank_context}{de_context}"
    row = [gene, chromosome, start, length, "#4d72aa", full_name, significance]
    rows.append("\t".join(row))
    i += 1

rows = "\n".join(rows)
header = "##name\tchromosome\tstart\tlength\tcolor\tfull_name\tsignificance\n"

content = header + rows
with open("gene_highlights.tsv", "w") as f:
    f.write(content)

print('Wrote content')
print(content)

# print(sorted_seen_genes)

