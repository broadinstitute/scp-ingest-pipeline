#using https://nbviewer.org/github/theislab/diffxpy_tutorials/blob/master/diffxpy_tutorials/test/multiple_tests_per_gene.ipynb
import pandas as pd
import scanpy as sc
import diffxpy.api as de



# TODO (SCP-5041): Extract to CLI arguments
#Pairwise tests between groups
#answers whether a given pair of groups shows differential expression for each gene

#arguments:
#data: Anndata object, data matrix with cells x genes

#grouping: str, column in data, the column that contains cell type labels. alternatively vector containing group labels

#test: type of test, default is z-test. other options: ’wald’, ’lrt’, ’t-test’, ’rank’

#lazy: bool, only possible if test is ztest, if true only evaluated once the user requests pval/coefficients for a specific pair of models

#noise_model: default is nb, specify NONE for wald and t test

h5ad_file = "ingest/B_Plasma.h5ad"
adata_file = sc.read_h5ad(h5ad_file)
test = de.test.pairwise(
    data=adata_file,
    grouping="ClusterName",
    #"rank" test correlates to Wilcoxon rank-sum
    test="rank",
    lazy=False,
    noise_model=None
 )

#all built in results returned from test function
cell_type = test.groups
pval = test.pval
lfc = test._logfc
names = test.gene_ids
mean = test.mean
qval = test.qval


# TODO change from hardcoding to parsing these arguments 

cluster_name = "umap_coords"
annotation_name = "ClusterName"
annot_scope = "B_cell_cluster" # TODO verify this is accurate for files
method = "wilcoxon"


#unfold p values array to go through in main loop
pval_clusters = {}
for i in range (len(pval)):
    pval_clusters[cell_type[i]] = pval[i]

lfc_clusters = {}
for i in range (len(lfc)):
    lfc_clusters[cell_type[i]] = lfc[i]

qval_clusters = {}
for i in range (len(qval)):
    qval_clusters[cell_type[i]] = qval[i]

inner_dict = {}

for j, i, q in zip(pval_clusters, lfc_clusters, qval_clusters):

    inner_pval = pval_clusters[j]
    inner_lfc = lfc_clusters[i]
    inner_qval = qval_clusters[i]

    count_inner = 0

    for k, p, o in zip(inner_pval, inner_lfc, inner_qval) :
        current_name = cell_type[count_inner]
        count_inner +=1
        count = 0
        inner_dict[j, current_name] = []
        #16 arrays, 4 for each
        for l, t, h in zip(k, p, o):
            if (j != current_name):
                #order: gene name, pvalue, qvalue, lfc
                current_list = [names[count], l, h, t]
                inner_dict[j, current_name].append(current_list)
                count += 1

col_vals = ["names", "pvalues", "qvalues", "log2FC"]
for i in inner_dict:
    arr = inner_dict[i]
    if len(arr) != 0:
        inner_df = pd.DataFrame(data = arr, columns = col_vals)
        #with open("ingest/{}--{}--{}--{}--{}--{}.tsv".format(cluster_name, annotation_name, i[0], i[1], annot_scope, method), "w") as external_file:
        #    print(inner_df.to_string(), file = external_file)
        inner_df.to_csv("ingest/{}--{}--{}--{}--{}--{}.tsv".format(cluster_name, annotation_name, i[0], i[1], annot_scope, method), sep ='\t')

'''
structure:
all array [B cell array, Plasma cell array, Plasma cell (long-lived; higher BLIMP1) array, Plasmablasts array]
B cell array [[B cells + B cells], [B cells + Plasma cells], [B cells + Plasma cell (long-lived; higher BLIMP1)], [B cells + Plasma cell other], [B cells + Plasmablasts]]
^each interior array is 30983

so to cycle through: first go through all array, then inner loop of each one

B cell array
    measurments:
        B cells + B cells
            for 30983 gene names
        B cells + Plasma cells
            for 30983 gene names
        B cells + Plasma cell (long-lived; higher BLIMP1)
            for 30983 gene names
        B cells + Plasmablasts
            for 30983 gene names

total = 495728
total without self loops = 371796


and so on for all others

Plasma cell array

Plasma cell (long-lived; higher BLIMP1) array

Plasmablasts array

'''

#CODE FOR RESULTS ALL IN ONE FILE
""" 
dict_results["cell_type"] = []
dict_results["compared_cell_type"] = []
dict_results["pvals"] = []
dict_results["log2FC"] = []

pval_clusters = {}
for i in range (len(pval)):
    pval_clusters[cell_type[i]] = pval[i]

lfc_clusters = {}
for i in range (len(lfc)):
    lfc_clusters[cell_type[i]] = lfc[i]

#unfold log fold change value array
lfc_list = []
for i in lfc:
    for j in i:
        for l in j:
            lfc_list.append(l)

total_loop = 0
for j, i in zip(pval_clusters, lfc_clusters):
    inner = pval_clusters[j]
    inner_lfc = lfc_clusters[i]
    count_inner = 0
    for k, p in zip(inner, inner_lfc) :
        current_name = cell_type[count_inner]
        count_inner +=1
        count = 0
        for l, t in zip(k, p):
            if (j != current_name):
                #print(j + "not equal to" + current_name)
                #total_loop +=1
                dict_results["names"].append(names[count])
                dict_results["cell_type"].append(j)
                dict_results["compared_cell_type"].append(current_name)
                dict_results["pvals"].append(l)
                dict_results["log2FC"].append(t)
                count += 1

#dict_results["log2FC"] = lfc_list
new_df = pd.DataFrame.from_dict(dict_results)

"""
