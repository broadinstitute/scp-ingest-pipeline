#using https://nbviewer.org/github/theislab/diffxpy_tutorials/blob/master/diffxpy_tutorials/test/multiple_tests_per_gene.ipynb
import anndata
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import numpy as np
import pandas as pd
import scipy.stats
import scanpy as sc
import os
import batchglm.api as glm
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
dict_results = {}

#column names in the final synthetic file results
dict_results["names"] = []
dict_results["cell_type"] = []
dict_results["compared_cell_type"] = []
dict_results["pvals"] = []

#unfold p values array to go through in main loop
pval_clusters = {}
for i in range (len(pval)):
    pval_clusters[cell_type[i]] = pval[i]

#unfold log fold change value array
lfc_list = []
for i in lfc:
    for j in i:
        for l in j:
            lfc_list.append(l)

for j in pval_clusters:
    inner = pval_clusters[j]
    count_inner = 0
    for k in inner:
        current_name = cell_type[count_inner]
        count_inner +=1
        count = 0
        #16 arrays, 4 for each
        for l in k:
            dict_results["names"].append(names[count])
            dict_results["cell_type"].append(j)
            dict_results["compared_cell_type"].append(current_name)
            dict_results["pvals"].append(l)
            count += 1

dict_results["log2FC"] = lfc_list

new_df = pd.DataFrame.from_dict(dict_results)

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

and so on for all others

Plasma cell array

Plasma cell (long-lived; higher BLIMP1) array

Plasmablasts array

'''

#create output file
with open("ingest/output_a_vs_b_de", "w") as external_file:
    print(new_df.to_string(), file = external_file)
=======
# Accessing results

#outputs shape of p values
np.set_printoptions(precision=3)
print("shape of p-values: %s" % str(test.pval.shape))


#brings up pvalue
print(test.pval[:,:,0])

#brings up plot 
test.plot_volcano()

