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

# TODO (SCP-5041): Extract these to CLI arguments

#Pairwise tests between groups
#answers whether a given pair of groups shows differential expression for each gene

#arguments:
#data: Anndata object, data matrix with cells x genes

#grouping: str, column in data, the column that contains cell type labels. alternatively vector containing group labels

#test: type of test, default is z-test. other options: ’wald’, ’lrt’, ’t-test’, ’rank’

#lazy: bool, only possible if test is ztest, if true only evaluated once the user requests pval/coefficients for a specific pair of models

#noise_model: default is nb, specify NONE for wald and t test

h5ad_file = "tests/data/anndata/trimmed_compliant_pbmc3K.h5ad"

adata_file = sc.read_h5ad(h5ad_file)

adata_file.var_names_make_unique()

print(adata_file.obs["louvain"])


test = de.test.pairwise(
    data=adata_file,
    grouping="louvain",
    test="rank",
    lazy=False,
    noise_model=None
 )

# Accessing results

#outputs shape of p values
np.set_printoptions(precision=3)
print("shape of p-values: %s" % str(test.pval.shape))


#brings up pvalue
print(test.pval[:,:,0])

#brings up plot 
test.plot_volcano()

