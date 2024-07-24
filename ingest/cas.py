"""Cell Annotation Service (CAS) ETL for Single Cell Portal

Context: https://github.com/broadinstitute/scp-ingest-pipeline/pull/353
"""

from cellarium.cas import CASClient

# TODO (SCP-5715):
# - Add CAS API token to Google Secrets Manager
# - Update block below to use GSM
# - Ensure team removes any local .cas-api-token files and .gitignore entry
with open(".cas-api-token") as f:
    api_token = f.read().strip()

cas = CASClient(api_token=api_token)

# Cavaet: processed expression values, as in this normalized / filtered 10x h5 file,
# are not suitable for CAS.  The current code is only intended as a "Hello World".
# The corresponding, commented-out "raw" 10x h5 file below would be best, but only
# if refined to not have ~600K cells analyzed by CAS.
h5_10x_path = "pbmc_unsorted_3k_filtered_feature_bc_matrix.h5"

# h5_10x_path = "pbmc_unsorted_3k_raw_feature_bc_matrix.h5"

response = cas.annotate_matrix_cell_type_summary_statistics_strategy(h5_10x_path)
print(response)
