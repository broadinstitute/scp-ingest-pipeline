"""Cell Annotation Service (CAS) ETL for Single Cell Portal

Context: https://github.com/broadinstitute/scp-ingest-pipeline/pull/353
"""

from cellarium.cas import CASClient
import json

from cellarium.cas._io import suppress_stderr
import cellarium.cas.postprocessing.ontology_aware as pp
from cellarium.cas.postprocessing.cell_ontology import CellOntologyCache
from cellarium.cas.postprocessing import insert_cas_ontology_aware_response_into_adata


import numpy as np
import pandas as pd
import scanpy as sc

with suppress_stderr():
    cl = CellOntologyCache()

input_path = "pbmc_10x_v3_4k.h5ad"
output_path = f"{input_path}__cas_response_2.json"

def run(input_path, output_path):
    # TODO (SCP-5715):
    # - Add CAS API token to Google Secrets Manager
    # - Update block below to use GSM
    # - Ensure team removes any local .cas-api-token files and .gitignore entry
    with open(".cas-api-token") as f:
        api_token = f.read().strip()

    print('api_token')
    print(api_token)
    cas = CASClient(api_token=api_token)
    print('cas')
    print(cas)

    adata = sc.read_h5ad(input_path)

    response = cas.annotate_matrix_cell_type_summary_statistics_strategy(
        matrix=adata
    )
    print("CAS response")
    print(response)
    print("")
    print("")
    print("")
    print("adata")
    print(adata)

    with open(output_path, "w") as f:
        f.write(json.dumps(response))

    print(f"Wrote CAS response to: {output_path}")
    return [response, adata]

def make_compliant(input_path):
    adata = sc.read_h5ad(input_path)
    # Uncomment only if running in debug with non-compliant file
    if 'biosample_id' not in adata.obs:
        adata.obs['biosample_id'] = "sample-1"
        adata.obs['donor_id'] = "donor-1"
        adata.obs['species'] = "NCBITaxon_9606"
        adata.obs['species__ontology_label'] = "Homo sapiens"
        adata.obs['disease'] = "PATO_0000461"
        adata.obs['disease__ontology_label'] = "normal"
        adata.obs['organ'] = "UBERON_0000178"
        adata.obs['organ__ontology_label'] = "blood"
        adata.obs['library_preparation_protocol'] = "EFO_0030059"
        adata.obs['library_preparation_protocol__ontology_label'] = "10x multiome"
        adata.obs['sex'] = "female"
        # adata.write(input_path)
        # print(f"Updated AnnData to be SCP-compliant: {cas_response_filepath}")
    return adata

def format_as_scp_metadatum(cas_response_path):
    with open(cas_response_filepath) as f:
        cas_response = json.loads(f.read())

    tsv_rows = []
    for cas_item in cas_response:
        cell = cas_item["query_cell_id"]
        first_match = cas_item["matches"][0]
        annotation_label = cas_item["matches"][0]

if __name__ == "__main__":
    run(input_path)
