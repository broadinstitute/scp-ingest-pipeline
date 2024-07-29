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
cas_response_output_path = f"{input_path}__cas_response_ontology_aware.json"

def cas_annotate(input_path, output_path):
    """Call CAS for an H5AD, write results; return input AnnData, CAS response
    """
    # TODO (SCP-5715):
    # - Add CAS API token to Google Secrets Manager
    # - Update block below to use GSM
    # - Ensure team removes any local .cas-api-token files and .gitignore entry
    with open(".cas-api-token") as f:
        api_token = f.read().strip()

    cas = CASClient(api_token=api_token)
    print('cas')
    print(cas)

    adata = sc.read_h5ad(input_path)

    # Returns summary results that have substantial post-processing support
    cas_response = cas.annotate_matrix_cell_type_ontology_aware_strategy(
        matrix=adata
    )

    # Returns dataset-specific results
    # cas_response = cas.annotate_matrix_cell_type_summary_statistics_strategy(
    #     matrix=adata
    # )

    with open(output_path, "w") as f:
        f.write(json.dumps(cas_response))

    print(f"Wrote CAS response to: {output_path}")

    return [adata, cas_response]

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

def run(input_path, cas_response_output_path):
    print("Running CAS ingest")
    # adata, cas_response = cas_annotate(input_path, cas_response_output_path)

    # TODO: Remove this block after debugging
    adata = sc.read_h5ad(input_path)
    with open(cas_response_output_path) as f:
        cas_ontology_aware_response = json.loads(f.read())

    insert_cas_ontology_aware_response_into_adata(cas_ontology_aware_response, adata, cl)

    pp.compute_most_granular_top_k_calls_cluster(
        adata=adata,
        cl=cl,
        min_acceptable_score=0.1,  # minimum acceptable score for a call
        cluster_label_obs_column='cluster_label',  # .obs column containing cluster labels
        top_k=3,  # how many top calls to make?
        obs_prefix='cas_cell_type_cluster'  # .obs column to write the top-k calls to
    )

    cas_anndata_output_path = f"cas_annotated__{input_path}"
    adata.write(cas_anndata_output_path)
    print(f"Wrote CAS-annotated AnnData: {cas_anndata_output_path}")

    final_output_path = f"scp_compliant_{cas_anndata_output_path}"
    adata = make_compliant(cas_anndata_output_path)
    adata.write(final_output_path)
    print(f"Wrote SCP-compliant, CAS-annotated AnnData: {final_output_path}")

if __name__ == "__main__":
    run(input_path, cas_response_output_path)
