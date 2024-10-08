"""Cell Annotation Service (CAS) ETL for Single Cell Portal

Context: https://github.com/broadinstitute/scp-ingest-pipeline/pull/353
"""

import json
import re

import scanpy as sc

from cellarium.cas import CASClient
from cellarium.cas._io import suppress_stderr
import cellarium.cas.postprocessing.ontology_aware as pp
from cellarium.cas.postprocessing.cell_ontology import CellOntologyCache
from cellarium.cas.postprocessing import insert_cas_ontology_aware_response_into_adata

# Comment out block below below unless debugging
# from cellarium.cas.visualization import CASCircularTreePlotUMAPDashApp

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
    print('CASClient')
    print(cas)

    adata = sc.read_h5ad(input_path)

    # Returns summary results that have substantial post-processing support
    cas_ontology_aware_response = cas.annotate_matrix_cell_type_ontology_aware_strategy(
        matrix=adata,
        chunk_size=500,
        feature_ids_column_name='gene_ids',
        feature_names_column_name='index'
    )

    # Returns dataset-specific results.
    # Ontology-aware strategy (above, not this) is preferred.
    # cas_response = cas.annotate_matrix_cell_type_summary_statistics_strategy(
    #     matrix=adata
    # )

    with open(output_path, "w") as f:
        f.write(json.dumps(cas_ontology_aware_response))

    adata.write(f"cas_annotated_before_postprocessing__{input_path}")

    print(f"Wrote CAS response to: {output_path}")

    return [adata, cas_ontology_aware_response]

def make_compliant_for_scp(adata):
    """Shim to make AnnData file compliant with SCP metadata schema

    Only use for demo / development purposes.  This is generally commented out
    upstream.
    """
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
        adata.obs['library_preparation_protocol__ontology_label'] = "10x 3' v3"
        adata.obs['sex'] = "female"
    return adata

# Stub for potential later expansion
# def format_as_scp_metadatum(cas_response_path):
#     with open(cas_response_filepath) as f:
#         cas_response = json.loads(f.read())

#     tsv_rows = []
#     for cas_item in cas_response:
#         cell = cas_item["query_cell_id"]
#         first_match = cas_item["matches"][0]
#         annotation_label = cas_item["matches"][0]

def merge_cas_results(adata, cas_ontology_aware_response):
    """Update AnnData with ontology-aware CAS results JSON
    """
    insert_cas_ontology_aware_response_into_adata(cas_ontology_aware_response, adata, cl)

    # Comment out block below unless debugging
    # DASH_SERVER_PORT = 8050
    # with suppress_stderr():
    #     CASCircularTreePlotUMAPDashApp(
    #         adata,  # the AnnData file
    #         cas_ontology_aware_response,  # CAS response
    #         cluster_label_obs_column="cluster_label",  # (optional) The .obs column name containing cluster labels
    #     ).run(port=DASH_SERVER_PORT, debug=False, jupyter_width="100%")

    # TODO: Consider conditionally using this `single` method if (and only if)
    # AnnData lacks "raw annotations" as returned from e.g. Scanpy or Seurat
    # ("0", "1", "2", etc.)
    # pp.compute_most_granular_top_k_calls_single(
    #     adata=adata,
    #     cl=cl,
    #     min_acceptable_score=0.1,  # minimum acceptable score for a call
    #     top_k=3,  # how many top calls to make?
    #     obs_prefix="cell_type_cas"  # .obs column to write the top-k calls to
    # )

    # Use this only if source AnnData has clusterings
    pp.compute_most_granular_top_k_calls_cluster(
        adata=adata,
        cl=cl,
        min_acceptable_score=0.1,  # minimum acceptable score for a call
        cluster_label_obs_column='cluster_label',  # .obs column containing cluster labels
        top_k=3,  # how many top calls to make?
        obs_prefix='cell_type_cas'  # .obs column to write the top-k calls to
    )

    return adata

def trim_cas_adata(adata):
    """Trim CAS AnnData to only annotation labels; omit IDs, scores

    This ensures the AnnData can trivially initialize a polished SCP study.
    The IDs, scores etc. can be easily repopulated via `merge_cas_results` for
    debugging, convenient fuller AnnData output if desired, etc.
    """
    # "Name" is ontology ID, e.g. CL_0000897.  Score is CAS confidence.
    annots_to_omit = ["name", "score"]

    columns = list(adata.obs)
    for to_omit in annots_to_omit:
        for col in columns:
            # E.g. re.match('cas.*_name_\d+', 'cas_cell_type_name_1')
            match = re.match(f".*_cas_.*{to_omit}_\d+", col)
            if match:
                del adata.obs[col]

    return adata

# TODO: Consider porting to SCP Core JS
# def friendlify_cas_adata(adata):
#     """Make CAS annotation names human-friendlier, e.g. "CAS
#     """

#     columns = list(adata.obs)
#     for old_col in columns:
#         if not '_cas_' in old_col:
#             continue
#         # E.g. cell_type_cas_label_42 -> cell_type_(CAS 42)
#         col = re.sub(r'_cas_label_(\d+)', r'_(CAS \1)_', old_col)

#         # E.g. cell_type_(CAS 42) -> cell type (CAS 42)
#         col = col.replace('_', ' ')

#         # E.g. cell type (CAS 42) -> Cell type (CAS 42)
#         col = col[0].upper() + col[1:]

#         # E.g. Cell type (CAS 1) -> Cell type (CAS)
#         col = col.replace('CAS 1)', 'CAS)')

#         adata.obs[col] = adata.obs[old_col]
#         del adata.obs[old_col]

#     return adata

def save_anndata(adata, stem, input_path):
    cas_anndata_output_path = f"{stem}__{input_path}"
    adata.write(cas_anndata_output_path)
    print(f"Wrote AnnData: {cas_anndata_output_path}")

def run(input_path, cas_response_output_path):
    print("Running CAS ingest")
    adata, cas_ontology_aware_response = cas_annotate(input_path, cas_response_output_path)

    # Comment out block below below unless debugging
    # adata = sc.read_h5ad(f"cas_annotated_before_postprocessing__{input_path}")
    # with open(cas_response_output_path) as f:
    #     cas_ontology_aware_response = json.loads(f.read())

    adata = merge_cas_results(adata, cas_ontology_aware_response)
    save_anndata(adata, "cas_annotated_", input_path)

    adata = trim_cas_adata(adata)
    save_anndata(adata, "cas_annotated_trimmed", input_path)

    # Comment out block below below unless debugging or generating demo data
    # TODO: Enable SCP metadata validation exemption for AnnData
    # adata = make_compliant_for_scp(adata)
    # save_anndata(adata, "cas_annotated_trimmed_compliant", input_path)


if __name__ == "__main__":
    run(input_path, cas_response_output_path)
