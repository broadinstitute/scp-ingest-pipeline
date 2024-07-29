import typing as t
from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from enum import Enum
from operator import itemgetter

import numpy as np
import scipy.sparse as sp
from anndata import AnnData

from cellarium.cas.postprocessing.cell_ontology.cell_ontology_cache import CellOntologyCache

# AnnData-related constants
CAS_CL_SCORES_ANNDATA_OBSM_KEY = "cas_cl_scores"
CAS_METADATA_ANNDATA_UNS_KEY = "cas_metadata"

# PhyloXML-related constants
CAS_SCORE_PROPERTY_REF = "CAS:score"
CAS_FRACTION_PROPERTY_REF = "CAS:fraction"
CAS_CL_LABEL_PROPERTY_REF = "CAS:cl_label"

def convert_cas_ontology_aware_response_to_score_matrix(
    adata: AnnData, cas_ontology_aware_response: list, cl: CellOntologyCache
) -> sp.csr_matrix:
    """
    Generate a sparse matrix of CAS ontology-aware scores.

    This function takes an AnnData object, a list of CAS ontology-aware responses, and a CellOntologyCache object
    and generates a sparse matrix of CAS ontology-aware scores. The sparse matrix represents the evidence scores
    of all ontology terms (columns) for each cell (rows).

    :param adata: An AnnData object containing the query cells.
    :type adata: AnnData

    :param cas_ontology_aware_response: A list of CAS ontology-aware responses.
    :type cas_ontology_aware_response: list

    :param cl: A CellOntologyCache object containing the cell ontology information.
    :type cl: CellOntologyCache

    :return: A sparse matrix of CAS ontology-aware scores.
    :rtype: sp.csr_matrix
    """
    row = []
    col = []
    data = []

    print('ok')
    map_length = len(cl.cl_names_to_idx_map)
    cl.cl_names_to_idx_map['unknown'] = map_length
    cl.cl_labels_to_names_map['unknown'] = 'unknown'

    obs_values = adata.obs.index.values
    for obs_idx, cas_cell_response in enumerate(cas_ontology_aware_response):
        # print('cas_cell_response', cas_cell_response)
        # print('cas_cell_response["query_cell_id"]', cas_cell_response["query_cell_id"])
        # print('obs_values[obs_idx]', obs_values[obs_idx])
        assert cas_cell_response["query_cell_id"] == obs_values[obs_idx]
        for match in cas_cell_response["matches"]:
            row.append(obs_idx)
            # eweitz 2024-07-26: Default CAS response lacks `cell_type_ontology_term_id`
            # col.append(cl.cl_names_to_idx_map[match["cell_type_ontology_term_id"]])

            col.append(cl.cl_names_to_idx_map[cl.cl_labels_to_names_map[match["cell_type"]]])
            # eweitz 2024-07-26: Default CAS response lacks `score`
            # data.append(match["score"])
            data.append(match["median_distance"])

    n_obs = len(cas_ontology_aware_response)
    n_cl_names = len(cl.cl_names)
    return sp.coo_matrix((data, (row, col)), shape=(n_obs, n_cl_names)).tocsr()


def insert_cas_ontology_aware_response_into_adata(
    cas_ontology_aware_response: list, adata: AnnData, cl: CellOntologyCache
) -> None:
    """
    Inserts Cellarium CAS ontology aware response into `obsm` property of a provided AnnData file as a
    :class:`scipy.sparse.csr_matrix` named `cas_cl_scores`.

    :param cas_ontology_aware_response: The Cellarium CAS ontology aware response.
    :type cas_ontology_aware_response: list

    :param adata: The AnnData object to insert the response into.
    :type adata: AnnData

    :param cl: The CellOntologyCache object containing cell ontology term names and labels.
    :type cl: CellOntologyCache

    :return: None
    """
    adata.obsm[CAS_CL_SCORES_ANNDATA_OBSM_KEY] = convert_cas_ontology_aware_response_to_score_matrix(
        adata, cas_ontology_aware_response, cl
    )
    adata.uns[CAS_METADATA_ANNDATA_UNS_KEY] = {"cl_names": cl.cl_names, "cl_labels": cl.cl_labels}
