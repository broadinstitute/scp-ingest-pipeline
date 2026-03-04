"""Tests for `cellxgene_to_scp.py` conversion and validation."""

import sys
from pathlib import Path

import anndata as ad
import pytest

# Ensure project root is on sys.path so we can import the package
sys.path.append(str(Path(__file__).resolve().parent.parent))
from ingest.cellxgene_to_scp import convert_cellxgene_to_scp, validate_scp_fields


TEST_H5AD = Path(__file__).parent / "data" / "anndata" / "cellxgene.human_liver_b_cells.h5ad"


def test_convert_creates_scp_h5ad(tmp_path):
    """Run conversion on a small cellxgene fixture and check output file and columns."""
    out_path = tmp_path / "cellxgene.scp.h5ad"

    # Run conversion
    result = convert_cellxgene_to_scp(str(TEST_H5AD), str(out_path), verbose=False)

    # Path should be returned and file should exist
    assert result == str(out_path)
    assert out_path.exists()

    # Open the result and check required columns exist
    adata = ad.read_h5ad(result)
    obs_cols = set(adata.obs.columns)

    required = {
        'biosample_id', 'donor_id', 'species', 'species__ontology_label',
        'disease', 'disease__ontology_label', 'organ', 'organ__ontology_label',
        'library_preparation_protocol', 'library_preparation_protocol__ontology_label',
        'sex',
    }

    missing = required - obs_cols
    assert not missing, f"Missing required SCP fields: {missing}"

    # sex must be one of allowed values
    assert set(adata.obs['sex'].unique()).issubset({'male', 'female', 'mixed', 'unknown'})


def test_validate_scp_fields_reports_valid(tmp_path):
    """Validate the converted file using the library function."""
    out_path = tmp_path / "cellxgene.scp.h5ad"
    convert_cellxgene_to_scp(str(TEST_H5AD), str(out_path), verbose=False)

    results = validate_scp_fields(str(out_path), verbose=False)

    assert results['valid'] is True
    assert results['n_cells'] > 0
    assert results['n_genes'] > 0
