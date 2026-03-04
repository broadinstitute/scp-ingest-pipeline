#!/usr/bin/env python3
"""
CELLxGENE to Single Cell Portal (SCP) H5AD Conversion Module

Converts H5AD files from CELLxGENE schema (v3.0.0) to SCP metadata convention.

CELLxGENE schema: https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md
SCP schema: https://singlecell.zendesk.com/hc/en-us/articles/360060609852-Required-metadata

Key differences and mappings:
- CELLxGENE uses `_ontology_term_id` suffix; SCP uses plain field or `__ontology_label`
- CELLxGENE: organism_ontology_term_id -> SCP: species (with __ontology_label)
- CELLxGENE: tissue_ontology_term_id -> SCP: organ (with __ontology_label)
- CELLxGENE: assay_ontology_term_id -> SCP: library_preparation_protocol (with __ontology_label)
- CELLxGENE: self_reported_ethnicity_ontology_term_id -> SCP: ethnicity (no direct SCP equivalent, but preserved)
- SCP requires: biosample_id, donor_id, species, disease, organ, library_preparation_protocol, sex
- SCP sex field uses controlled vocabulary: ["male", "female", "mixed", "unknown"]
"""

import anndata as ad
import pandas as pd
from pathlib import Path
from typing import Optional, Dict, Any
import sys


# CELLxGENE -> SCP field mappings
# Format: (cxg_ontology_id_field, cxg_label_field, scp_ontology_field, scp_label_field)
FIELD_MAPPINGS = {
    # organism -> species
    'organism': ('organism_ontology_term_id', 'organism', 'species', 'species__ontology_label'),
    # tissue -> organ
    'tissue': ('tissue_ontology_term_id', 'tissue', 'organ', 'organ__ontology_label'),
    # assay -> library_preparation_protocol
    'assay': ('assay_ontology_term_id', 'assay', 'library_preparation_protocol', 'library_preparation_protocol__ontology_label'),
    # disease (same name, different format)
    'disease': ('disease_ontology_term_id', 'disease', 'disease', 'disease__ontology_label'),
    # cell_type (recommended in SCP)
    'cell_type': ('cell_type_ontology_term_id', 'cell_type', 'cell_type', 'cell_type__ontology_label'),
}

# Sex mapping from PATO ontology terms to SCP controlled vocabulary
SEX_ONTOLOGY_TO_SCP = {
    'PATO:0000383': 'female',
    'PATO:0000384': 'male',
    'PATO:0001340': 'mixed',  # hermaphrodite - map to mixed
    'unknown': 'unknown',
}


def normalize_ontology_id(value: str) -> str:
    """
    Normalize ontology ID format.
    CELLxGENE uses 'PREFIX:ID' format (e.g., 'CL:0000000')
    SCP also uses 'PREFIX:ID' but sometimes with underscore in some contexts
    """
    if pd.isna(value) or value in ('', 'nan', 'None'):
        return ''
    value = str(value).strip()
    # Ensure colon format (not underscore)
    if '_' in value and ':' not in value:
        parts = value.split('_', 1)
        if parts[0].isalpha():
            return f"{parts[0]}:{parts[1]}"
    return value


def map_sex_to_scp(sex_ontology_term_id: str, sex_label: str) -> str:
    """
    Map CELLxGENE sex_ontology_term_id to SCP controlled vocabulary.
    SCP requires: ["male", "female", "mixed", "unknown"]
    """
    if pd.isna(sex_ontology_term_id) or sex_ontology_term_id in ('', 'nan', 'None', 'unknown'):
        return 'unknown'
    
    sex_id = str(sex_ontology_term_id).strip()
    
    # Direct mapping from known PATO terms
    if sex_id in SEX_ONTOLOGY_TO_SCP:
        return SEX_ONTOLOGY_TO_SCP[sex_id]
    
    # Try to infer from label
    if not pd.isna(sex_label):
        label_lower = str(sex_label).lower().strip()
        if label_lower in ('male', 'female', 'mixed', 'unknown'):
            return label_lower
        if 'female' in label_lower:
            return 'female'
        if 'male' in label_lower:
            return 'male'
    
    return 'unknown'


def derive_biosample_id(obs: pd.DataFrame) -> pd.Series:
    """
    Derive biosample_id from available columns.
    Priority: suspension_uuid > sample > source_dataset + donor_id
    """
    if 'biosample_id' in obs.columns:
        return obs['biosample_id']
    
    # Try common sample identifier columns
    for col in ['suspension_uuid', 'sample', 'sample_id', 'batch']:
        if col in obs.columns:
            return obs[col].astype(str)
    
    # Create from donor_id + tissue combination if available
    if 'donor_id' in obs.columns and 'tissue' in obs.columns:
        return obs['donor_id'].astype(str) + '_' + obs['tissue'].astype(str)
    
    # Fallback: use donor_id or generate placeholder
    if 'donor_id' in obs.columns:
        return obs['donor_id'].astype(str)
    
    # Last resort: create unique IDs
    return pd.Series(['biosample_' + str(i) for i in range(len(obs))], index=obs.index)


def convert_cellxgene_to_scp(
    input_path: str,
    output_path: Optional[str] = None,
    verbose: bool = True,
) -> str:
    """
    Convert a CELLxGENE H5AD file to SCP format.
    
    Args:
        input_path: Path to input CELLxGENE H5AD file
        output_path: Path to output SCP H5AD file (default: input.scp.h5ad)
        verbose: Print progress messages
    
    Returns:
        Path to output file
    """
    input_path = Path(input_path)
    if output_path is None:
        output_path = input_path.with_suffix('.scp.h5ad')
    else:
        output_path = Path(output_path)
    
    if verbose:
        print(f"Converting: {input_path}")
        print(f"Output: {output_path}")
    
    # Load the file
    if verbose:
        print("Loading H5AD file...")
    adata = ad.read_h5ad(input_path)
    
    if verbose:
        print(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")
        print(f"Original obs columns: {list(adata.obs.columns)}")
    
    # Create a mapping log
    changes = []
    
    # === REQUIRED SCP FIELDS ===
    
    # 1. species (from organism_ontology_term_id)
    # CELLxGENE v3+ may store organism in uns if uniform across all cells
    if 'organism_ontology_term_id' in adata.obs.columns:
        adata.obs['species'] = adata.obs['organism_ontology_term_id'].apply(normalize_ontology_id)
        changes.append("species <- obs.organism_ontology_term_id")
    elif 'organism_ontology_term_id' in adata.uns:
        adata.obs['species'] = normalize_ontology_id(adata.uns['organism_ontology_term_id'])
        changes.append(f"species <- uns.organism_ontology_term_id ({adata.uns['organism_ontology_term_id']})")
    
    if 'organism' in adata.obs.columns:
        adata.obs['species__ontology_label'] = adata.obs['organism'].astype(str)
        changes.append("species__ontology_label <- obs.organism")
    elif 'organism' in adata.uns:
        adata.obs['species__ontology_label'] = str(adata.uns['organism'])
        changes.append(f"species__ontology_label <- uns.organism ({adata.uns['organism']})")
    
    # 2. disease (already named disease_ontology_term_id in CELLxGENE)
    if 'disease_ontology_term_id' in adata.obs.columns:
        adata.obs['disease'] = adata.obs['disease_ontology_term_id'].apply(normalize_ontology_id)
        changes.append("disease <- disease_ontology_term_id (normalized)")
    if 'disease' in adata.obs.columns and 'disease__ontology_label' not in adata.obs.columns:
        # CELLxGENE has 'disease' as label, SCP needs 'disease__ontology_label'
        # But we just overwrote disease with ontology ID, so we need the original label
        pass
    # Get disease label from the original 'disease' column before we overwrote it
    # We need to handle this more carefully
    
    # 3. organ (from tissue_ontology_term_id)
    if 'tissue_ontology_term_id' in adata.obs.columns:
        adata.obs['organ'] = adata.obs['tissue_ontology_term_id'].apply(normalize_ontology_id)
        changes.append("organ <- tissue_ontology_term_id")
    if 'tissue' in adata.obs.columns:
        adata.obs['organ__ontology_label'] = adata.obs['tissue'].astype(str)
        changes.append("organ__ontology_label <- tissue")
    
    # 4. library_preparation_protocol (from assay_ontology_term_id)
    if 'assay_ontology_term_id' in adata.obs.columns:
        adata.obs['library_preparation_protocol'] = adata.obs['assay_ontology_term_id'].apply(normalize_ontology_id)
        changes.append("library_preparation_protocol <- assay_ontology_term_id")
    if 'assay' in adata.obs.columns:
        adata.obs['library_preparation_protocol__ontology_label'] = adata.obs['assay'].astype(str)
        changes.append("library_preparation_protocol__ontology_label <- assay")
    
    # 5. sex (from sex_ontology_term_id, convert to controlled vocab)
    if 'sex_ontology_term_id' in adata.obs.columns:
        sex_label_col = 'sex' if 'sex' in adata.obs.columns else None
        adata.obs['sex'] = adata.obs.apply(
            lambda row: map_sex_to_scp(
                row.get('sex_ontology_term_id', 'unknown'),
                row.get('sex', '') if sex_label_col else ''
            ),
            axis=1
        )
        changes.append("sex <- sex_ontology_term_id (mapped to SCP vocab)")
    elif 'sex' not in adata.obs.columns:
        adata.obs['sex'] = 'unknown'
        changes.append("sex <- 'unknown' (not present)")
    
    # 6. donor_id (should already exist in CELLxGENE 3.0)
    if 'donor_id' not in adata.obs.columns:
        adata.obs['donor_id'] = 'unknown'
        changes.append("donor_id <- 'unknown' (not present)")
    
    # 7. biosample_id (derive from available data)
    adata.obs['biosample_id'] = derive_biosample_id(adata.obs)
    changes.append("biosample_id <- derived from available columns")
    
    # === RECOMMENDED SCP FIELDS ===
    
    # cell_type (recommended)
    if 'cell_type_ontology_term_id' in adata.obs.columns:
        adata.obs['cell_type'] = adata.obs['cell_type_ontology_term_id'].apply(normalize_ontology_id)
        changes.append("cell_type <- cell_type_ontology_term_id")
    if 'cell_type' in adata.obs.columns:
        # CELLxGENE 'cell_type' is the label
        # We need to preserve it as __ontology_label before overwriting
        pass
    
    # === HANDLE DISEASE LABEL ===
    # The CELLxGENE 'disease' column contains the human-readable label
    # We need to save it before the overwrite
    # Re-read to get original disease label
    adata_orig = ad.read_h5ad(input_path)
    if 'disease' in adata_orig.obs.columns:
        adata.obs['disease__ontology_label'] = adata_orig.obs['disease'].astype(str)
        changes.append("disease__ontology_label <- original disease (label)")
    
    # Similarly for cell_type
    if 'cell_type' in adata_orig.obs.columns:
        adata.obs['cell_type__ontology_label'] = adata_orig.obs['cell_type'].astype(str)
        changes.append("cell_type__ontology_label <- original cell_type (label)")
    
    del adata_orig  # Free memory
    
    # === CLEAN UP COLUMN NAMES ===
    # Remove dots from column names (SCP doesn't allow them)
    dotted_cols = [c for c in adata.obs.columns if '.' in c]
    for col in dotted_cols:
        new_col = col.replace('.', '_')
        if new_col not in adata.obs.columns:
            adata.obs[new_col] = adata.obs[col]
            changes.append(f"{new_col} <- {col} (renamed, dot removed)")
        adata.obs = adata.obs.drop(columns=[col])
    
    # === ENSURE STRING TYPES FOR WRITE ===
    str_cols = [
        'species', 'species__ontology_label',
        'disease', 'disease__ontology_label', 
        'organ', 'organ__ontology_label',
        'library_preparation_protocol', 'library_preparation_protocol__ontology_label',
        'sex', 'donor_id', 'biosample_id',
        'cell_type', 'cell_type__ontology_label'
    ]
    for col in str_cols:
        if col in adata.obs.columns:
            adata.obs[col] = adata.obs[col].astype(object).fillna('').astype(str)
    
    # === SANITIZE VAR COLUMNS ===
    var_cols = list(adata.var.columns)
    new_var_cols = {}
    for c in var_cols:
        if c.startswith('_'):
            nc = c.lstrip('_')
            if nc in adata.var.columns:
                nc = 'var_' + nc
            new_var_cols[c] = nc
    if new_var_cols:
        adata.var = adata.var.rename(columns=new_var_cols)
        for old, new in new_var_cols.items():
            changes.append(f"var.{new} <- var.{old} (renamed)")
    
    # === PRESERVE RAW ===
    # SCP expects raw counts to be present if available. Always preserve `adata.raw` when present.
    if hasattr(adata, 'raw') and adata.raw is not None:
        changes.append("Preserved adata.raw")
    
    # === WRITE OUTPUT ===
    if verbose:
        print("\nChanges applied:")
        for c in changes:
            print(f"  - {c}")
        print(f"\nFinal obs columns: {list(adata.obs.columns)}")
        print(f"\nWriting to: {output_path}")
    
    # Write with gzip compression to match CELLxGENE file sizes
    adata.write_h5ad(output_path, compression='gzip')
    
    if verbose:
        print("Done!")
    
    return str(output_path)


def validate_scp_fields(h5ad_path: str, verbose: bool = True) -> Dict[str, Any]:
    """
    Validate that an H5AD file has required SCP metadata fields.
    
    Returns a dict with validation results.
    """
    adata = ad.read_h5ad(h5ad_path)
    obs = adata.obs
    
    required_fields = {
        'biosample_id': 'string - unique identifier for each sample',
        'donor_id': 'string - unique identifier for each donor',
        'species': 'ontology - NCBITaxon identifier',
        'species__ontology_label': 'ontology_label - NCBITaxon label',
        'disease': 'ontology - MONDO or PATO identifier',
        'disease__ontology_label': 'ontology_label - MONDO or PATO label',
        'organ': 'ontology - UBERON identifier',
        'organ__ontology_label': 'ontology_label - UBERON label',
        'library_preparation_protocol': 'ontology - EFO identifier',
        'library_preparation_protocol__ontology_label': 'ontology_label - EFO label',
        'sex': 'controlled list - male/female/mixed/unknown',
    }
    
    recommended_fields = {
        'cell_type': 'ontology - Cell Ontology identifier',
        'cell_type__ontology_label': 'ontology_label - Cell Ontology label',
    }
    
    results = {
        'file': h5ad_path,
        'n_cells': adata.n_obs,
        'n_genes': adata.n_vars,
        'required': {},
        'recommended': {},
        'valid': True,
    }
    
    # Check required fields
    for field, desc in required_fields.items():
        present = field in obs.columns
        if present:
            non_empty = obs[field].dropna().astype(str).str.strip().ne('').sum()
            sample_values = list(obs[field].dropna().unique()[:5])
        else:
            non_empty = 0
            sample_values = []
        
        results['required'][field] = {
            'present': present,
            'non_empty_count': int(non_empty),
            'sample_values': sample_values,
            'description': desc,
        }
        
        if not present or non_empty == 0:
            results['valid'] = False
    
    # Check recommended fields
    for field, desc in recommended_fields.items():
        present = field in obs.columns
        if present:
            non_empty = obs[field].dropna().astype(str).str.strip().ne('').sum()
            sample_values = list(obs[field].dropna().unique()[:5])
        else:
            non_empty = 0
            sample_values = []
        
        results['recommended'][field] = {
            'present': present,
            'non_empty_count': int(non_empty),
            'sample_values': sample_values,
            'description': desc,
        }
    
    if verbose:
        print(f"\n=== SCP Validation Results ===")
        print(f"File: {h5ad_path}")
        print(f"Cells: {results['n_cells']}, Genes: {results['n_genes']}")
        print(f"\nRequired fields:")
        for field, info in results['required'].items():
            status = '✅' if info['present'] and info['non_empty_count'] > 0 else '❌'
            print(f"  {status} {field}: {info['non_empty_count']} non-empty values")
            if info['sample_values']:
                print(f"      Sample: {info['sample_values'][:3]}")
        
        print(f"\nRecommended fields:")
        for field, info in results['recommended'].items():
            status = '✅' if info['present'] and info['non_empty_count'] > 0 else '⚠️'
            print(f"  {status} {field}: {info['non_empty_count']} non-empty values")
            if info['sample_values']:
                print(f"      Sample: {info['sample_values'][:3]}")
        
        print(f"\nOverall: {'✅ VALID' if results['valid'] else '❌ INVALID'}")
    
    return results


def main():
    """CLI entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Convert CELLxGENE H5AD files to SCP format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert a file
  python cellxgene_to_scp.py input.h5ad
  
  # Convert with custom output path
  python cellxgene_to_scp.py input.h5ad -o output.h5ad
  
  # Validate an existing SCP file
  python cellxgene_to_scp.py --validate file.h5ad
"""
    )
    parser.add_argument('input', help='Input H5AD file path')
    parser.add_argument('-o', '--output', help='Output H5AD file path (default: input.scp.h5ad)')
    parser.add_argument('--validate', action='store_true', help='Validate SCP fields only, do not convert')
    parser.add_argument('-q', '--quiet', action='store_true', help='Suppress progress messages')
    
    args = parser.parse_args()
    
    if args.validate:
        validate_scp_fields(args.input, verbose=not args.quiet)
    else:
        convert_cellxgene_to_scp(args.input, args.output, verbose=not args.quiet)


if __name__ == '__main__':
    main()
