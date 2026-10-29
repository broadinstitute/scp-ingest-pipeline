# nemo_scp_client

Scripts for turning NeMO (BICCN 10x 3′ v2) H5AD files into
[Single Cell Portal](https://singlecell.broadinstitute.org) (SCP)-ready AnnData,
plus an experiment applying the [Cellarium Cell Annotation Service](https://cellarium.ai)
(CAS) to cell type annotation.

## What's here

```
src/
  run_next_file.py              # Full pipeline for one file (or the next unprocessed one)
  nemo_h5ad_downloader.py       # Step 1: download a file's H5AD
  file_sample_crosstab.py       # Step 2: join file manifest -> sample manifest
  ontology_map_nemo_to_scp.py   # Step 3: map sample metadata to ontology IDs
  project_ontology_to_h5ad.py   # Step 4: write ontology columns into the H5AD's obs
  cluster_adata.py              # Optional: HVG -> PCA -> neighbors -> UMAP -> Leiden
cas.py                          # CAS quickstart tutorial + pseudobulk experiment
requirements.txt                # Python dependencies (see Installation)
data/                           # Manifests (inputs) and all pipeline outputs
docs/                           # Species breakdown and CAS compatibility notes
```

### The NeMO → SCP pipeline

Each NeMO file has an H5AD of raw counts but no SCP metadata. SCP requires
specific obs columns that hold ontology IDs. The pipeline adds them in four steps,
keyed by a NeMO `file_id` such as `nemo:fil-0njjzpd`:

1. **Download** (`nemo_h5ad_downloader.py`). Looks up the file in
   `data/file_manifest_*.tsv`, picks a `gs://` URL when one exists (otherwise
   HTTPS), and saves the file to `data/fil-<id>/<file_name>.h5ad`.
2. **Crosstab** (`file_sample_crosstab.py`). Expands the file's comma-separated
   `sample_id`s against `data/sample_manifest_*.tsv` and writes
   `data/fil-<id>.tsv`.
3. **Ontology mapping** (`ontology_map_nemo_to_scp.py`). Downloads SCP's
   minified ontologies from
   [scp-ingest-pipeline](https://github.com/broadinstitute/scp-ingest-pipeline)
   into `data/.ontology_cache/`. Then it matches each row's free-text values to
   ontology terms (exact label, then synonym, then fuzzy match at ≥ 0.75) and
   writes `data/fil-<id>.scp.tsv` with these added columns:

   | Column | Source | Ontology |
   |---|---|---|
   | `organ` | `sample_anatomies` (`monocyte` falls back to `blood`) | UBERON |
   | `library_preparation_protocol` | `sample_technique` | EFO |
   | `species` | `taxon` | NCBITaxon |
   | `disease` | always `normal` (NeMO has no disease annotation) | PATO |
   | `sex` | `subject_sex` | — |
   | `biosample_id`, `donor_id` | `sample_name` of the `nemo:smp-` row | — |

   Each ontology column also gets a matching `<column>__ontology_label` column.
4. **Projection** (`project_ontology_to_h5ad.py`). Copies the H5AD to
   `data/fil-<id>.scp.h5ad` and adds each column above as a single-category obs
   column. The obs group is edited in place with `h5py`, so `X` is never
   loaded into memory. The step fails if a column has more than one distinct
   value across the file's samples.

Every step skips work whose output already exists, unless you pass `--overwrite`.

### Clustering (optional)

`cluster_adata.py` runs a standard Scanpy workflow (seurat_v3 HVGs, normalize,
log1p, scale, PCA, cosine neighbors, UMAP, Leiden) on raw counts. The PCA
embedding is removed from `obsm` afterwards because SCP's validator treats every
`obsm` entry as 2D/3D visualization coordinates. Output: `data/fil-<id>.clustered.h5ad`.

### CAS experiment (`cas.py`)

This is a script export of the Cellarium CAS quickstart Colab notebook, run
against the 4k-cell PBMC tutorial dataset. It:

- annotates each cell with CAS's ontology-aware strategy and stores the top-k
  per-cell and per-cluster calls in `obs`;
- computes a **pseudobulk** profile per cluster (mean expression), annotates
  those profiles with CAS, and maps each label back to its cells as
  `pseudobulk_cas_label`.

Plots use the non-interactive `Agg` backend, so the script never opens windows.
See [docs/species_breakdown_and_annotation_notes.md](docs/species_breakdown_and_annotation_notes.md)
for why CAS currently only fits the human files in the NeMO manifest.

## Installation

You need Python 3.10+.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

- The steps that only process TSVs (download, crosstab, ontology mapping) use
  only the standard library; `requirements.txt` covers the rest.
  `project_ontology_to_h5ad.py` needs `h5py` and `numpy`. `cluster_adata.py`
  needs `scanpy`, `leidenalg`, and `scikit-misc` (for seurat_v3 HVG selection).
- **`gs://` downloads** use the [`gcloud` CLI](https://cloud.google.com/sdk/docs/install).
  The NeMO buckets are requester-pays, so you need a billing project. The
  downloader uses the first one it finds among: `--billing-project`, the
  `GCLOUD_BILLING_PROJECT` environment variable, and `gcloud config get-value project`.

  ```bash
  gcloud auth login
  export GCLOUD_BILLING_PROJECT=my-gcp-project
  ```
- **CAS:** `cas.py` installs `cellarium-cas[vis]` itself with `pip` every time it
  runs. It also needs a CAS API token in a `.cas-api-token` file in the project
  root. This file is a secret, so don't share or commit it.

## Usage

Run all scripts from the project root. Each one prints the path of its output
on stdout. By default, inputs and outputs live in `data/`.

### Process files end to end

```bash
# The next file in the manifest that has no .scp.h5ad yet
python src/run_next_file.py

# A specific file
python src/run_next_file.py nemo:fil-0njjzpd

# Re-run every step, replacing existing outputs
python src/run_next_file.py nemo:fil-0njjzpd --overwrite
```

To process the whole manifest, run `run_next_file.py` in a loop until it exits
non-zero with "every file_id in the manifest already has a .scp.h5ad".

### Run individual steps

```bash
python src/nemo_h5ad_downloader.py nemo:fil-0njjzpd [--billing-project P] [--timeout-seconds 60]
python src/file_sample_crosstab.py nemo:fil-0njjzpd      # omit the id to write all 217
python src/ontology_map_nemo_to_scp.py fil-0njjzpd.tsv   # accepts multiple TSVs
python src/project_ontology_to_h5ad.py nemo:fil-0njjzpd
```

### Cluster an SCP H5AD

```bash
python src/cluster_adata.py data/fil-0njjzpd.scp.h5ad
# For large files: keep the matrix sparse and skip UMAP if you only need labels
python src/cluster_adata.py data/fil-0njjzpd.scp.h5ad --no-zero-center --no-umap --resolution 2
```

Run `--help` on any script to see all of its options.

### Run the CAS experiment

```bash
python cas.py
```

The script downloads `pbmc_10x_v3_4k.h5ad` (~74 MB) to the current directory and
sends it to CAS. It saves no files: results stay on the in-memory `adata` /
`adata_bulk` objects. To keep them, run it interactively (for example
`python -i cas.py`) or add an `adata.write(...)` call.

## Output files

For each file `fil-<id>`, `data/` contains:

| File | Produced by |
|---|---|
| `fil-<id>/<file_name>.h5ad` | download (the original NeMO H5AD) |
| `fil-<id>.tsv` | crosstab |
| `fil-<id>.scp.tsv` | ontology mapping |
| `fil-<id>.scp.h5ad` | projection — **this is the file to upload to SCP** |
| `fil-<id>.clustered.h5ad` | `cluster_adata.py` (optional) |
