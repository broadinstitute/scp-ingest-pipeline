# Species/anatomy breakdown, CAS compatibility notes, and metadata projection

Covers the NeMO/SCORCH BICCN 10x 3' file manifest (`data/file_manifest_nemo_archive__biccn__10x_3prime_v2.tsv`,
217 files) and the sample manifest it joins against
(`data/sample_manifest_nemo_archive__biccn__10x_3prime_v2.tsv`).

## Species breakdown

Every one of the 217 files resolves to exactly one species (no missing or mixed-taxon files):

| Taxon | Files | % of manifest |
|---|---|---|
| rhesus macaque | 114 | 52.5% |
| Norway rat | 39 | 18.0% |
| human | 37 | 17.1% |
| house mouse | 27 | 12.4% |

## Macaque anatomy breakdown

The 114 macaque files split across these `sample_anatomies` values:

| Anatomy | Files | Basal ganglia? |
|---|---|---|
| nucleus accumbens | 22 | yes |
| caudate nucleus | 19 | yes |
| white matter of frontal lobe | 20 | no |
| spleen | 21 | no (non-brain) |
| caudate-putamen | 8 | yes |
| cerebrospinal fluid | 8 | no |
| white matter (unspecified) | 8 | no |
| monocyte | 5 | no (non-brain) |
| middle frontal gyrus | 3 | no |

**49 of 217 files (22.6%) are macaque basal ganglia** (nucleus accumbens + caudate nucleus +
caudate-putamen) — 43.0% of the 114 macaque files. The remaining 65 macaque files (57%) are
outside basal ganglia (frontal cortex/white matter, CSF, spleen, monocyte).

## Cellarium CAS (Cell Annotation Service) compatibility notes

CAS requires: (1) raw integer counts in `X`, and (2) both gene symbols and Ensembl gene IDs in `.var`.
Our SCP H5ADs satisfy (1) — Optimus output is raw, unnormalized counts (see `X` integer-valued check,
no `layers`, `uns['pipeline_version'] = 'Optimus_v8.0.4'`). (2) is satisfied structurally for the mouse
file but not the macaque file:

- `fil-0njjzpd` (mouse): has real Ensembl IDs (`ENSMUSG...`) alongside gene symbols.
- `fil-03vd8jx` (macaque): its `ensembl_ids` column is **not** real Ensembl IDs — it's the gene symbol
  duplicated (`PGBD2`, `ZNF692`, ...), because the macaque reference annotation doesn't provide Ensembl
  gene IDs. This file would fail CAS's ID requirement as submitted.

**Direct confirmation from the CAS team (personal communication):** Fedor Grab (CAS engineer) says
*"I think mouse is very likely, but we haven't included any other organisms in our immediate
roadmap."* Stephen Fleming (CAS scientist) agrees. They put a **3-6 month timeline** on mouse support
landing. So:

- Mouse support is a likely addition on a roughly 3-6 month timeline (i.e. sometime around late 2026
  to early 2027 from when this was communicated) — not shipped yet, and still an estimate rather than a
  committed date.
- Macaque and rat are explicitly **not** on the immediate roadmap — this isn't a gap that's about to
  close any time soon.

**Bottom line:** with CAS as it exists today, only the 37 human files (17.1% of the manifest) are
plausible candidates, and even those should be checked against `list_models()` once API access is set
up, since availability can change. Mouse (27 files, 12.4%) has a real, team-estimated path to becoming
viable within the next ~3-6 months, bringing potential coverage to 64 files (29.5%) once it ships. But
macaque and rat (153 files, 70.5% combined) have no shipped model today and, per the team directly,
aren't on the near-term roadmap at all — this majority of the manifest should be treated as out of
scope for CAS for the foreseeable future.

## Projecting sample metadata onto the SCORCH H5AD files

The file manifest and sample manifest are joined on `sample_id`, and that join is what lets sample-level
metadata (anatomy, technique, taxon, subject info, study info) be attached to a file's H5AD as new `obs`
columns. The linkage:

1. **File manifest → sample IDs.** Each file manifest row (`file_id`, `file_name`, `sample_id`, ...) has a
   `sample_id` column containing a **comma-separated list** of sample IDs tied to that file — typically a
   mix of `nemo:lib-...` (library/lane-level) and `nemo:smp-...` (biosample-level) records, since one
   H5AD can be built from several sequencing lanes of the same biological sample.
2. **Sample manifest → metadata.** The sample manifest has one row per `sample_id`, keyed on that same
   ID, with the actual metadata columns: `sample_anatomies`, `sample_modality`, `sample_technique`,
   `sample_specimen_type`, `subject_id`, `taxon`, `subject_age`, `subject_sex`, `subject_race`,
   `subject_ethnicity`, `study_id`, `study_name`, `study_program_name`, `project_grant`, `study_labs`.
3. **Crosstab join** (`src/file_sample_crosstab.py`): for a given `file_id`, split its `sample_id` list and
   look up each ID in the sample manifest, emitting one row per linked sample (`fil-<id>.tsv`). This
   preserves every lib- and smp-level record for that file, side by side.
4. **Ontology resolution** (`src/ontology_map_nemo_to_scp.py`): maps free-text values in that crosstab to
   controlled ontology terms — `sample_anatomies` → UBERON ID/label, `sample_technique` → EFO
   (`library_preparation_protocol`) ID/label, `taxon` → NCBITaxon ID/label — and derives `Biosample_ID`/
   `Donor_ID` from the biosample-level row's `sample_name`. Output is `fil-<id>.scp.tsv`.
5. **Projection onto the H5AD** (`src/project_ontology_to_h5ad.py`): because every cell in one H5AD comes
   from the same biosample (even when split across multiple lib- rows for different lanes), each of those
   8 ontology columns resolves to a **single constant value per file** — the resolver collects the
   distinct non-empty values for each column across all of a file's crosstab rows and raises if they
   disagree (that's the correctness check: it guarantees we aren't silently projecting a mismatched
   value). The constant is then written as a new categorical `obs` column directly via `h5py` on a copy
   of the H5AD (one category, all-zero codes for every cell) rather than loading and rewriting the full
   expression matrix — which matters at this data's scale (up to ~1.6M cells / tens of GB per file).
   Output is `fil-<id>.scp.h5ad`, an SCP-ready AnnData with ontology-mapped `obs` columns on every cell.
6. **Orchestration** (`src/run_next_file.py`): chains all of the above (download → crosstab → ontology
   mapping → projection) for the next file_id in manifest order that doesn't yet have a `.scp.h5ad`.

This join mechanism generalizes beyond the 8 ontology columns currently projected — any sample manifest
column (e.g. `subject_age`, `subject_sex`, `study_name`) could be projected the same way, since the same
`file_id → sample_id → sample metadata` linkage and same-value-across-a-file assumption holds for all of
them, not just the ontology-mapped ones.
