"""Project ontology columns from a `.scp.tsv` crosstab into its H5AD file.

For a given file_id, this reads `<file_id>.scp.tsv` (the output of
ontology_map_nemo_to_scp.py), resolves the constant value of each ontology
column added by that script (e.g. `organ`), and writes
those values as new per-cell obs columns onto a copy of the file's H5AD,
named `<file_id>.scp.h5ad`.

Each ontology column is expected to hold a single distinct non-empty value
across all crosstab rows for a file_id, since every cell in one H5AD comes
from the same sample/library. A column with conflicting non-empty values
raises an error rather than guessing which one is correct, and a column with
no non-empty value is skipped with a warning rather than written as an
empty-string category.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
from pathlib import Path

import h5py
import numpy as np

try:
    from nemo_h5ad_downloader import load_manifest_record, normalize_file_id
    from ontology_map_nemo_to_scp import ONTOLOGY_COLUMNS
# when importing into AnnDataIngestor, we need to use paths from ingest root
except ImportError:
    from nemo_scp_client.src.nemo_h5ad_downloader import load_manifest_record, normalize_file_id
    from nemo_scp_client.src.ontology_map_nemo_to_scp import ONTOLOGY_COLUMNS

DATA_DIR = Path(__file__).resolve().parent.parent / "data"

DEFAULT_MANIFEST = DATA_DIR / "file_manifest_nemo_archive__biccn__10x_3prime_v2.tsv"


def resolve_h5ad_path(
    file_id: str,
    manifest_path: str | Path = DEFAULT_MANIFEST,
    root: str | Path = DATA_DIR,
) -> Path:
    """Return the path to the downloaded H5AD for a file_id."""

    record = load_manifest_record(manifest_path, file_id)
    return Path(root) / normalize_file_id(record.file_id) / record.file_name


def load_ontology_values(scp_tsv_path: str | Path) -> dict[str, str]:
    """Resolve one constant value per ontology column from a `.scp.tsv`."""

    scp_tsv_path = Path(scp_tsv_path)
    values: dict[str, set[str]] = {column: set() for column in ONTOLOGY_COLUMNS}

    with scp_tsv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{scp_tsv_path} has no header row")

        missing = [column for column in ONTOLOGY_COLUMNS if column not in reader.fieldnames]
        if missing:
            raise ValueError(f"{scp_tsv_path} is missing ontology columns: {missing}")

        for row in reader:
            for column in ONTOLOGY_COLUMNS:
                value = (row.get(column) or "").strip()
                if value:
                    values[column].add(value)

    resolved: dict[str, str] = {}
    for column, distinct_values in values.items():
        if len(distinct_values) > 1:
            raise ValueError(
                f"Ambiguous value for {column!r} in {scp_tsv_path}: {sorted(distinct_values)}"
            )
        resolved[column] = next(iter(distinct_values), "")

    return resolved


def write_obs_categorical_column(obs_group: h5py.Group, column: str, value: str, n_obs: int) -> None:
    """Write a single-category obs column, matching AnnData's categorical encoding."""

    if column in obs_group:
        del obs_group[column]

    column_group = obs_group.create_group(column)
    column_group.attrs["encoding-type"] = "categorical"
    column_group.attrs["encoding-version"] = "0.2.0"
    column_group.attrs["ordered"] = False
    column_group.create_dataset("categories", data=np.array([value], dtype=h5py.string_dtype()))
    column_group.create_dataset("codes", data=np.zeros(n_obs, dtype="int8"))


def project_ontology_to_h5ad(
    file_id: str,
    manifest_path: str | Path = DEFAULT_MANIFEST,
    root: str | Path = DATA_DIR,
    scp_tsv_path: str | Path | None = None,
    h5ad_path: str | Path | None = None,
    output_path: str | Path | None = None,
    overwrite: bool = False,
) -> Path:
    """Copy a file's H5AD and add its ontology columns as new obs columns."""

    normalized_id = normalize_file_id(file_id)
    root = Path(root)

    if scp_tsv_path is None:
        scp_tsv_path = root / f"{normalized_id}.scp.tsv"
    if h5ad_path is None:
        h5ad_path = resolve_h5ad_path(file_id, manifest_path=manifest_path, root=root)
    if output_path is None:
        output_path = root / f"{normalized_id}.scp.h5ad"

    scp_tsv_path = Path(scp_tsv_path)
    h5ad_path = Path(h5ad_path)
    output_path = Path(output_path)

    if not scp_tsv_path.exists():
        raise FileNotFoundError(f"Crosstab TSV not found: {scp_tsv_path}")
    if not h5ad_path.exists():
        raise FileNotFoundError(f"H5AD not found: {h5ad_path}")

    if output_path.exists() and not overwrite:
        return output_path

    ontology_values = load_ontology_values(scp_tsv_path)

    empty_columns = sorted(column for column, value in ontology_values.items() if not value)
    if empty_columns:
        print(
            f"warning: no non-empty value for {', '.join(empty_columns)} in "
            f"{scp_tsv_path}; skipping those obs columns",
            file=sys.stderr,
        )
        for column in empty_columns:
            del ontology_values[column]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(h5ad_path, output_path)

    with h5py.File(output_path, "r+") as handle:
        obs_group = handle["obs"]
        n_obs = obs_group[obs_group.attrs["_index"]].shape[0]

        for column, value in ontology_values.items():
            write_obs_categorical_column(obs_group, column, value, n_obs)

        existing_order = list(obs_group.attrs["column-order"])
        new_order = existing_order + [
            c for c in ONTOLOGY_COLUMNS if c not in existing_order and c in ontology_values
        ]
        obs_group.attrs["column-order"] = np.array(new_order, dtype=object)

    return output_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Project ontology columns from a file's crosstab .scp.tsv into its H5AD."
    )
    parser.add_argument("file_id", help="File id such as nemo:fil-0njjzpd")
    parser.add_argument(
        "--manifest",
        default=str(DEFAULT_MANIFEST),
        help="Path to the file manifest TSV, used to locate the downloaded H5AD",
    )
    parser.add_argument(
        "--root",
        default=str(DATA_DIR),
        help="Directory containing the fil-<id>/ H5AD folder and fil-<id>.scp.tsv",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Rewrite the output H5AD if it already exists",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        output_path = project_ontology_to_h5ad(
            args.file_id,
            manifest_path=args.manifest,
            root=args.root,
            overwrite=args.overwrite,
        )
    except (KeyError, ValueError, FileNotFoundError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
