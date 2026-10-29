"""Run the full NEMO -> SCP pipeline for the next unprocessed file in the manifest.

For a file_id, this runs, in order:
1. nemo_h5ad_downloader.download_h5ad       -- fetch the H5AD
2. file_sample_crosstab.write_all_file_sample_tables -- build the sample crosstab TSV
3. ontology_map_nemo_to_scp.annotate_file_tsv         -- add ontology columns to the TSV
4. project_ontology_to_h5ad.project_ontology_to_h5ad  -- project those columns onto the H5AD

The result is `<file_id>.scp.h5ad` in the data directory: an SCP-ready AnnData
file with ontology-mapped obs columns.

If no file_id is given, the next one is picked by walking the file manifest in
order and choosing the first row whose `<file_id>.scp.h5ad` does not yet exist.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from file_sample_crosstab import (
    DEFAULT_FILE_MANIFEST as CROSSTAB_DEFAULT_FILE_MANIFEST,
    DEFAULT_SAMPLE_MANIFEST,
    write_all_file_sample_tables,
)
from nemo_h5ad_downloader import DEFAULT_MANIFEST, download_h5ad, normalize_file_id
from ontology_map_nemo_to_scp import annotate_file_tsv
from project_ontology_to_h5ad import project_ontology_to_h5ad


DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def iter_manifest_file_ids(manifest_path: str | Path = DEFAULT_MANIFEST) -> list[str]:
    """Return every file_id in the file manifest, in manifest order."""

    with Path(manifest_path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [(row.get("file_id") or "").strip() for row in reader if row.get("file_id")]


def find_next_file_id(
    manifest_path: str | Path = DEFAULT_MANIFEST,
    data_dir: str | Path = DATA_DIR,
) -> str | None:
    """Return the first file_id in manifest order without a `.scp.h5ad` yet."""

    data_dir = Path(data_dir)
    for file_id in iter_manifest_file_ids(manifest_path):
        output_path = data_dir / f"{normalize_file_id(file_id)}.scp.h5ad"
        if not output_path.exists():
            return file_id

    return None


def run_pipeline(
    file_id: str,
    file_manifest_path: str | Path = DEFAULT_MANIFEST,
    sample_manifest_path: str | Path = DEFAULT_SAMPLE_MANIFEST,
    data_dir: str | Path = DATA_DIR,
    overwrite: bool = False,
    timeout_seconds: int = 60,
    billing_project: str | None = None,
) -> Path:
    """Download, crosstab, annotate, and project a file_id into an SCP-ready H5AD."""

    data_dir = Path(data_dir)

    download_h5ad(
        file_id,
        manifest_path=file_manifest_path,
        output_root=data_dir,
        overwrite=overwrite,
        timeout_seconds=timeout_seconds,
        billing_project=billing_project,
    )

    crosstab_paths = write_all_file_sample_tables(
        file_manifest_path=file_manifest_path,
        sample_manifest_path=sample_manifest_path,
        output_root=data_dir,
        overwrite=overwrite,
        file_id=file_id,
    )
    crosstab_path = crosstab_paths[0]

    annotate_file_tsv(crosstab_path, overwrite=overwrite)

    return project_ontology_to_h5ad(
        file_id,
        manifest_path=file_manifest_path,
        root=data_dir,
        overwrite=overwrite,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the full NEMO -> SCP pipeline for a file_id, or the next "
        "unprocessed one if none is given."
    )
    parser.add_argument(
        "file_id",
        nargs="?",
        default=None,
        help="File id such as nemo:fil-03vd8jx. If omitted, the next file_id in "
        "the manifest without an existing <file_id>.scp.h5ad is used.",
    )
    parser.add_argument(
        "--file-manifest",
        default=str(CROSSTAB_DEFAULT_FILE_MANIFEST),
        help="Path to the file manifest TSV",
    )
    parser.add_argument(
        "--sample-manifest",
        default=str(DEFAULT_SAMPLE_MANIFEST),
        help="Path to the sample manifest TSV",
    )
    parser.add_argument(
        "--data-dir",
        default=str(DATA_DIR),
        help="Directory to read inputs from and write outputs to",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Rerun every stage even if its output already exists",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=60,
        help="Socket timeout for stalled H5AD downloads",
    )
    parser.add_argument(
        "--billing-project",
        default=None,
        help="Billing project for requester-pays gs:// downloads",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    file_id = args.file_id
    if file_id is None:
        file_id = find_next_file_id(manifest_path=args.file_manifest, data_dir=args.data_dir)
        if file_id is None:
            print("error: every file_id in the manifest already has a .scp.h5ad", file=sys.stderr)
            return 1
        print(f"next file_id: {file_id}", file=sys.stderr)

    try:
        output_path = run_pipeline(
            file_id,
            file_manifest_path=args.file_manifest,
            sample_manifest_path=args.sample_manifest,
            data_dir=args.data_dir,
            overwrite=args.overwrite,
            timeout_seconds=args.timeout_seconds,
            billing_project=args.billing_project,
        )
    except (KeyError, ValueError, FileNotFoundError, OSError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
