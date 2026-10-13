"""Build per-file sample metadata tables from the NEMO manifests.

For each file manifest row, this script reads the comma-separated sample IDs in
the `sample_id` column, looks up those sample IDs in the sample manifest, and
writes a TSV named `<file_id>.tsv` containing the matching sample metadata.

The output rows preserve the sample order from the file manifest and include the
file manifest context (`file_id` and `file_name`) alongside the sample metadata.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


DATA_DIR = Path(__file__).resolve().parent.parent / "data"

DEFAULT_FILE_MANIFEST = DATA_DIR / "file_manifest_nemo_archive__biccn__10x_3prime_v2.tsv"
DEFAULT_SAMPLE_MANIFEST = DATA_DIR / "sample_manifest_nemo_archive__biccn__10x_3prime_v2.tsv"


@dataclass(frozen=True)
class FileRecord:
    file_id: str
    file_name: str
    sample_ids: list[str]


def normalize_file_id(file_id: str) -> str:
    return file_id.split(":", 1)[-1]


def _split_ids(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def load_file_records(file_manifest_path: str | Path) -> list[FileRecord]:
    records: list[FileRecord] = []

    with Path(file_manifest_path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            file_id = (row.get("file_id") or "").strip()
            file_name = (row.get("file_name") or "").strip()
            sample_ids = _split_ids(row.get("sample_id") or "")

            if not file_id:
                raise ValueError("File manifest row is missing file_id")
            if not file_name:
                raise ValueError(f"File manifest row for {file_id} is missing file_name")
            if not sample_ids:
                raise ValueError(f"File manifest row for {file_id} is missing sample_id values")

            records.append(FileRecord(file_id=file_id, file_name=file_name, sample_ids=sample_ids))

    return records


def load_sample_manifest(sample_manifest_path: str | Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    with Path(sample_manifest_path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        sample_rows: dict[str, dict[str, str]] = {}

        for row in reader:
            sample_id = (row.get("sample_id") or "").strip()
            if not sample_id:
                continue
            sample_rows[sample_id] = row

    if not fieldnames:
        raise ValueError(f"Sample manifest {sample_manifest_path} has no header row")

    return fieldnames, sample_rows


def build_output_rows(
    file_record: FileRecord,
    sample_fieldnames: Iterable[str],
    sample_rows: dict[str, dict[str, str]],
) -> tuple[list[str], list[dict[str, str]]]:
    output_fieldnames = ["file_id", "file_name", *sample_fieldnames]
    output_rows: list[dict[str, str]] = []

    for sample_id in file_record.sample_ids:
        sample_row = sample_rows.get(sample_id)
        if sample_row is None:
            raise KeyError(
                f"Sample id {sample_id!r} from {file_record.file_id} was not found in the sample manifest"
            )

        merged_row = {"file_id": file_record.file_id, "file_name": file_record.file_name}
        merged_row.update(sample_row)
        output_rows.append(merged_row)

    return output_fieldnames, output_rows


def write_file_sample_tsv(
    file_record: FileRecord,
    sample_fieldnames: Iterable[str],
    sample_rows: dict[str, dict[str, str]],
    output_root: str | Path = DATA_DIR,
    overwrite: bool = False,
) -> Path:
    output_fieldnames, output_rows = build_output_rows(file_record, sample_fieldnames, sample_rows)
    output_path = Path(output_root) / f"{normalize_file_id(file_record.file_id)}.tsv"

    if output_path.exists() and not overwrite:
        return output_path

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(output_rows)

    return output_path


def write_all_file_sample_tables(
    file_manifest_path: str | Path = DEFAULT_FILE_MANIFEST,
    sample_manifest_path: str | Path = DEFAULT_SAMPLE_MANIFEST,
    output_root: str | Path = DATA_DIR,
    overwrite: bool = False,
    file_id: str | None = None,
) -> list[Path]:
    file_records = load_file_records(file_manifest_path)
    sample_fieldnames, sample_rows = load_sample_manifest(sample_manifest_path)

    if file_id is not None:
        normalized_id = normalize_file_id(file_id)
        candidate_ids = {file_id, normalized_id, f"nemo:{normalized_id}"}
        file_records = [
            record for record in file_records if record.file_id in candidate_ids
        ]
        if not file_records:
            raise KeyError(f"file_id {file_id!r} was not found in {file_manifest_path}")

    output_paths: list[Path] = []
    for file_record in file_records:
        output_paths.append(
            write_file_sample_tsv(
                file_record,
                sample_fieldnames,
                sample_rows,
                output_root=output_root,
                overwrite=overwrite,
            )
        )

    return output_paths


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Write one TSV per file_id with metadata for the samples referenced by that file."
    )
    parser.add_argument(
        "--file-manifest",
        default=str(DEFAULT_FILE_MANIFEST),
        help="Path to the file manifest TSV",
    )
    parser.add_argument(
        "--sample-manifest",
        default=str(DEFAULT_SAMPLE_MANIFEST),
        help="Path to the sample manifest TSV",
    )
    parser.add_argument(
        "--output-root",
        default=str(DATA_DIR),
        help="Directory where the per-file TSVs will be written",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Rewrite existing output files",
    )
    parser.add_argument(
        "file_id",
        nargs="?",
        default=None,
        help="Only write the TSV for this file id (e.g. nemo:fil-03vd8jx). "
        "If omitted, writes a TSV for every row in the file manifest.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        output_paths = write_all_file_sample_tables(
            file_manifest_path=args.file_manifest,
            sample_manifest_path=args.sample_manifest,
            output_root=args.output_root,
            overwrite=args.overwrite,
            file_id=args.file_id,
        )
    except (KeyError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    for output_path in output_paths:
        print(output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())