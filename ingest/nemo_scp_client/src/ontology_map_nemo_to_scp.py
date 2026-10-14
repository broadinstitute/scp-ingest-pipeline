"""Annotate NEMO per-file TSVs with UBERON identifiers.

This module downloads the minified ontology assets from the SCP ingest pipeline
repository, decompresses them into a local cache, and then uses the UBERON
ontology to map each row's `sample_anatomies` value to the best UBERON ID and
preferred label. It also maps `sample_technique` to the best EFO
`library_preparation_protocol` ID and preferred label.

The default workflow reads a file such as `fil-0igba00.tsv` and writes a new
TSV named `fil-0igba00.scp.tsv` with the extra obs columns required by the SCP
ingest pipeline: `organ`, `organ__ontology_label`, `library_preparation_protocol`,
`library_preparation_protocol__ontology_label`, `species`,
`species__ontology_label`, `disease`, `disease__ontology_label`, `sex`,
`biosample_id`, and `donor_id`.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import shutil
import sys
from dataclasses import dataclass
from difflib import SequenceMatcher
from pathlib import Path
from typing import Iterable
from urllib.request import urlopen


ONTOLOGY_BASE_URL = (
    "https://raw.githubusercontent.com/broadinstitute/scp-ingest-pipeline/"
    "development/ingest/validation/ontologies"
)

ONTOLOGY_FILES = {
    "cl": "cl.min.tsv.gz",
    "efo": "efo.min.tsv.gz",
    "hancestro": "hancestro.min.tsv.gz",
    "mondo": "mondo.min.tsv.gz",
    "ncbitaxon": "ncbitaxon.min.tsv.gz",
    "pato": "pato.min.tsv.gz",
    "uberon": "uberon.min.tsv.gz",
    "uo": "uo.min.tsv.gz",
}

DATA_DIR = Path(__file__).resolve().parent.parent / "data"

DEFAULT_CACHE_DIR = DATA_DIR / ".ontology_cache"

# Some NeMO `sample_anatomies` values name a cell type rather than an
# anatomical structure, so they have no confident UBERON match. Fall back to
# the tissue that cell type is drawn from, keyed by normalized text.
ANATOMY_FALLBACKS = {
    "monocyte": "blood",
}

# Columns appended to a crosstab TSV by annotate_file_tsv(), in the order written.
# Names match the SCP/AnnData obs keys required by the ingest pipeline.
ONTOLOGY_COLUMNS = [
    "organ",
    "organ__ontology_label",
    "library_preparation_protocol",
    "library_preparation_protocol__ontology_label",
    "species",
    "species__ontology_label",
    "disease",
    "disease__ontology_label",
    "sex",
    "biosample_id",
    "donor_id",
]


@dataclass(frozen=True)
class OntologyTerm:
    ontology_id: str
    label: str
    synonyms: tuple[str, ...]


class OntologyMatcher:
    def __init__(self, terms: Iterable[OntologyTerm]):
        self.terms = list(terms)
        self.label_index: dict[str, list[OntologyTerm]] = {}
        self.synonym_index: dict[str, list[OntologyTerm]] = {}

        for term in self.terms:
            self.label_index.setdefault(_normalize_text(term.label), []).append(term)
            for synonym in term.synonyms:
                self.synonym_index.setdefault(_normalize_text(synonym), []).append(term)

    def best_match(self, sample_anatomy: str) -> OntologyTerm | None:
        sample_anatomy = (sample_anatomy or "").strip()
        if not sample_anatomy:
            return None

        normalized_query = _normalize_text(sample_anatomy)

        exact_matches = self.label_index.get(normalized_query)
        if exact_matches:
            return sorted(exact_matches, key=lambda term: (len(term.label), term.ontology_id))[0]

        exact_matches = self.synonym_index.get(normalized_query)
        if exact_matches:
            return sorted(exact_matches, key=lambda term: (len(term.label), term.ontology_id))[0]

        best_term: OntologyTerm | None = None
        best_ratio = 0.0
        # Among equal ratios, prefer the shortest label and smallest ontology
        # id, matching the exact-match tie-breaks above.
        best_tie_break: tuple[int, str] = (10**9, "")

        for term in self.terms:
            ratio = _score_candidate(sample_anatomy, term)
            tie_break = (len(term.label), term.ontology_id)
            if ratio > best_ratio or (ratio == best_ratio and tie_break < best_tie_break):
                best_ratio = ratio
                best_tie_break = tie_break
                best_term = term

        if best_term is None or best_ratio < 0.75:
            return None

        return best_term


def _ontology_url(filename: str) -> str:
    return f"{ONTOLOGY_BASE_URL}/{filename}"


def _normalize_text(value: str) -> str:
    value = value.strip().lower()
    value = re.sub(r"[^a-z0-9]+", " ", value)
    return re.sub(r"\s+", " ", value).strip()


def _split_synonyms(value: str) -> tuple[str, ...]:
    if not value:
        return ()
    return tuple(item.strip() for item in value.split("||") if item.strip())


def ensure_minified_ontology(
    ontology_name: str,
    cache_dir: str | Path = DEFAULT_CACHE_DIR,
) -> Path:
    """Download and decompress a minified ontology TSV into the cache."""

    if ontology_name not in ONTOLOGY_FILES:
        raise KeyError(f"Unknown ontology {ontology_name!r}")

    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    gz_name = ONTOLOGY_FILES[ontology_name]
    tsv_name = gz_name.removesuffix(".gz")
    tsv_path = cache_dir / tsv_name

    if tsv_path.exists():
        return tsv_path

    gz_path = cache_dir / gz_name
    url = _ontology_url(gz_name)

    with urlopen(url) as response, gz_path.open("wb") as handle:
        shutil.copyfileobj(response, handle)

    with gzip.open(gz_path, "rb") as gz_handle, tsv_path.open("wb") as out_handle:
        shutil.copyfileobj(gz_handle, out_handle)

    gz_path.unlink(missing_ok=True)
    return tsv_path


def load_ontology_terms(ontology_name: str, cache_dir: str | Path = DEFAULT_CACHE_DIR) -> list[OntologyTerm]:
    tsv_path = ensure_minified_ontology(ontology_name, cache_dir=cache_dir)
    terms: list[OntologyTerm] = []

    with tsv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            ontology_id = row[0].strip()
            label = row[1].strip() if len(row) > 1 else ""
            synonyms = _split_synonyms(row[2].strip() if len(row) > 2 else "")
            if ontology_id and label:
                terms.append(OntologyTerm(ontology_id=ontology_id, label=label, synonyms=synonyms))

    return terms


def _score_candidate(query: str, term: OntologyTerm) -> float:
    """Return the best fuzzy similarity ratio between the query and the term."""

    normalized_query = _normalize_text(query)

    best_ratio = SequenceMatcher(None, normalized_query, _normalize_text(term.label)).ratio()
    for synonym in term.synonyms:
        best_ratio = max(best_ratio, SequenceMatcher(None, normalized_query, _normalize_text(synonym)).ratio())

    return best_ratio


def annotate_file_tsv(
    file_tsv_path: str | Path,
    cache_dir: str | Path = DEFAULT_CACHE_DIR,
    output_path: str | Path | None = None,
    overwrite: bool = False,
) -> Path:
    """Annotate a per-file TSV with a UBERON ID column."""

    file_tsv_path = Path(file_tsv_path)
    output_path = Path(output_path) if output_path is not None else file_tsv_path.with_name(
        f"{file_tsv_path.stem}.scp.tsv"
    )

    if output_path.exists() and not overwrite:
        return output_path

    uberon_matcher = OntologyMatcher(load_ontology_terms("uberon", cache_dir=cache_dir))

    # Pre-scan to find the biosample name from the nemo:smp- row.
    biosample_name = ""
    with file_tsv_path.open(newline="", encoding="utf-8") as pre_handle:
        for row in csv.DictReader(pre_handle, delimiter="\t"):
            if (row.get("sample_id", "") or "").startswith("nemo:smp-"):
                biosample_name = row.get("sample_name", "") or ""
                break

    with file_tsv_path.open(newline="", encoding="utf-8") as input_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{file_tsv_path} has no header row")

        if "sample_anatomies" not in reader.fieldnames:
            raise ValueError(f"{file_tsv_path} does not contain a sample_anatomies column")
        if "sample_technique" not in reader.fieldnames:
            raise ValueError(f"{file_tsv_path} does not contain a sample_technique column")
        if "taxon" not in reader.fieldnames:
            raise ValueError(f"{file_tsv_path} does not contain a taxon column")

        output_fieldnames = [*reader.fieldnames, *ONTOLOGY_COLUMNS]
        output_path.parent.mkdir(parents=True, exist_ok=True)

        efo_matcher = OntologyMatcher(load_ontology_terms("efo", cache_dir=cache_dir))
        ncbitaxon_matcher = OntologyMatcher(load_ontology_terms("ncbitaxon", cache_dir=cache_dir))
        pato_matcher = OntologyMatcher(load_ontology_terms("pato", cache_dir=cache_dir))

        # NeMO crosstabs carry no disease annotation, so every sample is treated
        # as unaffected/normal per the PATO convention for a missing disease state.
        normal_term = pato_matcher.best_match("normal")

        with output_path.open("w", newline="", encoding="utf-8") as output_handle:
            writer = csv.DictWriter(output_handle, fieldnames=output_fieldnames, delimiter="\t")
            writer.writeheader()

            for row in reader:
                sample_anatomy = row.get("sample_anatomies", "")
                if sample_anatomy == "None":
                    sample_anatomy = ""
                best_term = uberon_matcher.best_match(sample_anatomy)
                if best_term is None:
                    fallback_label = ANATOMY_FALLBACKS.get(_normalize_text(sample_anatomy))
                    if fallback_label:
                        best_term = uberon_matcher.best_match(fallback_label)
                row["organ"] = best_term.ontology_id if best_term else ""
                row["organ__ontology_label"] = best_term.label if best_term else ""

                sample_technique = row.get("sample_technique", "")
                if sample_technique == "None":
                    sample_technique = ""
                technique_term = efo_matcher.best_match(sample_technique)
                row["library_preparation_protocol"] = technique_term.ontology_id if technique_term else ""
                row["library_preparation_protocol__ontology_label"] = technique_term.label if technique_term else ""

                taxon = row.get("taxon", "")
                if taxon == "None":
                    taxon = ""
                taxon_term = ncbitaxon_matcher.best_match(taxon)
                row["species"] = taxon_term.ontology_id if taxon_term else ""
                row["species__ontology_label"] = taxon_term.label if taxon_term else ""

                row["disease"] = normal_term.ontology_id if normal_term else ""
                row["disease__ontology_label"] = normal_term.label if normal_term else ""

                subject_sex = row.get("subject_sex", "")
                row["sex"] = subject_sex if subject_sex != "None" else ""

                row["biosample_id"] = biosample_name
                row["donor_id"] = biosample_name
                writer.writerow(row)

    return output_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Annotate per-file TSVs with the best UBERON ID for sample_anatomies."
    )
    parser.add_argument(
        "file_tsv",
        nargs="+",
        help="One or more fil-<ID>.tsv filenames to annotate, resolved against "
        "--data-dir unless given as an absolute or otherwise-existing path",
    )
    parser.add_argument(
        "--data-dir",
        default=str(DATA_DIR),
        help="Directory to resolve bare file_tsv filenames against",
    )
    parser.add_argument(
        "--cache-dir",
        default=str(DEFAULT_CACHE_DIR),
        help="Directory used to cache downloaded ontology TSVs",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Optional directory for the annotated TSV outputs",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Rewrite existing annotated TSVs",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        for file_tsv in args.file_tsv:
            input_path = Path(file_tsv)
            if not input_path.exists() and not input_path.is_absolute():
                input_path = Path(args.data_dir) / input_path
            output_path = None
            if args.output_dir:
                output_path = Path(args.output_dir) / f"{input_path.stem}.scp.tsv"

            annotated_path = annotate_file_tsv(
                input_path,
                cache_dir=args.cache_dir,
                output_path=output_path,
                overwrite=args.overwrite,
            )
            print(annotated_path)
    except (KeyError, ValueError, OSError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())