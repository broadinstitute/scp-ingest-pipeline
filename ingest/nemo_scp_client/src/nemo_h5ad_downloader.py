"""Download H5AD files from a NEMO file manifest.

The manifest is expected to have at least these columns:
- file_id
- file_name
- urls

Example:
    download_h5ad("nemo:fil-03vd8jx")

This will create (relative to this module's ../data directory):
    data/fil-03vd8jx/UG58XXVD_112T_FrW.h5ad
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import os
import select
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
import re
import time
from typing import Iterable
from urllib.error import URLError
from urllib.request import urlopen
import socket


DATA_DIR = Path(__file__).resolve().parent.parent / "data"

DEFAULT_MANIFEST = DATA_DIR / "file_manifest_nemo_archive__biccn__10x_3prime_v2.tsv"


@dataclass(frozen=True)
class ManifestRecord:
    file_id: str
    file_name: str
    urls: list[str]


def normalize_file_id(file_id: str) -> str:
    """Return the folder-friendly id used on disk."""

    return file_id.split(":", 1)[-1]


def _candidate_file_ids(file_id: str) -> list[str]:
    normalized = normalize_file_id(file_id)
    if file_id == normalized:
        return [normalized, f"nemo:{normalized}"]
    return [file_id, normalized]


def _split_urls(urls_field: str) -> list[str]:
    return [url.strip() for url in urls_field.split(",") if url.strip()]


def load_manifest_record(manifest_path: str | Path, file_id: str) -> ManifestRecord:
    """Find a manifest row by file_id."""

    candidate_ids = set(_candidate_file_ids(file_id))

    with Path(manifest_path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            row_file_id = (row.get("file_id") or "").strip()
            if row_file_id in candidate_ids:
                file_name = (row.get("file_name") or "").strip()
                urls = _split_urls(row.get("urls") or "")
                if not file_name:
                    raise ValueError(f"Manifest row for {file_id} is missing file_name")
                if not urls:
                    raise ValueError(f"Manifest row for {file_id} is missing urls")
                return ManifestRecord(file_id=row_file_id, file_name=file_name, urls=urls)

    raise KeyError(f"file_id {file_id!r} was not found in {manifest_path}")


def pick_download_url(urls: Iterable[str]) -> str:
    """Prefer a GS URL over HTTPS if multiple sources are present."""

    urls = list(urls)
    if not urls:
        raise ValueError("No URLs were provided")

    for url in urls:
        if url.startswith("gs://"):
            return url

    for url in urls:
        if url.startswith("https://"):
            return url

    for url in urls:
        if url.startswith("http://"):
            return url

    return urls[0]


def _format_bytes(byte_count: int) -> str:
    units = ["B", "KB", "MB", "GB", "TB"]
    value = float(byte_count)
    for unit in units:
        if value < 1024.0 or unit == units[-1]:
            if unit == "B":
                return f"{int(value)} {unit}"
            return f"{value:.1f} {unit}"
        value /= 1024.0
    return f"{byte_count} B"


def _format_minutes(seconds: float | None) -> str:
    if seconds is None or seconds < 0:
        return "unknown"
    return f"{seconds / 60.0:.1f}"


def _resolve_gcloud_billing_project(explicit: str | None = None) -> str:
    if explicit:
        return explicit

    env_project = os.environ.get("GCLOUD_BILLING_PROJECT")
    if env_project:
        return env_project

    result = _run_gcloud(["gcloud", "config", "get-value", "project"])
    project = result.stdout.strip()
    if project and project != "(unset)":
        return project

    raise RuntimeError(
        "A billing project is required for gs:// downloads; pass --billing-project or set GCLOUD_BILLING_PROJECT"
    )


def _run_gcloud(command: list[str]) -> subprocess.CompletedProcess[str]:
    """Run a gcloud command, converting failures into RuntimeError for the CLI."""

    try:
        return subprocess.run(command, check=True, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            "gcloud is required for gs:// downloads but was not found on PATH"
        ) from exc
    except subprocess.CalledProcessError as exc:
        stderr_output = (exc.stderr or "").strip()
        raise RuntimeError(
            f"{' '.join(command)} failed with exit code {exc.returncode}: "
            f"{stderr_output or 'no stderr output'}"
        ) from exc


def _gcloud_object_size(url: str, billing_project: str) -> int | None:
    result = _run_gcloud(
        [
            "gcloud",
            "storage",
            "ls",
            "-L",
            f"--billing-project={billing_project}",
            url,
        ]
    )

    match = re.search(r"Content-Length:\s+(\d+)", result.stdout)
    if match:
        return int(match.group(1))

    return None


def download_gs_url(
    url: str,
    destination_path: Path,
    billing_project: str,
    timeout_seconds: int = 60,
) -> Path:
    """Download a gs:// URL with gcloud while preserving local progress reporting.

    Fails with RuntimeError if no data arrives for `timeout_seconds`, matching
    the stall behavior of the HTTPS path. gcloud's stderr goes to a temp file
    rather than a pipe, so a chatty transfer can never deadlock on a full pipe
    buffer.
    """

    destination_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = destination_path.with_suffix(destination_path.suffix + ".part")
    started_at = time.monotonic()
    bytes_downloaded = 0
    total_bytes_int = _gcloud_object_size(url, billing_project)

    command = [
        "gcloud",
        "storage",
        "cat",
        f"--billing-project={billing_project}",
        url,
    ]

    with tempfile.TemporaryFile() as stderr_handle:
        try:
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=stderr_handle,
                bufsize=0,
            )
        except FileNotFoundError as exc:
            raise RuntimeError(
                "gcloud is required for gs:// downloads but was not found on PATH"
            ) from exc

        try:
            assert process.stdout is not None

            _print_progress(
                downloaded=0,
                total=total_bytes_int,
                started_at=started_at,
                prefix=f"downloading {destination_path.name}",
            )

            with temp_path.open("wb") as handle:
                while True:
                    ready, _, _ = select.select([process.stdout], [], [], timeout_seconds)
                    if not ready:
                        raise RuntimeError(
                            f"Download stalled after "
                            f"{_format_minutes(time.monotonic() - started_at)} min: {url}"
                        )
                    chunk = process.stdout.read(1024 * 1024)
                    if not chunk:
                        break
                    handle.write(chunk)
                    bytes_downloaded += len(chunk)
                    _print_progress(
                        downloaded=bytes_downloaded,
                        total=total_bytes_int,
                        started_at=started_at,
                        prefix=f"downloading {destination_path.name}",
                    )

            return_code = process.wait()
            stderr_handle.seek(0)
            stderr_output = stderr_handle.read().decode("utf-8", errors="replace").strip()
            if return_code != 0:
                raise RuntimeError(
                    f"gcloud storage cat failed with exit code {return_code}: {stderr_output or 'no stderr output'}"
                )

        except Exception:
            if process.poll() is None:
                process.kill()
                process.wait()
            if temp_path.exists():
                temp_path.unlink()
            print(file=sys.stderr)
            raise

    temp_path.replace(destination_path)
    elapsed_seconds = max(time.monotonic() - started_at, 0.0)
    print(
        f"\rdownloaded {destination_path.name} "
        f"elapsed={elapsed_seconds / 60.0:.1f} min "
        f"downloaded={_format_bytes(bytes_downloaded)} "
        f"remaining=0 B eta=0.0 min time={datetime.now().strftime('%H:%M:%S')}"
        f"{' ' * 12}",
        file=sys.stderr,
        end="\n",
        flush=True,
    )
    return destination_path


def _print_progress(
    *,
    downloaded: int,
    total: int | None,
    started_at: float,
    prefix: str,
) -> None:
    elapsed_seconds = max(time.monotonic() - started_at, 0.0)
    clock_time = datetime.now().strftime("%H:%M:%S")

    if total is None or total <= 0:
        remaining_text = "unknown"
        eta_text = "unknown"
        total_text = "unknown"
    else:
        remaining_bytes = max(total - downloaded, 0)
        total_text = _format_bytes(total)
        remaining_text = _format_bytes(remaining_bytes)
        if downloaded > 0 and elapsed_seconds > 0:
            rate = downloaded / elapsed_seconds
            eta_seconds = remaining_bytes / rate if rate > 0 else None
        else:
            eta_seconds = None
        eta_text = _format_minutes(eta_seconds)

    message = (
        f"\r{prefix} "
        f"elapsed={elapsed_seconds / 60.0:.1f} min "
        f"downloaded={_format_bytes(downloaded)} "
        f"remaining={remaining_text} "
        f"eta={eta_text} min "
        f"time={clock_time}"
    )

    if total is not None and total > 0:
        message += f" total={total_text}"

    print(message, end="", file=sys.stderr, flush=True)


def download_url_with_timeout(url: str, destination_path: Path, timeout_seconds: int) -> Path:
    """Download a URL to destination_path, failing if the socket stalls too long."""

    destination_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = destination_path.with_suffix(destination_path.suffix + ".part")
    started_at = time.monotonic()
    bytes_downloaded = 0

    try:
        with urlopen(url, timeout=timeout_seconds) as response, temp_path.open("wb") as handle:
            total_bytes = response.headers.get("Content-Length")
            try:
                total_bytes_int = int(total_bytes) if total_bytes is not None else None
            except (TypeError, ValueError):
                total_bytes_int = None
            _print_progress(
                downloaded=0,
                total=total_bytes_int,
                started_at=started_at,
                prefix=f"downloading {destination_path.name}",
            )

            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                handle.write(chunk)
                bytes_downloaded += len(chunk)

                _print_progress(
                    downloaded=bytes_downloaded,
                    total=total_bytes_int,
                    started_at=started_at,
                    prefix=f"downloading {destination_path.name}",
                )
    except (TimeoutError, socket.timeout) as exc:
        print(file=sys.stderr)
        raise RuntimeError(
            f"Download stalled after {_format_minutes(time.monotonic() - started_at)} min: {url}"
        ) from exc
    except URLError as exc:
        print(file=sys.stderr)
        if temp_path.exists():
            temp_path.unlink()
        raise RuntimeError(f"Failed to download {url}: {exc}") from exc
    except Exception:
        print(file=sys.stderr)
        if temp_path.exists():
            temp_path.unlink()
        raise

    temp_path.replace(destination_path)
    elapsed_seconds = max(time.monotonic() - started_at, 0.0)
    print(
        f"\rdownloaded {destination_path.name} "
        f"elapsed={elapsed_seconds / 60.0:.1f} min "
        f"downloaded={_format_bytes(bytes_downloaded)} "
        f"remaining=0 B eta=0.0 min time={datetime.now().strftime('%H:%M:%S')}"
        f"{' ' * 12}",
        file=sys.stderr,
        end="\n",
        flush=True,
    )
    return destination_path


def download_h5ad(
    file_id: str,
    manifest_path: str | Path = DEFAULT_MANIFEST,
    output_root: str | Path = DATA_DIR,
    overwrite: bool = False,
    timeout_seconds: int = 60,
    billing_project: str | None = None,
) -> Path:
    """Download the H5AD for a file_id into a folder named after the id."""

    record = load_manifest_record(manifest_path, file_id)
    download_url_str = pick_download_url(record.urls)

    output_dir = Path(output_root) / normalize_file_id(record.file_id)
    destination_path = output_dir / record.file_name

    if destination_path.exists() and not overwrite:
        return destination_path

    if download_url_str.startswith("gs://"):
        return download_gs_url(
            download_url_str,
            destination_path,
            billing_project=_resolve_gcloud_billing_project(billing_project),
            timeout_seconds=timeout_seconds,
        )

    return download_url_with_timeout(download_url_str, destination_path, timeout_seconds)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download a NEMO H5AD file from a manifest row."
    )
    parser.add_argument("file_id", help="File id such as nemo:fil-03vd8jx")
    parser.add_argument(
        "--manifest",
        default=str(DEFAULT_MANIFEST),
        help="Path to the TSV manifest file",
    )
    parser.add_argument(
        "--output-root",
        default=str(DATA_DIR),
        help="Directory where the file_id folder will be created",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Redownload even if the destination file already exists",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=60,
        help="Socket timeout for stalled downloads",
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

    try:
        downloaded_path = download_h5ad(
            args.file_id,
            manifest_path=args.manifest,
            output_root=args.output_root,
            overwrite=args.overwrite,
            timeout_seconds=args.timeout_seconds,
            billing_project=args.billing_project,
        )
    except (KeyError, ValueError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(downloaded_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())