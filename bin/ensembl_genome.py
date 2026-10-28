#!/usr/bin/env python3
# Author: Victoria Begley (begley@ebi.ac.uk)
# License: MIT
"""Download soft-masked genome FASTA from Ensembl FTP for one species."""
from __future__ import annotations

import argparse
import gzip
import io
import shutil
import sys
import urllib.error
import urllib.request
import re


class EnsemblSpeciesNotFound(Exception):
    pass


_EG_FLAT_DIVISIONS = [
    ("metazoa",  "https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa"),
    ("fungi",    "https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi"),
    ("plants",   "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants"),
    ("protists", "https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists"),
]

# Preference order: primary assembly (no alt/patch contigs, reduces multi-mapping) over toplevel,
# and soft-masked over unmasked (prevents spurious alignments to repetitive regions).
_GENOME_PATTERNS = [
    r"\.dna_sm\.primary_assembly\.fa\.gz$",
    r"\.dna\.primary_assembly\.fa\.gz$",
    r"\.dna_sm\.toplevel\.fa\.gz$",
    r"\.dna\.toplevel\.fa\.gz$",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download genome FASTA from Ensembl for one species."
    )
    parser.add_argument("--species",        required=True)
    parser.add_argument("--release",        required=True)
    parser.add_argument("--base-url",       required=True)
    parser.add_argument("--output",         required=True, help="Output path for gzipped genome FASTA")
    parser.add_argument("--source-url",     required=True, help="Output file recording the download URL")
    parser.add_argument("--not-found-file", required=True,
                        help="Touched (empty) when the species is absent from all Ensembl FTPs; "
                             "process exits 0 and a fallback is triggered by the pipeline.")
    parser.add_argument("--base-url-override", default=None, help=argparse.SUPPRESS)
    return parser.parse_args()


def fetch_text(url: str, timeout: int = 60) -> str:
    try:
        with urllib.request.urlopen(url, timeout=timeout) as response:
            return response.read().decode("utf-8", errors="ignore")
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            raise EnsemblSpeciesNotFound(f"HTTP 404 at {url}") from exc
        # Some servers (and file:// dir URLs) don't serve directory listings directly; retry index.html.
        fallback_url = f"{url}index.html" if url.endswith("/") else f"{url}/index.html"
        try:
            with urllib.request.urlopen(fallback_url, timeout=timeout) as response:
                return response.read().decode("utf-8", errors="ignore")
        except urllib.error.HTTPError as fallback_exc:
            if fallback_exc.code == 404:
                raise EnsemblSpeciesNotFound(f"HTTP 404 at {fallback_url}") from fallback_exc
            raise
    except urllib.error.URLError:
        # file:// directory URLs raise URLError ("Is a directory"); fall back to the index.html listing.
        fallback_url = f"{url}index.html" if url.endswith("/") else f"{url}/index.html"
        try:
            with urllib.request.urlopen(fallback_url, timeout=timeout) as response:
                return response.read().decode("utf-8", errors="ignore")
        except urllib.error.HTTPError as fallback_exc:
            if fallback_exc.code == 404:
                raise EnsemblSpeciesNotFound(f"HTTP 404 at {fallback_url}") from fallback_exc
            raise urllib.error.URLError(str(fallback_exc)) from fallback_exc


def release_path_for_value(release: str, base_url: str) -> str:
    # The "current" alias differs by site: main Ensembl serves it at current_fasta/,
    # while EnsemblGenomes divisions (metazoa/fungi/plants/protists) serve it at current/fasta/.
    if release in ("current", "latest"):
        return "current/fasta" if "ensemblgenomes" in base_url else "current_fasta"
    if release.startswith("release-"):
        return f"{release}/fasta"
    return f"release-{release}/fasta"


def _find_genome_in_dna_dir(dna_dir: str) -> str | None:
    """Return the genome FASTA URL inside a species dna/ directory, or None."""
    try:
        listing = fetch_text(dna_dir, timeout=30)
    except (EnsemblSpeciesNotFound, urllib.error.URLError, urllib.error.HTTPError):
        return None
    hrefs = re.findall(r'href="([^"]+)"', listing)
    for pattern in _GENOME_PATTERNS:
        for href in hrefs:
            filename = href.split("/")[-1].strip()
            if re.search(pattern, filename):
                return f"{dna_dir}{filename}"
    return None


def _find_species_genome(base_url: str, release: str, species: str) -> tuple[str, str]:
    """Return (genome_fasta_url, division_label). Raises EnsemblSpeciesNotFound if absent."""

    # 1. Primary Ensembl (vertebrates + some others)
    release_path = release_path_for_value(release, base_url)
    url = _find_genome_in_dna_dir(f"{base_url}/{release_path}/{species}/dna/")
    if url:
        return url, "Ensembl"

    # 2. EnsemblGenomes flat divisions (metazoa, fungi, plants, protists)
    for division_name, division_base in _EG_FLAT_DIVISIONS:
        division_release_path = release_path_for_value(release, division_base)
        url = _find_genome_in_dna_dir(f"{division_base}/{division_release_path}/{species}/dna/")
        if url:
            return url, f"EnsemblGenomes/{division_name}"

    # 3. Prefix fallback — handles strain/subspecies suffixes (e.g. 'mus_musculus_c57bl6j')
    prefix = f"{species}_"
    all_sources = [("Ensembl", base_url)] + [
        (f"EnsemblGenomes/{n}", b) for n, b in _EG_FLAT_DIVISIONS
    ]
    for label, src_url in all_sources:
        src_release_path = release_path_for_value(release, src_url)
        root = f"{src_url}/{src_release_path}/"
        try:
            listing = fetch_text(root, timeout=30)
        except Exception:
            continue
        all_dirs = re.findall(r'href="([a-z][a-z0-9_]+)/"', listing)
        for d in sorted(d for d in all_dirs if d == species or d.startswith(prefix)):
            url = _find_genome_in_dna_dir(f"{src_url}/{src_release_path}/{d}/dna/")
            if url:
                print(
                    f"[ENSEMBL_GENOME] '{species}' matched '{d}' on {label}.",
                    file=sys.stderr,
                )
                return url, label

    raise EnsemblSpeciesNotFound(
        f"Species '{species}' genome not found on Ensembl or EnsemblGenomes"
    )


def main() -> int:
    args = parse_args()
    species  = args.species.strip().lower().replace(" ", "_")
    release  = args.release.strip()
    base_url = (args.base_url_override or args.base_url).rstrip("/")

    try:
        genome_url, division = _find_species_genome(base_url, release, species)
    except EnsemblSpeciesNotFound as exc:
        print(
            f"[ENSEMBL_GENOME] Species '{species}' not found on any Ensembl FTP "
            f"— will fall back: {exc}",
            file=sys.stderr,
        )
        with open(args.not_found_file, "w", encoding="utf-8"):
            pass
        return 0

    print(
        f"[ENSEMBL_GENOME] Resolved '{species}' via {division}: {genome_url}",
        file=sys.stderr,
    )

    with urllib.request.urlopen(genome_url, timeout=600) as response:
        with gzip.GzipFile(fileobj=io.BufferedReader(response)) as gz_in:
            with open(args.output, "wb") as fa_out:
                shutil.copyfileobj(gz_in, fa_out, length=1024 * 1024)

    with open(args.source_url, "w", encoding="utf-8") as fh:
        fh.write(f"{genome_url}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
