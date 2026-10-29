#!/usr/bin/env python3
# Author: Victoria Begley (begley@ebi.ac.uk)
# License: MIT
from __future__ import annotations

import argparse
import re
import sys
import urllib.error
import urllib.request


class EnsemblSpeciesNotFound(Exception):
    pass


_EG_FLAT_DIVISIONS = [
    ("metazoa",  "https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa"),
    ("fungi",    "https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi"),
    ("plants",   "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants"),
    ("protists", "https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists"),
]
_EG_BACTERIA_BASE = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download one Ensembl GTF annotation file for one species."
    )
    parser.add_argument("--species", required=True)
    parser.add_argument("--release", required=True)
    parser.add_argument("--base-url", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--source-urls", required=True)
    parser.add_argument(
        "--not-found-file",
        required=True,
        help="Path to write (empty) when the species is absent from all Ensembl FTPs; "
             "process exits 0 and GTF is silently skipped for this reference.",
    )
    return parser.parse_args()


def fetch_text(url: str, timeout: int = 60) -> str:
    try:
        with urllib.request.urlopen(url, timeout=timeout) as response:
            return response.read().decode("utf-8", errors="ignore")
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            raise EnsemblSpeciesNotFound(f"HTTP 404 at {url}") from exc
        fallback_url = f"{url}index.html" if url.endswith("/") else f"{url}/index.html"
        try:
            with urllib.request.urlopen(fallback_url, timeout=timeout) as response:
                return response.read().decode("utf-8", errors="ignore")
        except urllib.error.HTTPError as fallback_exc:
            if fallback_exc.code == 404:
                raise EnsemblSpeciesNotFound(f"HTTP 404 at {fallback_url}") from fallback_exc
            raise
    except urllib.error.URLError:
        fallback_url = f"{url}index.html" if url.endswith("/") else f"{url}/index.html"
        try:
            with urllib.request.urlopen(fallback_url, timeout=timeout) as response:
                return response.read().decode("utf-8", errors="ignore")
        except urllib.error.HTTPError as fallback_exc:
            if fallback_exc.code == 404:
                raise EnsemblSpeciesNotFound(f"HTTP 404 at {fallback_url}") from fallback_exc
            raise urllib.error.URLError(str(fallback_exc)) from fallback_exc


def release_path_for_value(release: str) -> str:
    if release in ("current", "latest"):
        return "current/gtf"
    if release.startswith("release-"):
        return f"{release}/gtf"
    return f"release-{release}/gtf"


def _find_gtf_in_listing(species_root: str, timeout: int = 60) -> str:
    listing = fetch_text(species_root, timeout=timeout)
    matches = [m for m in re.findall(r'href="([^"]+)"', listing) if re.search(r"\.gtf\.gz$", m)]
    non_abinitio = [m for m in matches if "abinitio" not in m.lower()]
    # Prefer chromosome-only GTF (matches primary_assembly FASTA, excludes patches/haplotypes),
    # then fall back to the full GTF.
    chr_only = [m for m in non_abinitio if re.search(r"\.chr\.gtf\.gz$", m)]
    gtf_name = (chr_only or non_abinitio or matches or [None])[0]
    if not gtf_name:
        raise RuntimeError(f"No .gtf.gz file found at {species_root}")
    return gtf_name


def _try_flat_division(base_url: str, release: str, species: str) -> tuple[str, str]:
    """Return (species_root_url, gtf_name) for a flat Ensembl/EnsemblGenomes division."""
    species_root = f"{base_url}/{release_path_for_value(release)}/{species}/"
    gtf_name = _find_gtf_in_listing(species_root)
    return species_root, gtf_name


def _try_bacteria(release: str, species: str) -> tuple[str, str]:
    """Search EnsemblBacteria collections for *species*, return (species_root_url, gtf_name)."""
    release_path = release_path_for_value(release)
    gtf_root = f"{_EG_BACTERIA_BASE}/{release_path}"

    try:
        listing = fetch_text(f"{gtf_root}/", timeout=30)
    except Exception as exc:
        raise EnsemblSpeciesNotFound(f"Cannot access EnsemblBacteria GTF FTP: {exc}") from exc

    collections = sorted(re.findall(r'href="(bacteria_\d+_collection/)"', listing))
    if not collections:
        raise EnsemblSpeciesNotFound(
            "No bacteria_N_collection directories found at EnsemblBacteria GTF FTP"
        )

    for collection in collections:
        collection_name = collection.rstrip("/")
        species_root = f"{gtf_root}/{collection_name}/{species}/"
        try:
            gtf_name = _find_gtf_in_listing(species_root, timeout=10)
            print(
                f"[ENSEMBL_GTF] Found '{species}' GTF in EnsemblBacteria "
                f"collection '{collection_name}'.",
                file=sys.stderr,
            )
            return species_root, gtf_name
        except (EnsemblSpeciesNotFound, RuntimeError):
            continue

    raise EnsemblSpeciesNotFound(
        f"Species '{species}' GTF not found in any EnsemblBacteria collection"
    )


def _find_species(base_url: str, release: str, species: str) -> tuple[str, str, str]:
    """Try every Ensembl source in priority order.

    Returns (species_root_url, gtf_name, division_label).
    """
    try:
        species_root, gtf_name = _try_flat_division(base_url, release, species)
        return species_root, gtf_name, "Ensembl"
    except EnsemblSpeciesNotFound:
        pass

    for division_name, division_base in _EG_FLAT_DIVISIONS:
        try:
            species_root, gtf_name = _try_flat_division(division_base, release, species)
            return species_root, gtf_name, f"EnsemblGenomes/{division_name}"
        except EnsemblSpeciesNotFound:
            continue

    try:
        species_root, gtf_name = _try_bacteria(release, species)
        return species_root, gtf_name, "EnsemblBacteria"
    except EnsemblSpeciesNotFound:
        pass

    raise EnsemblSpeciesNotFound(
        f"Species '{species}' not found on Ensembl, EnsemblGenomes, or EnsemblBacteria FTP"
    )


def main() -> int:
    args = parse_args()
    species = args.species.strip().lower().replace(" ", "_")
    release = args.release.strip()
    base_url = args.base_url.rstrip("/")

    try:
        species_root, gtf_name, division = _find_species(base_url, release, species)
    except EnsemblSpeciesNotFound as exc:
        print(
            f"[ENSEMBL_GTF] Species '{species}' not found on any Ensembl FTP "
            f"— dotplot-to-bp conversion will be skipped: {exc}",
            file=sys.stderr,
        )
        with open(args.not_found_file, "w", encoding="utf-8"):
            pass
        return 0
    except (urllib.error.URLError, OSError) as exc:
        print(
            f"[ENSEMBL_GTF] Network error reaching Ensembl FTP for '{species}' "
            f"— falling back to NCBI if available: {exc}",
            file=sys.stderr,
        )
        with open(args.not_found_file, "w", encoding="utf-8"):
            pass
        return 0

    print(
        f"[ENSEMBL_GTF] Resolved '{species}' GTF via {division}.",
        file=sys.stderr,
    )

    gtf_url = f"{species_root}{gtf_name}"
    urllib.request.urlretrieve(gtf_url, args.output)

    with open(args.source_urls, "w", encoding="utf-8") as handle:
        handle.write(f"{gtf_url}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
