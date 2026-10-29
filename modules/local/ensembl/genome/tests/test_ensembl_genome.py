import importlib.util
import sys
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest

_SCRIPT = Path(__file__).resolve().parents[5] / "bin" / "ensembl_genome.py"
spec = importlib.util.spec_from_file_location("ensembl_genome", _SCRIPT)
eg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(eg)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _dir_html(*filenames: str) -> str:
    """Minimal FTP-style HTML listing."""
    links = "".join(f'<a href="{f}">{f}</a>\n' for f in filenames)
    return f"<html><body>{links}</body></html>"


def _species_listing(*species: str) -> str:
    links = "".join(f'<a href="{s}/">{s}/</a>\n' for s in species)
    return f"<html><body>{links}</body></html>"


# ---------------------------------------------------------------------------
# release_path_for_value
# ---------------------------------------------------------------------------

def test_release_path_current():
    assert eg.release_path_for_value("current") == "current_fasta"

def test_release_path_latest():
    assert eg.release_path_for_value("latest") == "current_fasta"

def test_release_path_numbered():
    assert eg.release_path_for_value("112") == "release-112/fasta"

def test_release_path_prefixed():
    assert eg.release_path_for_value("release-112") == "release-112/fasta"


# ---------------------------------------------------------------------------
# _find_genome_in_dna_dir
# ---------------------------------------------------------------------------

def test_find_genome_prefers_dna_sm():
    html = _dir_html(
        "Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz",
        "Homo_sapiens.GRCh38.dna.toplevel.fa.gz",
    )
    with patch.object(eg, "fetch_text", return_value=html):
        url = eg._find_genome_in_dna_dir("https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/")
    assert url and "dna_sm" in url


def test_find_genome_falls_back_to_unmasked():
    html = _dir_html("Homo_sapiens.GRCh38.dna.toplevel.fa.gz")
    with patch.object(eg, "fetch_text", return_value=html):
        url = eg._find_genome_in_dna_dir("https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/")
    assert url and "dna.toplevel" in url


def test_find_genome_returns_none_when_absent():
    html = _dir_html("README.txt", "CHECKSUMS")
    with patch.object(eg, "fetch_text", return_value=html):
        url = eg._find_genome_in_dna_dir("https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/")
    assert url is None


def test_find_genome_returns_none_on_404():
    with patch.object(eg, "fetch_text", side_effect=eg.EnsemblSpeciesNotFound("404")):
        url = eg._find_genome_in_dna_dir("https://ftp.ensembl.org/pub/current_fasta/unknown/dna/")
    assert url is None


# ---------------------------------------------------------------------------
# _find_species_genome — primary Ensembl hit
# ---------------------------------------------------------------------------

def test_find_species_genome_primary_ensembl():
    genome_file = "Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz"

    def fake_fetch(url, timeout=60):
        if "/homo_sapiens/dna/" in url:
            return _dir_html(genome_file)
        raise eg.EnsemblSpeciesNotFound(f"404 {url}")

    with patch.object(eg, "fetch_text", side_effect=fake_fetch):
        url, division = eg._find_species_genome(
            "https://ftp.ensembl.org/pub", "current", "homo_sapiens"
        )

    assert "dna_sm" in url
    assert division == "Ensembl"


# ---------------------------------------------------------------------------
# _find_species_genome — EnsemblGenomes fallback
# ---------------------------------------------------------------------------

def test_find_species_genome_ensemblgenomes_fallback():
    genome_file = "Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa.gz"

    def fake_fetch(url, timeout=60):
        if "ensemblgenomes" in url and "/caenorhabditis_elegans/dna/" in url:
            return _dir_html(genome_file)
        raise eg.EnsemblSpeciesNotFound(f"404 {url}")

    with patch.object(eg, "fetch_text", side_effect=fake_fetch):
        url, division = eg._find_species_genome(
            "https://ftp.ensembl.org/pub", "current", "caenorhabditis_elegans"
        )

    assert "dna_sm" in url
    assert "EnsemblGenomes" in division


# ---------------------------------------------------------------------------
# _find_species_genome — prefix fallback
# ---------------------------------------------------------------------------

def test_find_species_genome_prefix_fallback():
    genome_file = "Mus_musculus_c57bl6j.C57BL_6J_v1.dna_sm.toplevel.fa.gz"

    def fake_fetch(url, timeout=60):
        # Exact species dir doesn't exist
        if "/mus_musculus/dna/" in url:
            raise eg.EnsemblSpeciesNotFound("404")
        # Root listing has strain variant
        if url.endswith("current_fasta/") and "ensemblgenomes" not in url:
            return _species_listing("mus_musculus_c57bl6j", "homo_sapiens")
        # Strain variant dna/ dir exists
        if "/mus_musculus_c57bl6j/dna/" in url:
            return _dir_html(genome_file)
        raise eg.EnsemblSpeciesNotFound(f"404 {url}")

    with patch.object(eg, "fetch_text", side_effect=fake_fetch):
        url, division = eg._find_species_genome(
            "https://ftp.ensembl.org/pub", "current", "mus_musculus"
        )

    assert "dna_sm" in url


# ---------------------------------------------------------------------------
# _find_species_genome — not found
# ---------------------------------------------------------------------------

def test_find_species_genome_not_found():
    with patch.object(eg, "fetch_text", side_effect=eg.EnsemblSpeciesNotFound("404")):
        with pytest.raises(eg.EnsemblSpeciesNotFound):
            eg._find_species_genome(
                "https://ftp.ensembl.org/pub", "current", "unknown_species"
            )
