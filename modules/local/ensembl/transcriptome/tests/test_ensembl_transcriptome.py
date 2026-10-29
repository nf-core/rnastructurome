import importlib.util
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

# Import the standalone script as a module
_SCRIPT = Path(__file__).resolve().parents[5] / "bin" / "ensembl_transcriptome.py"
spec = importlib.util.spec_from_file_location("ensembl_transcriptome", _SCRIPT)
et = importlib.util.module_from_spec(spec)
spec.loader.exec_module(et)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _html(*hrefs: str) -> str:
    """Minimal FTP-style HTML listing with the given directory hrefs."""
    links = "".join(f'<a href="{h}/">{h}/</a>\n' for h in hrefs)
    return f"<html><body>{links}</body></html>"


def _cdna_html(filename: str) -> str:
    return f'<html><body><a href="{filename}">{filename}</a></body></html>'


# ---------------------------------------------------------------------------
# _prefix_candidates — flat layout
# ---------------------------------------------------------------------------

def test_prefix_candidates_flat_exact_match():
    """Species whose name exactly equals the prefix stem is returned."""
    def fake_fetch(url, timeout=60):
        return _html("homo_sapiens", "mus_musculus", "homo_sapiens_alt")

    with patch.object(et, "fetch_text", side_effect=fake_fetch):
        results = et._prefix_candidates("https://ftp.ensembl.org/pub", "current", "homo_sapiens_")

    species = [s for _, s in results]
    assert "homo_sapiens_alt" in species


def test_prefix_candidates_flat_prefix_match():
    """Strain-suffixed name (e.g. escherichia_coli_k_12) is matched by prefix."""
    def fake_fetch(url, timeout=60):
        return _html("homo_sapiens", "escherichia_coli_k_12", "escherichia_coli_o157")

    with patch.object(et, "fetch_text", side_effect=fake_fetch):
        results = et._prefix_candidates("https://ftp.example.org/pub", "current", "escherichia_coli_")

    species = [s for _, s in results]
    assert "escherichia_coli_k_12" in species
    assert "escherichia_coli_o157" in species
    assert "homo_sapiens" not in species


def test_prefix_candidates_flat_no_match():
    def fake_fetch(url, timeout=60):
        return _html("homo_sapiens", "mus_musculus")

    with patch.object(et, "fetch_text", side_effect=fake_fetch):
        results = et._prefix_candidates("https://ftp.example.org/pub", "current", "escherichia_coli_")

    assert results == []


# ---------------------------------------------------------------------------
# _prefix_candidates — collection layout (bacteria-style)
# ---------------------------------------------------------------------------

def test_prefix_candidates_collection_layout():
    """Strain in a bacteria_N_collection is found via the collection scan."""
    def fake_fetch(url, timeout=60):
        if url.endswith("/current_fasta/"):
            return _html("bacteria_1_collection", "bacteria_2_collection")
        if "bacteria_1_collection/" in url and url.endswith("/"):
            return _html("bacillus_subtilis_168", "escherichia_coli_k_12")
        if "bacteria_2_collection/" in url and url.endswith("/"):
            return _html("salmonella_typhimurium")
        raise et.EnsemblSpeciesNotFound(f"unexpected url: {url}")

    with patch.object(et, "fetch_text", side_effect=fake_fetch):
        results = et._prefix_candidates(
            "https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria", "current", "escherichia_coli_"
        )

    collections = {c for c, _ in results}
    species = [s for _, s in results]
    assert "escherichia_coli_k_12" in species
    assert "bacteria_1_collection" in collections
    assert "bacillus_subtilis_168" not in species


def test_prefix_candidates_collection_unreachable_silently_skipped():
    """A collection that throws a network error is skipped, not fatal."""
    def fake_fetch(url, timeout=60):
        if url.endswith("/current_fasta/"):
            return _html("bacteria_1_collection", "bacteria_2_collection")
        if "bacteria_1_collection/" in url:
            raise Exception("network error")
        if "bacteria_2_collection/" in url:
            return _html("escherichia_coli_k_12")
        raise et.EnsemblSpeciesNotFound(f"unexpected url: {url}")

    with patch.object(et, "fetch_text", side_effect=fake_fetch):
        results = et._prefix_candidates(
            "https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria", "current", "escherichia_coli_"
        )

    species = [s for _, s in results]
    assert "escherichia_coli_k_12" in species


# ---------------------------------------------------------------------------
# _find_species — end-to-end prefix fallback
# ---------------------------------------------------------------------------

def test_find_species_bacteria_prefix_fallback():
    """'escherichia_coli' resolves to a strain via the prefix fallback.

    _prefix_candidates is mocked to isolate _find_species logic from
    the collection-scan details (which are covered by the tests above).
    fetch_text is mocked to serve the cdna directory listing once the
    correct strain URL is constructed.
    """
    def fake_prefix_candidates(base_url, release, prefix):
        if "bacteria" in base_url:
            return [("bacteria_1_collection", "escherichia_coli_k_12")]
        return []

    def fake_fetch(url, timeout=60):
        if "escherichia_coli_k_12/cdna/" in url:
            return _cdna_html("Escherichia_coli_K_12.test.cdna.all.fa.gz")
        raise et.EnsemblSpeciesNotFound(f"not found: {url}")

    with patch.object(et, "_prefix_candidates", side_effect=fake_prefix_candidates):
        with patch.object(et, "fetch_text", side_effect=fake_fetch):
            cdna_dir, cdna_name, ncrna_dir, label = et._find_species(
                "https://ftp.ensembl.org/pub", "current", "escherichia_coli"
            )

    assert "escherichia_coli_k_12" in cdna_dir
    assert cdna_name == "Escherichia_coli_K_12.test.cdna.all.fa.gz"
    assert label == "EnsemblBacteria"


def test_find_species_not_found_raises():
    """Species absent from all sources raises EnsemblSpeciesNotFound."""
    def fake_fetch(url, timeout=60):
        raise et.EnsemblSpeciesNotFound(f"not found: {url}")

    with patch.object(et, "fetch_text", side_effect=fake_fetch):
        with pytest.raises(et.EnsemblSpeciesNotFound):
            et._find_species("https://ftp.ensembl.org/pub", "current", "made_up_species")
