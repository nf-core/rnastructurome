import importlib.util
from pathlib import Path


_SCRIPT = Path(__file__).resolve().parents[5] / "bin" / "ensembl_gtf.py"
spec = importlib.util.spec_from_file_location("ensembl_gtf", _SCRIPT)
eg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(eg)


def test_release_path_current_uses_main_ensembl_current_gtf_layout():
    assert eg.release_path_for_value("current") == "current/gtf"


def test_release_path_latest_uses_main_ensembl_current_gtf_layout():
    assert eg.release_path_for_value("latest") == "current/gtf"


def test_release_path_numbered_release_uses_release_gtf_layout():
    assert eg.release_path_for_value("115") == "release-115/gtf"


def test_release_path_prefixed_release_uses_release_gtf_layout():
    assert eg.release_path_for_value("release-115") == "release-115/gtf"
