import shutil
import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[6] / "bin" / "rctools_covered_bed.py"
FIXTURES = Path(__file__).resolve().parent / "fixtures"


def _run(tmp_path, covered_ids, *, attr_name="transcript_id"):
    gtf = tmp_path / "test.gtf"
    shutil.copy(FIXTURES / "test.gtf", gtf)
    (tmp_path / "covered_ids.txt").write_text("\n".join(covered_ids) + "\n")

    result = subprocess.run(
        [sys.executable, str(SCRIPT), str(gtf), "exon", attr_name],
        cwd=tmp_path,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    return (tmp_path / "covered.bed").read_text()


def test_covered_bed_uses_transcript_id_and_spliced_exon_length(tmp_path):
    out = _run(tmp_path, ["YP_009724389.1"])

    assert out == "YP_009724389.1\t0\t21290\tYP_009724389.1\n"


def test_covered_bed_can_use_gene_id_attribute(tmp_path):
    out = _run(tmp_path, ["S"], attr_name="gene_id")

    assert out == "S\t0\t3822\tS\n"


def test_covered_bed_omits_uncovered_transcripts(tmp_path):
    out = _run(tmp_path, ["missing_transcript"])

    assert out == ""
