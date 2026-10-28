"""Guard tests for bin/rnaframework_rfeval_window.py (rf-eval reactivity windowing)."""
import re
import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[5] / "bin" / "rnaframework_rfeval_window.py"

SEQ = "ACGTACGTAC" "GGGGTTTTAA" "CCCCAAAATT"  # 30 nt
REACT = [f"{i / 100:.3f}" for i in range(30)]


def _xml(seq, react, tid="chr1"):
    return (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<data combined="FALSE" max="10" norm="raw" tool="rf-norm">\n'
        f'\t<transcript id="{tid}" length="{len(seq)}">\n'
        f"\t\t<sequence>{seq}</sequence>\n"
        f'\t\t<reactivity>{",".join(react)}</reactivity>\n'
        "\t</transcript>\n</data>\n"
    )


def _run(tmp_path, windows_text, xmls=None):
    xmls = xmls if xmls is not None else {"chr1": _xml(SEQ, REACT)}
    staged = tmp_path / "xml_input1"
    staged.mkdir()
    for name, text in xmls.items():
        (staged / f"{name}.xml").write_text(text)
    manifest = tmp_path / "windows.tsv"
    manifest.write_text(windows_text)
    return subprocess.run(
        [
            sys.executable, str(SCRIPT),
            "--windows", str(manifest),
            "--xml-glob", str(tmp_path / "xml_input*/*.xml"),
            "--outdir", str(tmp_path / "out"),
        ],
        check=False,
        capture_output=True,
        text=True,
    )


def _load(path):
    text = Path(path).read_text()
    seq = "".join(re.search(r"<sequence>(.*?)</sequence>", text, re.S).group(1).split())
    body = re.search(r"<reactivity>(.*?)</reactivity>", text, re.S).group(1)
    react = [v for v in re.split(r"[,\s]+", body.strip()) if v]
    tid, length = re.search(r'id="([^"]+)" length="(\d+)"', text).groups()
    return tid, int(length), seq, react


def test_plus_strand_window_is_a_one_based_inclusive_slice(tmp_path):
    result = _run(tmp_path, "chr1\t5\t12\telement_a\t+\n")
    assert result.returncode == 0, result.stderr
    tid, length, seq, react = _load(tmp_path / "out" / "element_a.xml")
    assert tid == "element_a"  # id is rewritten to match the .db entry
    assert (seq, react) == (SEQ[4:12], REACT[4:12])
    assert length == len(seq) == 8


def test_minus_strand_window_is_reverse_complemented(tmp_path):
    result = _run(tmp_path, "chr1\t5\t12\telement_b\t-\n")
    assert result.returncode == 0, result.stderr
    _, _, seq, react = _load(tmp_path / "out" / "element_b.xml")
    complement = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")
    assert seq == SEQ[4:12].translate(complement)[::-1]
    assert react == REACT[4:12][::-1]


def test_strand_column_is_optional_and_defaults_to_plus(tmp_path):
    result = _run(tmp_path, "chr1 1 4 element_c\n")
    assert result.returncode == 0, result.stderr
    _, _, seq, _ = _load(tmp_path / "out" / "element_c.xml")
    assert seq == SEQ[0:4]


def test_header_row_and_comments_are_ignored(tmp_path):
    manifest = "# a comment\nref_seq_id\tstart\tend\tstructure_id\nchr1\t1\t4\telement_d\n"
    result = _run(tmp_path, manifest)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "out" / "element_d.xml").exists()


def test_nan_reactivities_are_carried_through(tmp_path):
    react = ["NaN", "NaN"] + REACT[2:]
    result = _run(tmp_path, "chr1\t1\t4\telement_e\t+\n", {"chr1": _xml(SEQ, react)})
    assert result.returncode == 0, result.stderr
    _, _, _, out_react = _load(tmp_path / "out" / "element_e.xml")
    assert out_react == ["NaN", "NaN", "0.020", "0.030"]


def test_out_of_range_window_is_skipped_not_truncated(tmp_path):
    result = _run(tmp_path, "chr1\t20\t99\toob\t+\nchr1\t1\t4\tgood\t+\n")
    assert result.returncode == 0, result.stderr
    assert not (tmp_path / "out" / "oob.xml").exists()
    assert (tmp_path / "out" / "good.xml").exists()
    assert "skip oob" in result.stderr


def test_reference_with_no_matching_xml_is_reported(tmp_path):
    result = _run(tmp_path, "missing\t1\t4\tnowhere\t+\nchr1\t1\t4\tgood\t+\n")
    assert result.returncode == 0, result.stderr
    assert "no XML found" in result.stderr
    assert (tmp_path / "out" / "good.xml").exists()


def test_reactivity_length_mismatch_drops_the_reference(tmp_path):
    bad = _xml(SEQ, REACT[:20])  # 30 nt sequence, only 20 reactivities
    result = _run(tmp_path, "chr1\t1\t4\ttruncated\t+\n", {"chr1": bad})
    assert result.returncode == 1
    assert "20 reactivities vs 30 nt" in result.stderr


def test_exits_nonzero_when_no_window_is_written(tmp_path):
    result = _run(tmp_path, "chr1\t99\t120\tnope\t+\n")
    assert result.returncode == 1
    assert "no windows written" in result.stderr


def test_short_manifest_row_is_a_hard_error(tmp_path):
    result = _run(tmp_path, "chr1\t1\t4\n")
    assert result.returncode == 1
    assert "need >=4 columns" in result.stderr
