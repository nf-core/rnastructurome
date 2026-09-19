"""Guard tests for bin/rfnormfactor_align_rc.sh (RC alignment ahead of rf-normfactor)."""
import shutil
import struct
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[5] / "bin" / "rfnormfactor_align_rc.sh"

pytestmark = pytest.mark.skipif(shutil.which("rf-rctools") is None, reason="rf-rctools not on PATH")


def _write_rc(path, transcripts):
    """Minimal RC file: per-transcript header, 2-bit-per-nibble sequence, stops, coverage, reads."""
    out = bytearray()
    for tid, seq in transcripts.items():
        hexseq = seq.translate(str.maketrans("ACGTN", "01234"))
        out += struct.pack("<L", len(tid) + 1) + tid.encode() + b"\0" + struct.pack("<L", len(seq))
        out += bytes.fromhex(hexseq + "0" * (len(hexseq) % 2))
        out += struct.pack(f"<{2 * len(seq)}L", *([1] * len(seq)), *([10] * len(seq)))
        out += struct.pack("<L", 10)
    out += struct.pack("<Q", 10 * len(transcripts)) + struct.pack("<H", 1) + b"[eofrc]"
    path.write_bytes(out)


def _view_ids(rc):
    lines = subprocess.run(["rf-rctools", "view", str(rc)], check=True, capture_output=True, text=True).stdout.splitlines()
    return [lines[i] for i in range(len(lines)) if i % 5 == 0 and lines[i]]


def test_aligns_rc_files_whose_names_start_with_a_digit(tmp_path):
    _write_rc(tmp_path / "125ng_NAIN3_r1.rc", {"tx1": "ACGTACGT", "tx2": "GGGCCC"})
    _write_rc(tmp_path / "125ng_DMSO_r1.rc", {"tx1": "ACGTACGT", "tx3": "TTTAAA"})

    subprocess.run(
        [str(SCRIPT), "homo_sapiens", "aligned", "125ng_NAIN3_r1.rc", "125ng_DMSO_r1.rc"],
        cwd=tmp_path, check=True, capture_output=True, text=True,
    )

    assert _view_ids(tmp_path / "aligned" / "125ng_NAIN3_r1.rc") == ["tx1"]
    assert _view_ids(tmp_path / "aligned" / "125ng_DMSO_r1.rc") == ["tx1"]
