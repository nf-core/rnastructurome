import math
import subprocess
import sys
import textwrap
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[4]
SCRIPT = REPO_ROOT / "bin" / "r2dt_colour_svg.py"

_SVG_NS = "http://www.w3.org/2000/svg"

sys.path.insert(0, str(REPO_ROOT / "bin"))
from r2dt_colour_svg import (  # noqa: E402
    _parse_rfnorm_xml,
    _r2dt_svg_tid,
    _shape_colour,
    _svg_nucleotide_stats,
    colour_svg,
    load_reactivities,
)


def _make_svg(
    tmp_path: Path,
    name: str,
    positions: list[int],
    *,
    hidden_arcs: list[int] | None = None,
    inserted_positions: list[int] | None = None,
) -> Path:
    """Create a minimal R2DT-style SVG with nucleotide groups."""
    groups = [
        (
            "  <g>\n"
            f"    <title>{pos} (position.label in template: {pos}.A)</title>\n"
            '    <text x="0" y="0">A</text>\n'
            "  </g>"
        )
        for pos in positions
    ]
    for pos in inserted_positions or []:
        groups.append(
            "  <g>\n"
            f"    <title>{pos} (inserted)</title>\n"
            '    <text x="0" y="0">A</text>\n'
            "  </g>"
        )
    for count in hidden_arcs or []:
        groups.append(
            '  <g class="insertion-arc">\n'
            f"    <title>{count:,} nucleotides not shown</title>\n"
            '    <path d="M 0,0 L 1,1" />\n'
            "  </g>"
        )

    groups_text = "\n".join(groups)
    svg = (
        '<?xml version="1.0" encoding="utf-8"?>\n'
        f'<svg xmlns="{_SVG_NS}" xmlns:xlink="http://www.w3.org/1999/xlink">\n'
        f"{groups_text}\n"
        "</svg>\n"
    )
    p = tmp_path / name
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(svg)
    return p


def _make_xml(tmp_path: Path, tid: str, values: list) -> Path:
    """Create a minimal rf-norm XML file."""
    reactivity_str = ",".join(
        "nan" if (v is None or (isinstance(v, float) and math.isnan(v))) else str(v)
        for v in values
    )
    xml = textwrap.dedent(f"""\
        <?xml version="1.0" encoding="utf-8"?>
        <data>
            <transcript id="{tid}">
                <length>{len(values)}</length>
                <reactivity>{reactivity_str}</reactivity>
            </transcript>
        </data>
    """)
    p = tmp_path / f"{tid}.xml"
    p.write_text(xml)
    return p


def _legend_labels(svg_path: Path) -> list[str]:
    root = ET.parse(svg_path).getroot()
    labels = []
    for g in root.iter(f"{{{_SVG_NS}}}g"):
        if "legend" not in (g.get("class") or "").split():
            continue
        labels.extend(
            text.text
            for text in g.iter(f"{{{_SVG_NS}}}text")
            if text.text is not None
        )
    return labels


def test_shape_colour_nan():
    assert _shape_colour(float("nan")) == "#B1B3B6"


def test_shape_colour_none():
    assert _shape_colour(None) == "#B1B3B6"


def test_shape_colour_negative():
    assert _shape_colour(-0.1) == "#B1B3B6"


def test_shape_colour_low():
    assert _shape_colour(0.0) == "#000000"
    assert _shape_colour(0.40) == "#000000"


def test_shape_colour_medium():
    assert _shape_colour(0.41) == "#FFCD2F"
    assert _shape_colour(0.69) == "#FFCD2F"


def test_shape_colour_high():
    assert _shape_colour(0.70) == "#9A2322"
    assert _shape_colour(1.0) == "#9A2322"


def test_parse_rfnorm_xml_basic(tmp_path):
    xml_path = _make_xml(tmp_path, "URS001", [0.1, float("nan"), 0.9])
    result = _parse_rfnorm_xml(xml_path)
    assert "URS001" in result
    vals = result["URS001"]
    assert vals[0] == pytest.approx(0.1)
    assert math.isnan(vals[1])
    assert vals[2] == pytest.approx(0.9)


def test_parse_rfnorm_xml_invalid_file(tmp_path):
    bad = tmp_path / "bad.xml"
    bad.write_text("not valid xml <<<")
    result = _parse_rfnorm_xml(bad)
    assert result == {}


def test_parse_rfnorm_xml_multiple_transcripts(tmp_path):
    xml_content = textwrap.dedent("""\
        <?xml version="1.0" encoding="utf-8"?>
        <data>
            <transcript id="A"><length>2</length><reactivity>0.1,0.2</reactivity></transcript>
            <transcript id="B"><length>2</length><reactivity>nan,0.5</reactivity></transcript>
        </data>
    """)
    p = tmp_path / "multi.xml"
    p.write_text(xml_content)
    result = _parse_rfnorm_xml(p)
    assert set(result.keys()) == {"A", "B"}
    assert result["A"] == pytest.approx([0.1, 0.2])
    assert math.isnan(result["B"][0])
    assert result["B"][1] == pytest.approx(0.5)


def test_load_reactivities_single_replicate(tmp_path):
    xml_dir = tmp_path / "xml_input_rep1"
    xml_dir.mkdir()
    _make_xml(xml_dir, "URS001", [0.1, 0.5, 0.9])
    result = load_reactivities(tmp_path, {"URS001"})
    assert "URS001" in result
    assert result["URS001"] == pytest.approx([0.1, 0.5, 0.9])


def test_load_reactivities_averages_replicates(tmp_path):
    for i, vals in enumerate(([0.0, 0.6], [0.2, 0.8])):
        d = tmp_path / f"xml_input_rep{i}"
        d.mkdir()
        _make_xml(d, "URS001", vals)
    result = load_reactivities(tmp_path, {"URS001"})
    assert result["URS001"] == pytest.approx([0.1, 0.7])


def test_load_reactivities_nan_ignored_in_average(tmp_path):
    d1 = tmp_path / "xml_input_rep1"
    d1.mkdir()
    _make_xml(d1, "URS001", [float("nan"), 0.2])
    d2 = tmp_path / "xml_input_rep2"
    d2.mkdir()
    _make_xml(d2, "URS001", [0.4, 0.6])
    result = load_reactivities(tmp_path, {"URS001"})
    assert result["URS001"][0] == pytest.approx(0.4)
    assert result["URS001"][1] == pytest.approx(0.4)


def test_load_reactivities_empty_dir(tmp_path):
    result = load_reactivities(tmp_path, {"URS001"})
    assert result == {}


def test_load_reactivities_ignores_non_xml_input_dirs(tmp_path):
    d = tmp_path / "some_other_dir"
    d.mkdir()
    _make_xml(d, "URS001", [0.5])
    result = load_reactivities(tmp_path, {"URS001"})
    assert result == {}


def test_r2dt_svg_tid_preserves_hyphenated_transcript_id(tmp_path):
    svg_path = tmp_path / "URS-001-template-name.colored.svg"

    tid = _r2dt_svg_tid(svg_path, {"URS-001"})

    assert tid == "URS-001"


def test_svg_nucleotide_stats_counts_insertions_and_ignores_labels(tmp_path):
    svg_path = _make_svg(
        tmp_path / "in",
        "URS001.svg",
        list(range(1, 21)),
        hidden_arcs=[1002],
        inserted_positions=[21, 22],
    )
    root = ET.parse(svg_path).getroot()

    stats = _svg_nucleotide_stats(root)

    assert stats.drawn == 20
    assert stats.inserted == 2
    assert stats.hidden == 1002
    assert stats.max_hidden == 1002
    assert stats.insertion_fraction == pytest.approx(1004 / 1024)


def test_colour_svg_colours_nucleotides(tmp_path):
    svg_path = _make_svg(tmp_path / "in", "URS001.svg", list(range(1, 21)))
    reactivities = [0.2, 0.6, 0.9] + [0.1] * 17
    out_path = tmp_path / "out" / "URS001.svg"

    n = colour_svg(svg_path, reactivities, out_path)

    assert n == 20
    tree = ET.parse(out_path)
    texts = [
        text for text in tree.getroot().iter(f"{{{_SVG_NS}}}text")
        if text.text == "A"
    ]
    assert len(texts) == 20
    styles = [t.get("style", "") for t in texts]
    assert any("#000000" in s for s in styles)
    assert any("#FFCD2F" in s for s in styles)
    assert any("#9A2322" in s for s in styles)
    assert _legend_labels(out_path) == ["0.7+", "0.4-0.7", "0-0.4", "No data"]


def test_colour_svg_all_nan_reactivity_skips(tmp_path):
    svg_path = _make_svg(tmp_path / "in", "URS001.svg", list(range(1, 21)))
    out_path = tmp_path / "out" / "URS001.svg"

    n = colour_svg(svg_path, [float("nan")] * 20, out_path)

    assert n == 0
    assert not out_path.exists()


def test_colour_svg_out_of_range_position_skipped(tmp_path):
    svg_path = _make_svg(tmp_path / "in", "URS001.svg", list(range(101, 121)))
    out_path = tmp_path / "out" / "URS001.svg"

    n = colour_svg(svg_path, [0.1, 0.2, 0.3], out_path)

    assert n == 0
    assert not out_path.exists()


def test_colour_svg_invalid_svg(tmp_path):
    bad = tmp_path / "bad.svg"
    bad.write_text("not xml <<<")
    n = colour_svg(bad, [0.5], tmp_path / "out.svg")
    assert n == 0


def test_colour_svg_creates_output_dir(tmp_path):
    svg_path = _make_svg(tmp_path / "in", "URS001.svg", list(range(1, 21)))
    nested_out = tmp_path / "deep" / "nested" / "URS001.svg"
    colour_svg(svg_path, [0.5] * 20, nested_out)
    assert nested_out.exists()


def test_colour_svg_skips_near_empty_template_match(tmp_path):
    svg_path = _make_svg(tmp_path / "in", "URS001.svg", [1])
    out_path = tmp_path / "out" / "URS001.svg"

    n = colour_svg(svg_path, [0.5], out_path)

    assert n == 0
    assert not out_path.exists()


def test_colour_svg_skips_large_hidden_insertion_arc(tmp_path):
    svg_path = _make_svg(
        tmp_path / "in",
        "URS001.svg",
        list(range(1, 201)),
        hidden_arcs=[101],
    )
    out_path = tmp_path / "out" / "URS001.svg"

    n = colour_svg(svg_path, [0.5] * 200, out_path)

    assert n == 0
    assert not out_path.exists()


def test_colour_svg_skips_high_insertion_fraction(tmp_path):
    svg_path = _make_svg(
        tmp_path / "in",
        "URS001.svg",
        list(range(1, 31)),
        inserted_positions=list(range(31, 51)),
    )
    out_path = tmp_path / "out" / "URS001.svg"

    n = colour_svg(svg_path, [0.5] * 50, out_path)

    assert n == 0
    assert not out_path.exists()


def test_colour_svg_allows_modest_insertions(tmp_path):
    svg_path = _make_svg(
        tmp_path / "in",
        "URS001.svg",
        list(range(1, 31)),
        inserted_positions=list(range(31, 36)),
    )
    out_path = tmp_path / "out" / "URS001.svg"

    n = colour_svg(svg_path, [0.5] * 35, out_path)

    assert n == 35
    assert out_path.exists()


def test_colour_svg_replaces_existing_legend(tmp_path):
    svg_path = _make_svg(tmp_path / "in", "URS001.svg", list(range(1, 21)))
    out_path = tmp_path / "out" / "URS001.svg"

    colour_svg(svg_path, [0.5] * 20, out_path)
    colour_svg(out_path, [0.5] * 20, out_path)

    root = ET.parse(out_path).getroot()
    legends = [
        g for g in root.iter(f"{{{_SVG_NS}}}g")
        if "legend" in (g.get("class") or "").split()
    ]
    assert len(legends) == 1
    assert _legend_labels(out_path) == ["0.7+", "0.4-0.7", "0-0.4", "No data"]


def test_cli_end_to_end(tmp_path):
    svg_dir = tmp_path / "svg"
    svg_dir.mkdir()
    _make_svg(svg_dir, "URS001-template.colored.svg", list(range(1, 21)))

    xml_input_dir = tmp_path / "xml_input_rep1"
    xml_input_dir.mkdir()
    _make_xml(xml_input_dir, "URS001", [0.1, 0.5, 0.9] + [0.1] * 17)

    out_dir = tmp_path / "coloured"

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--svg-dir",
            str(svg_dir),
            "--xml-search-dir",
            str(tmp_path),
            "--out-dir",
            str(out_dir),
        ],
        check=False,
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0, result.stderr
    assert (out_dir / "URS001.svg").exists()
    assert _legend_labels(out_dir / "URS001.svg") == ["0.7+", "0.4-0.7", "0-0.4", "No data"]


def test_cli_exits_nonzero_when_no_reactivities(tmp_path):
    svg_dir = tmp_path / "svg"
    svg_dir.mkdir()
    _make_svg(svg_dir, "URS001-template.colored.svg", list(range(1, 21)))
    out_dir = tmp_path / "coloured"

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--svg-dir",
            str(svg_dir),
            "--xml-search-dir",
            str(tmp_path),
            "--out-dir",
            str(out_dir),
        ],
        check=False,
        text=True,
        capture_output=True,
    )

    assert result.returncode != 0


def test_cli_skips_svg_with_no_matching_transcript(tmp_path):
    svg_dir = tmp_path / "svg"
    svg_dir.mkdir()
    _make_svg(svg_dir, "URS001-template.colored.svg", list(range(1, 21)))
    _make_svg(svg_dir, "UNKNOWN-template.colored.svg", list(range(1, 21)))

    xml_input_dir = tmp_path / "xml_input_rep1"
    xml_input_dir.mkdir()
    _make_xml(xml_input_dir, "URS001", [0.5] * 20)

    out_dir = tmp_path / "coloured"

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--svg-dir",
            str(svg_dir),
            "--xml-search-dir",
            str(tmp_path),
            "--out-dir",
            str(out_dir),
        ],
        check=False,
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0
    assert (out_dir / "URS001.svg").exists()
    assert not (out_dir / "UNKNOWN.svg").exists()


def test_cli_matches_hyphenated_transcript_id(tmp_path):
    svg_dir = tmp_path / "svg"
    svg_dir.mkdir()
    _make_svg(svg_dir, "URS-001-template.colored.svg", list(range(1, 21)))

    xml_input_dir = tmp_path / "xml_input_rep1"
    xml_input_dir.mkdir()
    _make_xml(xml_input_dir, "URS-001", [0.5] * 20)

    out_dir = tmp_path / "coloured"

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--svg-dir",
            str(svg_dir),
            "--xml-search-dir",
            str(tmp_path),
            "--out-dir",
            str(out_dir),
        ],
        check=False,
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0, result.stderr
    assert (out_dir / "URS-001.svg").exists()
