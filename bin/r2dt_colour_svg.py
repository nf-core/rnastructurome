#!/usr/bin/env python3
# Author: Victoria Begley (begley@ebi.ac.uk)
# License: MIT
"""
Overlay SHAPE / chemical-probing reactivities onto R2DT SVG diagrams.

Reads RNAframework rf-norm XML files and R2DT SVG outputs, then colours
each nucleotide according to its normalised reactivity value using the
standard SHAPE colouring scheme:

  NaN / missing / negative  →  grey     (#B1B3B6)
  0.00 – 0.40               →  black    (#000000)
  0.40 – 0.70               →  yellow   (#FFCD2F)
  > 0.70                    →  dark red (#9A2322)

When multiple XML files contain data for the same transcript (e.g. biological
replicates), per-position reactivities are averaged, ignoring NaN values.

Usage
-----
    r2dt_colour_svg.py \\
        --svg-dir       /path/to/r2dt/results/svg \\
        --xml-search-dir /path/to/work/dir \\
        --out-dir       /path/to/coloured_svg
"""

import argparse
import math
import re
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

# SVG namespace used by R2DT / Traveler
_SVG_NS = 'http://www.w3.org/2000/svg'
_TITLE_RE = re.compile(r'^(\d+)\s')   # leading integer = 1-based position
_INSERTION_RE = re.compile(r'([\d,]+)\s+nucleotides?\s+not\s+shown', re.IGNORECASE)
_INSERTED_RE  = re.compile(r'^\d+\s+\(inserted\)')   # individually inserted nucleotide
_NUCLEOTIDE_RE = re.compile(r'^[ACGUT]$')
_MAX_INSERTION_FRACTION = 0.35        # skip SVGs where >35% of nts are insertions
_MAX_HIDDEN_INSERTION_NT = 100         # skip SVGs with a single large hidden arc
_MIN_DRAWN_NUCLEOTIDES = 20            # skip near-empty template matches
_LEGEND_SCALE = 0.75
_LEGEND_WIDTH = 80
_LEGEND_HEIGHT = 54
_LEGEND_X_MARGIN = 5
_LEGEND_Y_MARGIN = 10
_LEGEND_ITEMS = (
    ('#9A2322', '0.7+'),
    ('#FFCD2F', '0.4-0.7'),
    ('#000000', '0-0.4'),
    ('#B1B3B6', 'No data'),
)


@dataclass(frozen=True)
class _SvgNucleotideStats:
    drawn: int = 0
    inserted: int = 0
    hidden: int = 0
    max_hidden: int = 0

    @property
    def total(self) -> int:
        return self.drawn + self.inserted + self.hidden

    @property
    def insertion_fraction(self) -> float:
        return (self.inserted + self.hidden) / self.total if self.total else 0.0


# ── SHAPE colour scheme ──────────────────────────────────────────────────────

def _shape_colour(r) -> str:
    if r is None or (isinstance(r, float) and math.isnan(r)) or r < 0:
        return '#B1B3B6'   # grey     – no data
    if r <= 0.40:
        return '#000000'   # black    – low reactivity
    if r < 0.70:
        return '#FFCD2F'   # yellow   – medium reactivity
    return '#9A2322'       # dark red – high reactivity


# ── rf-norm XML parser ───────────────────────────────────────────────────────

def _parse_rfnorm_xml(path: Path) -> dict:
    """Return {transcript_id: [float|NaN, ...]} from one rf-norm XML file."""
    try:
        tree = ET.parse(path)
    except ET.ParseError as exc:
        print(f'[WARN] Could not parse {path}: {exc}', file=sys.stderr)
        return {}

    result = {}
    for transcript in tree.getroot().iter('transcript'):
        tid = transcript.get('id') or transcript.findtext('id')
        if not tid:
            continue
        raw = transcript.findtext('reactivity') or transcript.findtext('values') or ''
        values = []
        for token in raw.strip().split(','):
            token = token.strip()
            if token.lower() in ('nan', 'na', '', 'null', 'none'):
                values.append(float('nan'))
            else:
                try:
                    values.append(float(token))
                except ValueError:
                    values.append(float('nan'))
        if values:
            result[tid] = values
    return result


def _index_reactivity_xml(xml_search_dir: Path) -> dict:
    """Return {transcript_id: [xml_path, ...]} using staged XML file stems."""
    xml_index: dict = defaultdict(list)
    for xml_path in xml_search_dir.rglob('xml_input*/*.xml'):
        xml_index[xml_path.stem].append(xml_path)
    return xml_index


def _r2dt_svg_tid(svg_path: Path, known_tids: Optional[set] = None) -> str:
    """Return transcript ID for an R2DT SVG, preserving IDs that contain hyphens."""
    name = svg_path.name
    for suffix in ('.colored.svg', '.enriched.svg', '.svg'):
        if name.endswith(suffix):
            base = name[:-len(suffix)]
            break
    else:
        base = svg_path.stem

    if known_tids:
        matches = [
            tid for tid in known_tids
            if base == tid or base.startswith(f'{tid}-')
        ]
        if matches:
            return max(matches, key=len)

    return base.split('-')[0]


def load_reactivities(xml_search_dir: Path, tids_needed: set, xml_index: Optional[dict] = None) -> dict:
    """
    Load rf-norm XML reactivities only for transcripts in tids_needed.
    Builds a filename index first (fast), then parses only the relevant files.
    Multiple XML files for the same transcript are averaged position-by-position.
    """
    if xml_index is None:
        xml_index = _index_reactivity_xml(xml_search_dir)

    if not xml_index:
        print('[WARN] No reactivity data found in xml_input* directories.', file=sys.stderr)
        return {}

    raw: dict = defaultdict(list)   # tid -> list of reactivity vectors
    for tid in tids_needed:
        for xml_path in xml_index.get(tid, []):
            for _, values in _parse_rfnorm_xml(xml_path).items():
                raw[tid].append(values)

    averaged = {}
    for tid, reps in raw.items():
        if len(reps) == 1:
            averaged[tid] = reps[0]
        else:
            length = max(len(r) for r in reps)
            avg = []
            for i in range(length):
                vals = [
                    r[i] for r in reps
                    if i < len(r) and r[i] is not None and not math.isnan(r[i])
                ]
                avg.append(sum(vals) / len(vals) if vals else float('nan'))
            averaged[tid] = avg

    return averaged


# ── SVG colouring ────────────────────────────────────────────────────────────

def _numeric_svg_dimension(value: Optional[str], default: float) -> float:
    """Return the leading numeric part of an SVG dimension."""
    if not value:
        return default
    m = re.match(r'\s*([0-9.]+)', value)
    return float(m.group(1)) if m else default


def _remove_existing_legend(root) -> None:
    """Remove any existing legend group so repeated colouring is idempotent."""
    parent_map = {child: parent for parent in root.iter() for child in parent}
    for g in list(root.iter(f'{{{_SVG_NS}}}g')):
        classes = set((g.get('class') or '').split())
        if 'legend' in classes:
            parent = parent_map.get(g)
            if parent is not None:
                parent.remove(g)


def _append_legend(root) -> None:
    """Append a SHAPE reactivity colour legend matching ViennaRNA outputs."""
    _remove_existing_legend(root)

    svg_w = _numeric_svg_dimension(root.get('width'), 452.0)
    svg_h = _numeric_svg_dimension(root.get('height'), 452.0)
    tx = svg_w - _LEGEND_WIDTH * _LEGEND_SCALE - _LEGEND_X_MARGIN
    ty = svg_h - _LEGEND_HEIGHT * _LEGEND_SCALE - _LEGEND_Y_MARGIN

    legend = ET.SubElement(root, f'{{{_SVG_NS}}}g')
    legend.set('transform', f'translate({tx:.1f},{ty:.1f}) scale({_LEGEND_SCALE})')
    legend.set('class', 'legend')

    for idx, (colour, label) in enumerate(_LEGEND_ITEMS):
        yo = idx * 14
        rect = ET.SubElement(legend, f'{{{_SVG_NS}}}rect')
        rect.set('x', '0')
        rect.set('y', str(yo))
        rect.set('width', '12')
        rect.set('height', '12')
        rect.set('fill', colour)
        rect.set('stroke', 'black')

        text = ET.SubElement(legend, f'{{{_SVG_NS}}}text')
        text.set('x', '18')
        text.set('y', str(yo + 10))
        text.set('class', 'nucleotide')
        text.set('style', 'font-family: Arial, Helvetica, sans-serif; font-size: 11px; fill: #000000;')
        text.text = label


def _nucleotide_text(g):
    """Return the nucleotide <text> element for an R2DT position group, if any."""
    text_el = g.find(f'{{{_SVG_NS}}}text')
    if text_el is None or text_el.text is None:
        return None
    return text_el if _NUCLEOTIDE_RE.match(text_el.text.strip()) else None


def _svg_nucleotide_stats(root) -> _SvgNucleotideStats:
    """Return nucleotide counts used to decide whether an R2DT SVG is publishable.

    Counts both bulk insertion arcs (<g class="insertion-arc"> with
    title "N nucleotides not shown") and individually inserted nucleotides
    (title "N (inserted)"). Template-matched positions contribute to the
    denominator but not the numerator. 5'/3' labels are ignored.
    """
    inserted = 0
    hidden = 0
    max_hidden = 0
    drawn = 0
    for g in root.iter(f'{{{_SVG_NS}}}g'):
        title_el = g.find(f'{{{_SVG_NS}}}title')
        if title_el is None or not title_el.text:
            continue
        title = title_el.text.strip()
        if g.get('class') == 'insertion-arc':
            m = _INSERTION_RE.search(title)
            if m:
                count = int(m.group(1).replace(',', ''))
                hidden += count
                max_hidden = max(max_hidden, count)
            continue

        if _nucleotide_text(g) is None:
            continue
        if _INSERTED_RE.match(title):
            inserted += 1
        elif _TITLE_RE.match(title):
            drawn += 1

    return _SvgNucleotideStats(
        drawn=drawn,
        inserted=inserted,
        hidden=hidden,
        max_hidden=max_hidden,
    )


def _svg_skip_reason(root) -> Optional[str]:
    """Return a human-readable reason to skip a low-quality R2DT SVG, if any."""
    stats = _svg_nucleotide_stats(root)
    if stats.drawn < _MIN_DRAWN_NUCLEOTIDES:
        return f'only {stats.drawn} template-position nucleotides drawn'
    if stats.max_hidden > _MAX_HIDDEN_INSERTION_NT:
        return f'largest hidden insertion arc is {stats.max_hidden:,} nt'
    if stats.insertion_fraction > _MAX_INSERTION_FRACTION:
        return f'{stats.insertion_fraction:.0%} of nucleotides are insertions'
    return None


def _insertion_fraction(root) -> float:
    """Return fraction of total nucleotides that are insertions (0–1)."""
    return _svg_nucleotide_stats(root).insertion_fraction


def colour_svg(svg_path: Path, reactivities: list, out_path: Path) -> int:
    """
    Colour nucleotide text elements in an R2DT SVG by SHAPE reactivity.
    Nucleotides are identified by the leading integer in their <title> text,
    which gives the 1-based sequence position.
    Returns the number of non-grey nucleotides coloured, or 0 if the SVG
    was skipped (insertion arc threshold, all-NaN at SVG positions, or parse error).
    """
    ET.register_namespace('', _SVG_NS)
    ET.register_namespace('xlink', 'http://www.w3.org/1999/xlink')

    try:
        tree = ET.parse(svg_path)
    except ET.ParseError as exc:
        print(f'[WARN] Could not parse SVG {svg_path}: {exc}', file=sys.stderr)
        return 0

    skip_reason = _svg_skip_reason(tree.getroot())
    if skip_reason:
        print(
            f'[R2DT colour] Skipping {svg_path.name}: {skip_reason} — '
            f'will be drawn by ViennaRNA instead',
            file=sys.stderr,
        )
        return 0

    # First pass: resolve colours for every SVG nucleotide position.
    # The reactivity array may be longer than the SVG (e.g. a tRNA embedded in a
    # larger transcript); only the positions actually drawn in the SVG matter.
    pending: list[tuple] = []   # (text_el, colour)
    for g in tree.getroot().iter(f'{{{_SVG_NS}}}g'):
        title_el = g.find(f'{{{_SVG_NS}}}title')
        if title_el is None or not title_el.text:
            continue
        m = _TITLE_RE.match(title_el.text.strip())
        if not m:
            continue
        pos = int(m.group(1)) - 1          # 0-based index
        if pos < 0 or pos >= len(reactivities):
            continue
        text_el = g.find(f'{{{_SVG_NS}}}text')
        if text_el is None or _nucleotide_text(g) is None:
            continue
        pending.append((text_el, _shape_colour(reactivities[pos])))

    if not pending:
        print(
            f'[R2DT colour] Skipping {svg_path.name}: no drawable nucleotide '
            f'positions with reactivity data — will be drawn by ViennaRNA instead',
            file=sys.stderr,
        )
        return 0

    # Skip if every drawn position would be grey — the reactivity array may have
    # real values elsewhere (e.g. flanking exons) but not at the template positions.
    if pending and not any(colour != '#B1B3B6' for _, colour in pending):
        print(
            f'[R2DT colour] Skipping {svg_path.name}: all-NaN reactivity at '
            f'SVG positions — will be drawn by ViennaRNA instead',
            file=sys.stderr,
        )
        return 0

    # Second pass: apply colours.
    for text_el, colour in pending:
        existing = re.sub(r'fill\s*:[^;]+;?\s*', '', text_el.get('style', '')).strip()
        new_style = f'fill: {colour};'
        if existing:
            new_style = f'{new_style} {existing}'
        text_el.set('style', new_style)

    _append_legend(tree.getroot())

    root = tree.getroot()
    bg = ET.Element(f'{{{_SVG_NS}}}rect')
    bg.set('width', '100%')
    bg.set('height', '100%')
    bg.set('fill', 'white')
    root.insert(0, bg)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    tree.write(str(out_path), xml_declaration=True, encoding='unicode')
    return sum(1 for _, colour in pending if colour != '#B1B3B6')


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    """Overlay SHAPE reactivities onto R2DT SVGs and report colouring statistics."""
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--svg-dir',        required=True, type=Path,
                    help='Directory containing R2DT SVG output files')
    ap.add_argument('--xml-search-dir', required=True, type=Path,
                    help='Root directory to search recursively for rf-norm XML files '
                         '(looks for xml_input*/*.xml pattern)')
    ap.add_argument('--out-dir',        required=True, type=Path,
                    help='Output directory for reactivity-coloured SVGs')
    args = ap.parse_args()

    svg_paths = sorted(args.svg_dir.glob('*.colored.svg'))
    if not svg_paths:
        print('[R2DT colour] 0 SVGs coloured, 0 skipped (no SVGs found)', file=sys.stderr)
        return

    xml_index = _index_reactivity_xml(args.xml_search_dir)
    known_tids = set(xml_index.keys())
    svg_tids = {
        svg_path: _r2dt_svg_tid(svg_path, known_tids)
        for svg_path in svg_paths
    }
    tids_needed = set(svg_tids.values())

    reactivities = load_reactivities(args.xml_search_dir, tids_needed, xml_index=xml_index)
    if not reactivities:
        print('[ERROR] No reactivity data loaded — cannot colour SVGs.', file=sys.stderr)
        sys.exit(1)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    n_coloured = 0
    n_skipped  = 0

    for svg_path in svg_paths:
        # R2DT names SVGs as {transcript_id}-{template}.colored.svg.
        tid = svg_tids[svg_path]
        if tid not in reactivities:
            n_skipped += 1
            continue
        vals = reactivities[tid]
        if not any(not math.isnan(v) and v >= 0 for v in vals):
            print(
                f'[R2DT colour] Skipping {tid}: all-NaN reactivity — '
                f'no experimental support for template model, will fall back to ViennaRNA',
                file=sys.stderr,
            )
            n_skipped += 1
            continue
        n = colour_svg(svg_path, vals, args.out_dir / (tid + '.svg'))
        if n:
            n_coloured += 1
        else:
            n_skipped += 1

    print(
        f'[R2DT colour] {n_coloured} SVGs coloured, {n_skipped} skipped '
        f'(insertion arc threshold exceeded, all-NaN reactivity, or empty SVG)',
        file=sys.stderr,
    )


if __name__ == '__main__':
    main()
