#!/usr/bin/env python3
# Author: Victoria Begley (begley@ebi.ac.uk)
# License: MIT
"""Colour a ViennaRNA RNAplot SVG with SHAPE reactivity bubbles.

Matches the RNAframework (Interface/ViennaRNA.pm) visual style:
  - Circles at cx=x+4, cy=y-4, r=8, fill-opacity=0.75
  - Grey (#B1B3B6) for no data, black (#000000) for 0-0.4,
    yellow (#FFCD2F) for 0.4-0.7, dark red (#9A2322) for >=0.7
  - Legend appended at bottom-right

Usage: viennarna_colour_svg.py <shape_file> <svg_file>
The SVG is modified in-place.
"""
from __future__ import annotations

import math
import re
import sys
import xml.etree.ElementTree as ET

_SVG_NS = 'http://www.w3.org/2000/svg'


def bubble_colour(r: float | None) -> str:
    if r is None or math.isnan(r) or r < 0:
        return '#B1B3B6'
    if r <= 0.40:
        return '#000000'
    if r < 0.70:
        return '#FFCD2F'
    return '#9A2322'


def text_colour(bg: str) -> str:
    return '#ffffff' if bg in ('#000000', '#9A2322') else '#000000'


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <shape_file> <svg_file>", file=sys.stderr)
        sys.exit(1)

    shape_file, svg_file = sys.argv[1], sys.argv[2]

    reactivities: dict[int, float] = {}
    with open(shape_file) as fh:
        for line in fh:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                try:
                    reactivities[int(parts[0]) - 1] = float(parts[1])
                except ValueError:
                    pass

    ET.register_namespace('', _SVG_NS)
    ET.register_namespace('xlink', 'http://www.w3.org/1999/xlink')
    tree = ET.parse(svg_file)
    root = tree.getroot()

    for style_el in root.iter(f'{{{_SVG_NS}}}style'):
        if style_el.text:
            css = style_el.text
            css = css.replace('SansSerif', 'Arial, Helvetica, sans-serif')
            css = re.sub(r'(\.basepairs\b[^}]*?)stroke\s*:\s*red', r'\1stroke: #808080', css)
            style_el.text = css

    parent_map = {child: parent for parent in root.iter() for child in parent}
    texts = [e for e in root.iter(f'{{{_SVG_NS}}}text')
             if e.get('class', '').strip() == 'nucleotide']

    for i, text_el in enumerate(texts):
        val = reactivities.get(i)
        fill = bubble_colour(val if val is not None else float('nan'))

        x = float(text_el.get('x', 0))
        y = float(text_el.get('y', 0))

        circle = ET.Element(f'{{{_SVG_NS}}}circle')
        circle.set('cx', f'{x + 4:.3f}')
        circle.set('cy', f'{y - 4:.3f}')
        circle.set('r', '8')
        circle.set('fill', fill)
        circle.set('fill-opacity', '0.75')
        circle.set('stroke', fill)
        circle.set('stroke-width', '0.8')

        parent = parent_map.get(text_el)
        if parent is not None:
            parent.insert(list(parent).index(text_el), circle)

        text_el.set('style', f'fill: {text_colour(fill)};')

    svg_w = float(root.get('width', '452'))
    svg_h = float(root.get('height', '452'))
    sc = 0.75
    tx = svg_w - 80 * sc - 5
    ty = svg_h - 54 * sc - 10
    legend = ET.SubElement(root, f'{{{_SVG_NS}}}g')
    legend.set('transform', f'translate({tx:.1f},{ty:.1f}) scale({sc})')
    for idx, (color, label) in enumerate([
        ('#9A2322', '0.7+'),
        ('#FFCD2F', '0.4-0.7'),
        ('#000000', '0-0.4'),
        ('#B1B3B6', 'No data'),
    ]):
        yo = idx * 13
        r = ET.SubElement(legend, f'{{{_SVG_NS}}}rect')
        r.set('x', '0'); r.set('y', str(yo)); r.set('width', '12'); r.set('height', '12')
        r.set('fill', color); r.set('stroke', '#000000'); r.set('stroke-width', '0.5')
        t = ET.SubElement(legend, f'{{{_SVG_NS}}}text')
        t.set('x', '15'); t.set('y', str(yo + 10))
        t.set('style', 'font-family: Arial, Helvetica, sans-serif; font-size: 11px; fill: #000000;')
        t.text = label

    tree.write(svg_file, xml_declaration=True, encoding='unicode')


if __name__ == "__main__":
    main()
