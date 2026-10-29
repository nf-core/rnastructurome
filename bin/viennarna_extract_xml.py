#!/usr/bin/env python3
# Author: Victoria Begley (begley@ebi.ac.uk)
# License: MIT
"""Extract per-position reactivities for a transcript from RNAframework XML files.

Averages values across replicates; positions with no data across all replicates are
output as -999 (the RNAplot sentinel for "no data").

Usage: viennarna_extract_xml.py <transcript_id> <xml_file> [<xml_file> ...]
Output: tab-separated <position>\\t<reactivity> lines written to stdout.
"""
from __future__ import annotations

import math
import sys
import xml.etree.ElementTree as ET


def main():
    """Extract and average per-position reactivities from XML files, writing to stdout."""
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <transcript_id> <xml_file> [...]", file=sys.stderr)
        sys.exit(1)

    tid = sys.argv[1]
    reps: list[list[float]] = []

    for xml_path in sys.argv[2:]:
        try:
            root = ET.parse(xml_path).getroot()
        except (OSError, ET.ParseError):
            continue
        for t in root.iter('transcript'):
            if (t.get('id') or t.findtext('id')) != tid:
                continue
            raw = t.findtext('reactivity') or t.findtext('values') or ''
            vals: list[float] = []
            for v in raw.strip().split(','):
                v = v.strip()
                try:
                    vals.append(float(v))
                except ValueError:
                    vals.append(float('nan'))
            if vals:
                reps.append(vals)
            break

    if not reps:
        sys.exit(0)

    length = max(len(r) for r in reps)
    for i in range(length):
        col = [r[i] for r in reps if i < len(r) and not math.isnan(r[i])]
        v = sum(col) / len(col) if col else float('nan')
        out = -999 if math.isnan(v) or v < 0 else v
        print(f"{i + 1}\t{out}")


if __name__ == "__main__":
    main()
