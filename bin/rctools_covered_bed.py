#!/usr/bin/env python3
# Author: Victoria Begley (begley@ebi.ac.uk)
# License: MIT
"""Build a 4-column BED of covered transcripts and their spliced lengths from a GTF.

Reads transcript IDs (one per line) from covered_ids.txt in the working directory and
computes each transcript's spliced length as the sum of its feature (default exon) lengths
in the GTF. Writes covered.bed as "<id>\\t0\\t<length>\\t<id>" — the 4th column keeps the
clean transcript ID so rf-rctools extract does not rename the region to <id>_0-<end>.

Usage: rctools_covered_bed.py <gtf> <feature_name> <attr_name>
"""
import re
import sys

gtf_path, feature_name, attr_name = sys.argv[1:4]
attr_re = re.compile(r'%s\s+"([^"]+)"' % re.escape(attr_name))

with open("covered_ids.txt") as fh:
    covered = {line.strip() for line in fh if line.strip()}
lengths = {}
with open(gtf_path) as fh:
    for line in fh:
        if not line or line.startswith('#'):
            continue
        cols = line.rstrip('\n').split('\t')
        if len(cols) < 9 or cols[2] != feature_name:
            continue
        m = attr_re.search(cols[8])
        if not m or m.group(1) not in covered:
            continue
        try:
            start = int(cols[3]); end = int(cols[4])
        except ValueError:
            continue
        lengths[m.group(1)] = lengths.get(m.group(1), 0) + (end - start + 1)
with open("covered.bed", "w") as out:
    for tx_id, length in lengths.items():
        if length > 0:
            out.write(f"{tx_id}\t0\t{length}\t{tx_id}\n")
