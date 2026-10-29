#!/usr/bin/env python
# Author: Victoria Begley (begley@ebi.ac.uk)
# License: MIT
"""Build a circular-rotation baseline reactivity set for rf-eval.

Each decoy keeps a structure and its reactivity values but rotates the profile
against the sequence, so scoring them gives what that structure scores by chance.
Rotation beats shuffling: reactivity is autocorrelated, and a shuffle gives a
baseline ~1.1-1.5x too narrow. Written as <id>__rot<NNNN>.xml plus a decoy .db.
"""
import argparse
import glob
import os
import re
import sys


def load_db(path):
    """Parse a dot-bracket .db into {id: (sequence, structure)}."""
    entries, current = {}, None
    for line in open(path, encoding="utf-8"):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            current = line[1:].split()[0]
            entries[current] = []
        elif current is not None:
            entries[current].append(line)
    bad = [k for k, v in entries.items() if len(v) < 2]
    if bad:
        sys.exit(f"{path}: entries missing sequence or structure: {', '.join(bad)}")
    return {k: (v[0], v[1]) for k, v in entries.items()}


def load_xml(path):
    text = open(path, encoding="utf-8").read()
    header = re.search(r"<data[^>]*>", text).group(0)
    seq = "".join(re.search(r"<sequence>(.*?)</sequence>", text, re.S).group(1).split())
    body = re.search(r"<reactivity>(.*?)</reactivity>", text, re.S).group(1)
    return header, seq, [v for v in re.split(r"[,\s]+", body.strip()) if v]


def write_xml(path, header, sid, seq, react):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fh.write(header + "\n")
        fh.write(f'\t<transcript id="{sid}" length="{len(seq)}">\n')
        fh.write("\t\t<sequence>" + seq + "</sequence>\n")
        fh.write("\t\t<reactivity>" + ",".join(react) + "</reactivity>\n")
        fh.write("\t</transcript>\n</data>\n")


def pick_shifts(length, min_shift, limit):
    """All rotations far enough from the identity, thinned evenly to <= limit."""
    shifts = [s for s in range(1, length) if min(s, length - s) >= min_shift]
    if len(shifts) <= limit:
        return shifts
    step = len(shifts) / float(limit)
    return [shifts[int(i * step)] for i in range(limit)]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--xml-glob", required=True, help="glob for per-structure reactivity XMLs")
    ap.add_argument("--db", required=True, help="reference .db matching those XMLs")
    ap.add_argument("--outdir", required=True, help="dir to write decoy XMLs into")
    ap.add_argument("--decoy-db", required=True, help="path for the decoy .db")
    ap.add_argument("-k", "--replicates", type=int, default=200,
                    help="max decoys per structure (default: %(default)s)")
    ap.add_argument("--min-shift", type=int, default=5,
                    help="smallest rotation treated as a null (default: %(default)s)")
    args = ap.parse_args()

    if args.replicates < 1:
        sys.exit("--replicates must be >= 1")

    ref = load_db(args.db)
    os.makedirs(args.outdir, exist_ok=True)
    entries, written, skipped = [], 0, []

    for path in sorted(glob.glob(args.xml_glob)):
        sid = os.path.splitext(os.path.basename(path))[0]
        if sid not in ref:
            skipped.append(f"{sid}: no matching structure in {os.path.basename(args.db)}")
            continue
        header, seq, react = load_xml(path)
        structure = ref[sid][1]
        if len(structure) != len(react):
            skipped.append(f"{sid}: structure {len(structure)} nt vs {len(react)} reactivities")
            continue
        shifts = pick_shifts(len(react), args.min_shift, args.replicates)
        if not shifts:
            skipped.append(f"{sid}: {len(react)} nt too short to rotate")
            continue
        for i, shift in enumerate(shifts):
            did = f"{sid}__rot{i:04d}"
            write_xml(os.path.join(args.outdir, did + ".xml"), header, did, seq,
                      react[-shift:] + react[:-shift])
            entries.append(f">{did}\n{seq}\n{structure}")
            written += 1

    for msg in skipped:
        print(f"[rfeval-baseline] skip {msg}", file=sys.stderr)
    if not written:
        sys.exit("no decoys written; check that XML names match structure ids")

    with open(args.decoy_db, "w", encoding="utf-8") as fh:
        fh.write("\n".join(entries) + "\n")
    print(f"[rfeval-baseline] wrote {written} decoy(s) to {args.outdir}")


if __name__ == "__main__":
    main()
