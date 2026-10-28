#!/usr/bin/env python
# Author: Victoria Begley (begley@ebi.ac.uk)
# License: MIT
"""Window rf-norm XML reactivities to user-supplied reference-structure regions.

rf-eval compares a reference structure to reactivities position-by-position over
the whole transcript, so a sub-region structure (e.g. an Rfam element on a whole
chromosome) can only be scored fairly against a sliced XML. For each row in the
windows manifest this cuts the matching per-transcript XML to [start, end] and
writes <structure_id>.xml, whose id matches the structure entry in --rfeval_reference.

Windows manifest: whitespace/tab-separated, '#' comments ignored, columns:
    ref_seq_id  start  end  structure_id  [strand]
ref_seq_id matches the rf-norm XML filename stem (= transcript/reference id).
Coordinates are 1-based inclusive, in the XML's own sequence space (genome
coordinates on the genome route, transcript coordinates on the transcriptome
route). strand is optional ('+' default); '-' reverse-complements the window.
"""
import argparse
import glob
import os
import re
import sys

COMPLEMENT = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")


def parse_windows(path):
    rows = []
    with open(path, encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            parts = line.split()
            if parts[0].lower() in ("ref_seq_id", "chrom", "chromosome"):
                continue  # header row
            if len(parts) < 4:
                sys.exit(f"windows line {lineno}: need >=4 columns, got {len(parts)}")
            ref, start, end, sid = parts[0], parts[1], parts[2], parts[3]
            strand = parts[4] if len(parts) > 4 else "+"
            rows.append((ref, int(start), int(end), sid, strand))
    if not rows:
        sys.exit(f"no usable rows in {path}")
    return rows


def load_xml(path):
    text = open(path, encoding="utf-8").read()
    data = re.search(r"<data[^>]*>", text).group(0)
    seq = "".join(re.search(r"<sequence>(.*?)</sequence>", text, re.S).group(1).split())
    rb = re.search(r"<reactivity>(.*?)</reactivity>", text, re.S).group(1)
    react = [v for v in re.split(r"[,\s]+", rb.strip()) if v != ""]
    return data, seq, react


def load_ref(path):
    """Read and validate one reference XML once; returns (payload, reason-unusable)."""
    data, seq, react = load_xml(path)
    if len(react) != len(seq):
        return None, f"{len(react)} reactivities vs {len(seq)} nt"
    return (data, seq, react), None


def write_window(out, data, sid, seq, react, strand):
    if strand == "-":
        seq = seq.translate(COMPLEMENT)[::-1]
        react = react[::-1]
    length = len(seq)
    with open(out, "w", encoding="utf-8") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fh.write(data + "\n")
        fh.write(f'\t<transcript id="{sid}" length="{length}">\n')
        fh.write("\t\t<sequence>" + seq + "</sequence>\n")
        fh.write("\t\t<reactivity>" + ",".join(react) + "</reactivity>\n")
        fh.write("\t</transcript>\n</data>\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--windows", required=True, help="coordinate manifest")
    ap.add_argument("--xml-glob", default="xml_input*/*.xml",
                    help="glob for staged rf-norm XMLs (default: %(default)s)")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    # index staged XMLs by filename stem (rf-norm names each file by transcript id)
    xmls = {os.path.splitext(os.path.basename(p))[0]: p
            for p in glob.glob(args.xml_glob)}
    if not xmls:
        sys.exit(f"no XMLs matched {args.xml_glob}")

    os.makedirs(args.outdir, exist_ok=True)
    written, skipped = 0, []
    refs, bad = {}, {}  # ref -> payload or None; unusable ref -> [reason, dropped]
    for ref, start, end, sid, strand in parse_windows(args.windows):
        if ref not in refs:
            # each reference is read and length-checked once, however many windows use it
            if ref in xmls:
                refs[ref], reason = load_ref(xmls[ref])
            else:
                refs[ref], reason = None, "no XML found"
            if reason:
                bad[ref] = [reason, 0]
        if refs[ref] is None:
            bad[ref][1] += 1
            continue
        data, seq, react = refs[ref]
        if start < 1 or end > len(seq) or start > end:
            skipped.append(f"{sid}: window {start}-{end} outside {ref} (len {len(seq)})")
            continue
        write_window(os.path.join(args.outdir, f"{sid}.xml"), data, sid,
                     seq[start - 1:end], react[start - 1:end], strand)
        written += 1

    for ref, (reason, dropped) in bad.items():
        print(f"[rfeval-window] skip reference '{ref}': {reason} "
              f"({dropped} window(s) dropped)", file=sys.stderr)
    for msg in skipped:
        print(f"[rfeval-window] skip {msg}", file=sys.stderr)
    if not written:
        sys.exit("no windows written; check that ref_seq_id values match XML names")
    print(f"[rfeval-window] wrote {written} window(s) to {args.outdir}")


if __name__ == "__main__":
    main()
