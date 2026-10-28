#!/usr/bin/env python3
# Author: Victoria Begley (begley@ebi.ac.uk)
# License: MIT
"""Remap RNAframework transcript-coordinate WIG tracks to genomic WIG tracks.

RNAframework emits reactivity and Shannon entropy WIG files keyed by transcript
ID.  This script uses exon coordinates from a GTF to project those per-base
values onto genomic coordinates, averaging values when multiple transcripts
map to the same genomic base.
"""

from __future__ import annotations

import argparse
import gzip
import math
import re
import sys
from collections import defaultdict
from collections.abc import Iterator
from pathlib import Path
from typing import TextIO
from typing import TypedDict


# Keep this in sync with rnaframework_dotplot2bp.py: yeast RNAframework
# In yeast transcript IDs can use "-" where the GTF uses "_" for isoform suffixes.
YEAST_ISOFORM_PATTERN = re.compile(
    r"^(Y[A-P][LR][0-9]{3}[CW])-([A-Z])(_(?:mRNA|ncRNA|snRNA|snoRNA|rRNA|tRNA))$"
)
TRANSCRIPT_ID_PATTERN = re.compile(r'transcript_id "([^"]+)"')

Exon = tuple[int, int]


class TranscriptEntry(TypedDict):
    """Transcript metadata needed to map transcript positions to genome positions."""

    seqname: str
    strand: str
    exons: list[Exon]


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the remapping script."""

    parser = argparse.ArgumentParser(
        description="Remap transcript-coordinate WIG values to genomic-coordinate WIG values."
    )
    parser.add_argument("--wig", required=True, help="Input transcript-coordinate WIG file.")
    parser.add_argument(
        "--gtf",
        required=True,
        help="GTF annotation with transcript exon coordinates.",
    )
    parser.add_argument("--output-wig", required=True, help="Output genomic-coordinate WIG file.")
    parser.add_argument(
        "--chrom-sizes",
        required=True,
        help="Output genomic chromosome sizes file.",
    )
    parser.add_argument(
        "--organism",
        default="",
        help="Organism name used for transcript ID normalisation.",
    )
    parser.add_argument(
        "--ucsc-common-chrom-names",
        action="store_true",
        help="Convert common Ensembl chromosome names such as 1, X, MT to chr1, chrX, chrM.",
    )
    return parser.parse_args()


def normalize_transcript_id(organism: str, transcript_id: str) -> str:
    """Normalize organism-specific transcript IDs before GTF lookup."""

    if organism == "saccharomyces_cerevisiae":
        return YEAST_ISOFORM_PATTERN.sub(r"\1_\2\3", transcript_id)
    return transcript_id


def transcript_id_candidates(organism: str, transcript_id: str) -> list[str]:
    """Return possible GTF transcript IDs for an observed WIG transcript ID."""

    normalized = normalize_transcript_id(organism, transcript_id)
    candidates = [normalized]
    if "." in normalized:
        candidates.append(normalized.split(".", 1)[0])
    return candidates


def open_text(path: Path) -> TextIO:
    """Open plain-text or gzip-compressed text files."""

    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def to_ucsc_common_chrom_name(seqname: str) -> str:
    """Convert common Ensembl chromosome names to UCSC names used by IGV."""

    if seqname.startswith("chr"):
        return seqname
    if seqname == "MT":
        return "chrM"
    if seqname == "M":
        return "chrM"
    if seqname in {"X", "Y"}:
        return f"chr{seqname}"
    if seqname.isdigit():
        return f"chr{seqname}"
    return seqname


def chrom_sort_key(seqname: str) -> tuple[int, int | str]:
    """Return a stable chromosome sort key with primary chromosomes first."""

    chrom = seqname[3:] if seqname.startswith("chr") else seqname
    if chrom == "M":
        chrom = "MT"
    if chrom.isdigit():
        return (0, int(chrom))
    if chrom == "X":
        return (1, 23)
    if chrom == "Y":
        return (1, 24)
    if chrom == "MT":
        return (1, 25)
    return (2, chrom)


def load_transcripts(
    organism: str,
    gtf_path: Path,
    ucsc_common_chrom_names: bool,
) -> tuple[dict[str, TranscriptEntry], dict[str, int]]:
    """Load exon models and inferred chromosome sizes from a GTF file."""

    transcripts: dict[str, TranscriptEntry] = {}
    chrom_sizes: dict[str, int] = {}

    with open_text(gtf_path) as handle:
        for raw_line in handle:
            if not raw_line or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            match = TRANSCRIPT_ID_PATTERN.search(fields[8])
            if not match:
                continue
            transcript_id = normalize_transcript_id(organism, match.group(1))
            seqname = fields[0]
            if ucsc_common_chrom_names:
                seqname = to_ucsc_common_chrom_name(seqname)
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]

            chrom_sizes[seqname] = max(chrom_sizes.get(seqname, 0), end)
            if transcript_id not in transcripts:
                transcripts[transcript_id] = {"seqname": seqname, "strand": strand, "exons": []}
            entry = transcripts[transcript_id]
            entry["exons"].append((start, end))
            if "." in transcript_id:
                transcripts.setdefault(transcript_id.split(".", 1)[0], entry)

    for entry in transcripts.values():
        exons = entry["exons"]
        exons.sort(key=lambda exon: exon[0], reverse=entry["strand"] == "-")

    return transcripts, chrom_sizes


def map_transcript_pos(entry: TranscriptEntry, position: int) -> int:
    """Map a 1-based transcript coordinate to a 1-based genomic coordinate."""

    cursor = 1
    for start, end in entry["exons"]:
        exon_len = end - start + 1
        if cursor <= position < cursor + exon_len:
            offset = position - cursor
            if entry["strand"] == "-":
                return end - offset
            return start + offset
        cursor += exon_len
    raise ValueError(f"Transcript position {position} exceeds transcript length")


def wig_records(wig_path: Path) -> Iterator[tuple[str, int, int, float]]:
    """Yield transcript ID, 1-based position, span, and value from a WIG file."""

    chrom = None
    is_fixed = False
    fixed_start = 1
    fixed_step = 1
    span = 1
    pos_counter = 0

    with open_text(wig_path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("track"):
                continue
            if line.startswith("variableStep"):
                parts = dict(part.split("=", 1) for part in line.split() if "=" in part)
                chrom = parts.get("chrom")
                span = int(parts.get("span", 1))
                is_fixed = False
                continue
            if line.startswith("fixedStep"):
                parts = dict(part.split("=", 1) for part in line.split() if "=" in part)
                chrom = parts.get("chrom")
                fixed_start = int(parts.get("start", 1))
                fixed_step = int(parts.get("step", 1))
                span = int(parts.get("span", 1))
                pos_counter = 0
                is_fixed = True
                continue
            if chrom is None:
                raise ValueError(f"Encountered WIG data line before chrom declaration: {line!r}")
            if is_fixed:
                pos = fixed_start + pos_counter * fixed_step
                value = float(line)
                pos_counter += 1
            else:
                pos_str, value_str = line.split()[:2]
                pos = int(pos_str)
                value = float(value_str)
            yield chrom, pos, span, value


def format_value(value: float) -> str:
    """Format WIG numeric values compactly without unnecessary trailing digits."""

    return f"{value:.6g}"


def remap_wig(
    organism: str,
    wig_path: Path,
    gtf_path: Path,
    output_wig_path: Path,
    chrom_sizes_path: Path,
    ucsc_common_chrom_names: bool,
) -> int:
    """Remap transcript WIG values into genomic coordinates and write WIG outputs."""

    transcripts, chrom_sizes = load_transcripts(organism, gtf_path, ucsc_common_chrom_names)
    if not transcripts:
        raise ValueError(f"No exon transcript records found in {gtf_path}")

    genomic_values: dict[str, dict[int, list[float]]] = defaultdict(
        lambda: defaultdict(lambda: [0.0, 0.0])
    )
    missing_transcripts: set[str] = set()
    skipped_positions = 0

    for transcript_id, pos, span, value in wig_records(wig_path):
        if not math.isfinite(value):
            continue
        entry = None
        for candidate in transcript_id_candidates(organism, transcript_id):
            entry = transcripts.get(candidate)
            if entry:
                break
        if not entry:
            missing_transcripts.add(transcript_id)
            continue

        seqname = entry["seqname"]
        for transcript_pos in range(pos, pos + span):
            try:
                genome_pos = map_transcript_pos(entry, transcript_pos)
            except ValueError:
                skipped_positions += 1
                continue
            bucket = genomic_values[seqname][genome_pos]
            bucket[0] += value
            bucket[1] += 1.0

    for transcript_id in sorted(missing_transcripts):
        print(
            f"WARNING: missing transcript_id '{transcript_id}' in annotation; skipped.",
            file=sys.stderr,
        )
    if skipped_positions:
        print(
            "WARNING: skipped "
            f"{skipped_positions} transcript position(s) outside annotated exon length.",
            file=sys.stderr,
        )
    if not genomic_values:
        print(
            "WARNING: No genomic WIG records were generated "
            "(all reactivity values are NaN or zero). "
            "Writing empty output files.",
            file=sys.stderr,
        )
        chrom_sizes_path.open("wt", encoding="utf-8").close()
        with output_wig_path.open("wt", encoding="utf-8") as handle:
            handle.write("track type=wiggle_0\n")
        return 0

    with chrom_sizes_path.open("wt", encoding="utf-8") as handle:
        for seqname in sorted(chrom_sizes, key=chrom_sort_key):
            handle.write(f"{seqname}\t{chrom_sizes[seqname]}\n")

    with output_wig_path.open("wt", encoding="utf-8") as handle:
        handle.write("track type=wiggle_0\n")
        for seqname in sorted(genomic_values, key=chrom_sort_key):
            handle.write(f"variableStep chrom={seqname}\n")
            for genome_pos in sorted(genomic_values[seqname]):
                total, count = genomic_values[seqname][genome_pos]
                handle.write(f"{genome_pos} {format_value(total / count)}\n")

    return 0


def main() -> int:
    """Run the command-line remapping workflow."""

    args = parse_args()
    organism = args.organism.strip().lower().replace(" ", "_")
    try:
        return remap_wig(
            organism=organism,
            wig_path=Path(args.wig),
            gtf_path=Path(args.gtf),
            output_wig_path=Path(args.output_wig),
            chrom_sizes_path=Path(args.chrom_sizes),
            ucsc_common_chrom_names=args.ucsc_common_chrom_names,
        )
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
