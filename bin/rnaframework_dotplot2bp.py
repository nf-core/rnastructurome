#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import math
import re
import sys
from pathlib import Path


YEAST_ISOFORM_PATTERN = re.compile(r"^(Y[A-P][LR][0-9]{3}[CW])-([A-Z])(_(?:mRNA|ncRNA|snRNA|snoRNA|rRNA|tRNA))$")
TRANSCRIPT_ID_PATTERN = re.compile(r'transcript_id "([^"]+)"')
COLOR_BINS = [
    (189, 189, 189, "0-10% probability"),
    (242, 204, 84, "10-40% probability"),
    (120, 182, 220, "40-70% probability"),
    (126, 198, 143, "70-100% probability"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert RNAFramework dotplot .dp files to IGV .bp arc files."
    )
    parser.add_argument("--organism", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--fold-dir", required=True)
    parser.add_argument("--gtf", default=None)
    parser.add_argument(
        "--transcript-coords",
        action="store_true",
        help="Output transcript-coordinate arcs (transcript_id as chromosome, raw 1-based positions). "
             "GTF is not required when this flag is set.",
    )
    parser.add_argument(
        "--ucsc-common-chrom-names",
        action="store_true",
        help="Convert common Ensembl chromosome names such as 1, X, MT to chr1, chrX, chrM.",
    )
    parser.add_argument(
        "--min-color-index",
        type=int,
        default=0,
        choices=[0, 1, 2, 3],
        help="Minimum color index to include (0=all, 1=>=10%%, 2=>=40%%, 3=>=70%% probability).",
    )
    args = parser.parse_args()
    if not args.transcript_coords and args.gtf is None:
        parser.error("--gtf is required unless --transcript-coords is set")
    return args


def to_ucsc_common_chrom_name(seqname: str) -> str:
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


def normalize_transcript_id(organism: str, transcript_id: str) -> str:
    if organism == "saccharomyces_cerevisiae":
        return YEAST_ISOFORM_PATTERN.sub(r"\1_\2\3", transcript_id)
    return transcript_id


def transcript_id_candidates(organism: str, transcript_id: str) -> list[str]:
    normalized = normalize_transcript_id(organism, transcript_id)
    candidates = [normalized]
    if "." in normalized:
        candidates.append(normalized.split(".", 1)[0])
    return candidates


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def load_transcripts(organism: str, gtf_path: Path) -> dict[str, dict[str, object]]:
    transcripts: dict[str, dict[str, object]] = {}
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
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            entry = transcripts.setdefault(transcript_id, {"seqname": seqname, "strand": strand, "exons": []})
            entry["exons"].append((start, end))
            if "." in transcript_id:
                transcripts.setdefault(transcript_id.split(".", 1)[0], entry)

    for entry in transcripts.values():
        exons = entry["exons"]
        exons.sort(key=lambda exon: exon[0], reverse=(entry["strand"] == "-"))
    return transcripts


def map_transcript_pos(entry: dict[str, object], position: int) -> int:
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


def color_index_for_probability(score: float) -> int:
    probability = math.pow(10.0, -score)
    if probability >= 0.7:
        return 3
    if probability >= 0.4:
        return 2
    if probability >= 0.1:
        return 1
    return 0


def main() -> int:
    args = parse_args()
    organism = args.organism.strip().lower().replace(" ", "_")
    transcript_coords = args.transcript_coords
    ucsc_chr = args.ucsc_common_chrom_names
    min_color_index = args.min_color_index
    fold_dir = Path(args.fold_dir)
    out_dir = Path(f"{args.prefix}_bp")
    dotplot_out_dir = out_dir / "dotplot"
    warnings_path = out_dir / "conversion_warnings.log"
    dotplot_in_dir = fold_dir / "dotplot"

    dotplot_out_dir.mkdir(parents=True, exist_ok=True)
    warnings: list[str] = []

    if transcript_coords:
        transcripts = {}
    else:
        gtf_path = Path(args.gtf)
        transcripts = load_transcripts(organism, gtf_path)

    dotplot_paths = sorted(dotplot_in_dir.glob("*.dp")) if dotplot_in_dir.is_dir() else []
    bp_count = 0

    for dotplot_path in dotplot_paths:
        transcript_id = normalize_transcript_id(organism, dotplot_path.stem)

        if not transcript_coords:
            entry = None
            matched_transcript_id = None
            for candidate in transcript_id_candidates(organism, transcript_id):
                entry = transcripts.get(candidate)
                if entry:
                    matched_transcript_id = candidate
                    break
            if not entry:
                warnings.append(f"Missing transcript_id '{transcript_id}' in annotation for {dotplot_path.name}; skipping.")
                continue
            if matched_transcript_id != transcript_id:
                warnings.append(
                    f"Matched dotplot transcript_id '{transcript_id}' to annotation transcript_id "
                    f"'{matched_transcript_id}' for {dotplot_path.name}."
                )

        output_path = dotplot_out_dir / f"{transcript_id}.bp"
        with (
            dotplot_path.open("rt", encoding="utf-8") as reader,
            output_path.open("wt", encoding="utf-8") as writer,
        ):
            first_line = reader.readline()
            header_line = reader.readline()
            if not first_line or not header_line:
                warnings.append(f"Malformed dotplot file {dotplot_path.name}; skipping.")
                output_path.unlink(missing_ok=True)
                continue
            for red, green, blue, label in COLOR_BINS:
                writer.write(f"color:\t{red}\t{green}\t{blue}\t{label}\n")
            arcs = []
            for raw_line in reader:
                line = raw_line.strip()
                if not line:
                    continue
                fields = line.split("\t")
                if len(fields) < 3:
                    continue
                left_pos = int(fields[0])
                right_pos = int(fields[1])
                color_index = color_index_for_probability(float(fields[2]))
                if color_index < min_color_index:
                    continue

                if transcript_coords:
                    start, end = sorted((left_pos, right_pos))
                    seqname = transcript_id
                else:
                    try:
                        left_genome = map_transcript_pos(entry, left_pos)
                        right_genome = map_transcript_pos(entry, right_pos)
                    except ValueError as exc:
                        warnings.append(
                            f"Skipping base-pair record in {dotplot_path.name}: {exc}."
                        )
                        continue
                    start, end = sorted((left_genome, right_genome))
                    seqname = entry["seqname"]
                    if not seqname or seqname.lower() == "none":
                        warnings.append(
                            f"Skipping base-pair record in {dotplot_path.name}: "
                            f"empty or null seqname for transcript '{transcript_id}'."
                        )
                        continue
                    if ucsc_chr:
                        seqname = to_ucsc_common_chrom_name(seqname)

                arcs.append((seqname, start, end, color_index))

            # rf-fold writes dot-plot pairs in no fixed order; sort so reruns give identical files.
            for seqname, start, end, color_index in sorted(arcs):
                writer.write(f"{seqname}\t{start}\t{start}\t{end}\t{end}\t{color_index}\n")

        if arcs:
            bp_count += 1
        else:
            output_path.unlink(missing_ok=True)
            warnings.append(f"No convertible base-pair records found in {dotplot_path.name}; skipping.")

    if warnings:
        warnings_path.write_text("\n".join(warnings) + "\n", encoding="utf-8")
    elif warnings_path.exists():
        warnings_path.unlink()

    if dotplot_paths and bp_count == 0:
        print(
            "WARNING: No .bp files were generated from the available .dp files. "
            "Check conversion_warnings.log for per-transcript details.",
            file=sys.stderr,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
