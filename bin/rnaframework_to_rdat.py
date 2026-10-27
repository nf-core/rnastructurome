#!/usr/bin/env python3
# Author: Victoria Begley (begley@ebi.ac.uk)
# License: MIT
"""Convert RNAframework XML + rf-fold .db files to RDAT format."""
from __future__ import annotations

import argparse
import math
import re
import statistics
import sys
from dataclasses import dataclass, field
from pathlib import Path
from xml.etree import ElementTree


SCORING_LABELS: dict[str, str] = {
    "1": "Ding",
    "2": "Rouskin",
    "3": "Siegfried",
    "4": "Zubradt",
}

NORM_LABELS: dict[str, str] = {
    "2": "90% Winsorizing",
    "3": "Box-plot",
}


@dataclass
class RdatRecord:
    """All data needed to write one RDAT entry."""

    transcript_id: str
    rna_sequence: str
    dot_bracket: str
    reactivities: list[float | None]
    reactivity_errors: list[float | None] | None = None
    comments: list[str] = field(default_factory=list)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--xml-dir", required=True, help="Directory containing rf-norm XML files"
    )
    parser.add_argument(
        "--structures-dir", required=True, help="Directory containing rf-fold .db files"
    )
    parser.add_argument(
        "--prefix", required=True, help="Output prefix; files written to <prefix>_rdat/"
    )
    parser.add_argument(
        "--fasta", default=None, help="Transcript FASTA used for alignment and counting"
    )
    parser.add_argument("--pipeline-version", default=None, help="Pipeline version string")
    parser.add_argument(
        "--principle", default=None, help="Probing principle (RT-stop or MaP)"
    )
    parser.add_argument(
        "--rfnorm-scoring-method", default=None, help="rf-norm scoring method code (1-4)"
    )
    parser.add_argument(
        "--rfnorm-norm-method", default=None, help="rf-norm normalisation method code (2-3)"
    )
    parser.add_argument(
        "--gtf", default=None, help="GTF annotation file used for coordinate mapping"
    )
    parser.add_argument(
        "--dotplot-dir",
        default=None,
        help="Directory containing rf-fold .dp files; writes <id>.bpp.tsv beside each RDAT",
    )
    return parser.parse_args()


def parse_xml(xml_path: Path) -> tuple[str, str, list[float | None]]:
    """Return (transcript_id, dna_sequence, reactivities).

    Reactivities are per-base floats; None indicates a missing/NaN value.
    """
    tree = ElementTree.parse(xml_path)
    root = tree.getroot()

    transcript_el = root.find("transcript")
    if transcript_el is None:
        raise ValueError(f"No <transcript> element found in {xml_path}")

    transcript_id = transcript_el.get("id", xml_path.stem)

    seq_el = transcript_el.find("sequence")
    dna_sequence = "".join((seq_el.text or "").split()) if seq_el is not None else ""

    react_el = transcript_el.find("reactivity")
    reactivities: list[float | None] = []
    if react_el is not None and react_el.text:
        for token in react_el.text.split(","):
            token = token.strip()
            if not token or token.lower() == "nan":
                reactivities.append(None)
            else:
                try:
                    reactivities.append(float(token))
                except ValueError:
                    reactivities.append(None)

    return transcript_id, dna_sequence, reactivities


def parse_db(db_path: Path) -> tuple[str, str]:
    """Return (rna_sequence, dot_bracket) from an rf-fold .db file.

    .db format:
        Line 1: >transcript_id
        Line 2: RNA sequence (may use T or U)
        Line 3: dot-bracket string, optionally followed by whitespace and MFE score
    """
    raw = db_path.read_text(encoding="utf-8").splitlines()
    lines = [line.rstrip("\n") for line in raw if line.strip()]
    if len(lines) < 3:
        raise ValueError(f"Malformed .db file (fewer than 3 non-empty lines): {db_path}")

    rna_sequence = lines[1].strip().upper().replace("T", "U")

    # Strip optional trailing MFE score like " (-770.27)"
    dot_bracket = re.split(r"\s", lines[2].strip())[0]

    return rna_sequence, dot_bracket


def parse_dotplot(dp_path: Path) -> dict[tuple[int, int], float]:
    """Return {(i, j): probability} from an rf-fold .dp file (1-based, -log10 probabilities)."""
    probabilities: dict[tuple[int, int], float] = {}
    for line in dp_path.read_text(encoding="utf-8").splitlines()[2:]:
        fields = line.split()
        if len(fields) == 3:
            probabilities[(int(fields[0]), int(fields[1]))] = 10 ** -float(fields[2])
    return probabilities


def structure_pairs(dot_bracket: str) -> list[tuple[int, int]]:
    """Return the 1-based (i, j) base pairs of a dot-bracket string, pseudoknot brackets included."""
    closers = {")": "(", "]": "[", "}": "{", ">": "<"}
    stacks: dict[str, list[int]] = {opener: [] for opener in closers.values()}
    pairs: list[tuple[int, int]] = []
    for pos, char in enumerate(dot_bracket, start=1):
        if char in stacks:
            stacks[char].append(pos)
        elif char in closers:
            pairs.append((stacks[closers[char]].pop(), pos))
    return sorted(pairs)


def write_bpp(out_path: Path, dot_bracket: str, probabilities: dict[tuple[int, int], float]) -> None:
    """Write one row per structure base pair with its rf-fold pairing probability."""
    with out_path.open("wt", encoding="utf-8") as fh:
        fh.write("i\tj\tprobability\n")
        for i, j in structure_pairs(dot_bracket):
            # rf-fold leaves very unlikely pairs out of the dot plot, so their probability is unknown.
            probability = probabilities.get((i, j))
            fh.write(f"{i}\t{j}\t{'NaN' if probability is None else f'{probability:.6g}'}\n")


def resolve_db_path(
    structures_dir: Path, transcript_id: str, xml_stems: list[str]
) -> Path | None:
    """Return the first matching .db file for transcript_id, or None."""
    candidates: list[str] = [transcript_id, transcript_id.rsplit(".", 1)[0]]
    for stem in xml_stems:
        candidates.extend([stem, stem.rsplit(".", 1)[0]])

    seen: set[str] = set()
    for candidate in candidates:
        if not candidate or candidate in seen:
            continue
        seen.add(candidate)
        db_path = structures_dir / f"{candidate}.db"
        if db_path.exists():
            return db_path

    return None


def aggregate_reactivities(
    reactivity_sets: list[list[float | None]],
    length: int,
) -> tuple[list[float | None], list[float | None] | None]:
    """Compute per-position mean (and SEM when >1 replicate) across reactivity sets."""
    aggregated: list[float | None] = []
    errors: list[float | None] | None = [] if len(reactivity_sets) > 1 else None

    for idx in range(length):
        values = [
            reactivities[idx]
            for reactivities in reactivity_sets
            if idx < len(reactivities) and reactivities[idx] is not None
        ]
        aggregated.append(statistics.fmean(values) if values else None)

        if errors is not None:
            if len(values) >= 2:
                errors.append(statistics.stdev(values) / math.sqrt(len(values)))
            else:
                errors.append(None)

    return aggregated, errors


def write_rdat(out_path: Path, record: RdatRecord) -> None:
    """Write a single RDAT file from a RdatRecord."""
    length = len(record.rna_sequence)

    def fmt(values: list[float | None]) -> str:
        return " ".join("NaN" if v is None else f"{v:.6g}" for v in values)

    padded_reactivities = list(record.reactivities) + [None] * max(
        0, length - len(record.reactivities)
    )

    with out_path.open("wt", encoding="utf-8") as fh:
        fh.write(f"NAME\t{record.transcript_id}\n")
        fh.write(f"SEQUENCE\t{record.rna_sequence}\n")
        fh.write(f"STRUCTURE\t{record.dot_bracket}\n")
        fh.write("ANNOTATION_DATA:1\tmodifier:DMS\n")
        fh.write(f"REACTIVITY\t{fmt(padded_reactivities[:length])}\n")
        if record.reactivity_errors is not None:
            padded_errors = list(record.reactivity_errors) + [None] * max(
                0, length - len(record.reactivity_errors)
            )
            fh.write(f"REACTIVITY_ERROR\t{fmt(padded_errors[:length])}\n")
        for comment in record.comments:
            fh.write(f"COMMENT\t{comment}\n")


def build_comment(args: argparse.Namespace) -> str | None:
    """Build a single COMMENT string from pipeline metadata args, or None if empty."""
    parts: list[str] = []
    if args.pipeline_version:
        parts.append(f"Generated by nf-core/rnastructurome v{args.pipeline_version}")
    if args.fasta:
        parts.append(f"FASTA: {args.fasta}")
    if args.gtf:
        parts.append(f"GTF: {args.gtf}")
    if args.principle:
        parts.append(f"Principle: {args.principle}")
    if args.rfnorm_scoring_method:
        sm = str(args.rfnorm_scoring_method)
        parts.append(f"rf-norm scoring: sm={sm} ({SCORING_LABELS.get(sm, sm)})")
    if args.rfnorm_norm_method:
        nm = str(args.rfnorm_norm_method)
        parts.append(f"normalisation: nm={nm} ({NORM_LABELS.get(nm, nm)})")
    return "; ".join(parts) if parts else None


def main() -> int:
    """Entry point."""
    args = parse_args()
    xml_dir = Path(args.xml_dir)
    structures_dir = Path(args.structures_dir)
    out_dir = Path(f"{args.prefix}_rdat")
    out_dir.mkdir(parents=True, exist_ok=True)

    comment = build_comment(args)
    comments = [comment] if comment else []

    xml_files = sorted(path for path in xml_dir.rglob("*.xml") if path.is_file())
    if not xml_files:
        print(f"No XML files found in {xml_dir}", file=sys.stderr)
        return 1

    written = 0
    skipped = 0
    transcript_reactivities: dict[str, list[list[float | None]]] = {}
    transcript_xml_stems: dict[str, list[str]] = {}

    for xml_path in xml_files:
        try:
            transcript_id, _dna_sequence, reactivities = parse_xml(xml_path)
        except Exception as exc:  # noqa: BLE001
            print(f"WARNING: Could not parse {xml_path.name}: {exc}", file=sys.stderr)
            skipped += 1
            continue

        transcript_reactivities.setdefault(transcript_id, []).append(reactivities)
        transcript_xml_stems.setdefault(transcript_id, []).append(xml_path.stem)

    for transcript_id in sorted(transcript_reactivities):
        xml_stems = transcript_xml_stems.get(transcript_id, [])
        db_path = resolve_db_path(structures_dir, transcript_id, xml_stems)
        if db_path is None:
            print(
                f"WARNING: No .db file found for {transcript_id} in {structures_dir}; skipping.",
                file=sys.stderr,
            )
            skipped += 1
            continue

        try:
            rna_sequence, dot_bracket = parse_db(db_path)
        except Exception as exc:  # noqa: BLE001
            print(f"WARNING: Could not parse {db_path.name}: {exc}", file=sys.stderr)
            skipped += 1
            continue

        aggregated_reactivities, reactivity_errors = aggregate_reactivities(
            transcript_reactivities[transcript_id],
            len(rna_sequence),
        )
        out_path = out_dir / f"{transcript_id}.rdat"
        try:
            write_rdat(
                out_path,
                RdatRecord(
                    transcript_id=transcript_id,
                    rna_sequence=rna_sequence,
                    dot_bracket=dot_bracket,
                    reactivities=aggregated_reactivities,
                    reactivity_errors=reactivity_errors,
                    comments=comments,
                ),
            )
            written += 1
        except Exception as exc:  # noqa: BLE001
            print(f"WARNING: Could not write RDAT for {transcript_id}: {exc}", file=sys.stderr)
            skipped += 1
            continue

        dp_path = Path(args.dotplot_dir) / f"{db_path.stem}.dp" if args.dotplot_dir else None
        if dp_path is not None and dp_path.exists():
            try:
                write_bpp(out_dir / f"{transcript_id}.bpp.tsv", dot_bracket, parse_dotplot(dp_path))
            except Exception as exc:  # noqa: BLE001
                print(f"WARNING: Could not write base-pair probabilities for {transcript_id}: {exc}", file=sys.stderr)

    print(f"Wrote {written} RDAT file(s), skipped {skipped}.", file=sys.stderr)

    if written == 0 and xml_files:
        print("ERROR: No RDAT files were produced.", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
