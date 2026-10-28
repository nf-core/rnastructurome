#!/usr/bin/env python
# Author: Victoria Begley (begley@ebi.ac.uk)
# License: MIT
"""Merge rf-eval metrics with its rotation baseline into one publishable table.

Only AUROC has a fixed chance level (0.5); DSCI and the unpaired coefficient
depend on each structure's paired/unpaired ratio and helix layout. So this adds
the mean and sd across the decoy scores, per structure and per metric: read the
gap between a score and its own baseline mean, in units of the baseline sd.
"""
import argparse
import collections
import math
import os
import re
import sys

DECOY_ID = re.compile(r"^(?P<sid>.+)__rot\d+$")
COLUMNS = ("coeff_unpaired", "dsci", "auroc")


def read_metrics(path):
    """Parse an rf-eval metrics.txt into {id: {column: value}}."""
    rows = {}
    with open(path, encoding="utf-8") as fh:
        header = fh.readline()
        if not header.strip():
            return rows
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            try:
                values = [float(v) for v in fields[1:4]]
            except ValueError:
                continue
            rows[fields[0]] = dict(zip(COLUMNS, values))
    return rows


def mean_sd(values):
    n = len(values)
    mean = sum(values) / n
    sd = math.sqrt(sum((v - mean) ** 2 for v in values) / (n - 1)) if n > 1 else 0.0
    return mean, sd


def fmt(value):
    return "" if value is None else f"{value:.4f}"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--metrics", required=True, help="rf-eval metrics.txt")
    ap.add_argument("--baseline-metrics", help="metrics.txt from the decoy run")
    ap.add_argument("--output", required=True, help="output TSV")
    args = ap.parse_args()

    observed = read_metrics(args.metrics)
    if not observed:
        sys.exit(f"no usable rows in {args.metrics}")

    baselines = collections.defaultdict(lambda: collections.defaultdict(list))
    if args.baseline_metrics and os.path.exists(args.baseline_metrics):
        for did, row in read_metrics(args.baseline_metrics).items():
            match = DECOY_ID.match(did)
            if not match:
                continue
            for column in COLUMNS:
                baselines[match.group("sid")][column].append(row[column])

    # each metric sits next to its own baseline, so a score is read against the
    # right pair of columns rather than one block away
    header = ["structure"]
    for column in COLUMNS:
        header.append(column)
        if baselines:
            header += [f"{column}_baseline_mean", f"{column}_baseline_sd"]
    if baselines:
        header.append("baseline_num")

    # rf-eval's "Overall" row pools every structure into one score. Structures have
    # different baselines, so a pooled number is not interpretable — drop it.
    ordered = sorted(k for k in observed if k != "Overall")

    with open(args.output, "w", encoding="utf-8") as fh:
        fh.write("\t".join(header) + "\n")
        for sid in ordered:
            decoys = baselines.get(sid) if baselines else None
            row = [sid]
            for column in COLUMNS:
                row.append(fmt(observed[sid][column]))
                if not baselines:
                    continue
                if not decoys:
                    row += ["", ""]
                    continue
                mean, sd = mean_sd(decoys[column])
                row += [fmt(mean), fmt(sd)]
            if baselines:
                row.append(str(len(decoys["dsci"])) if decoys else "")
            fh.write("\t".join(row) + "\n")

    print(f"[rfeval-metrics] wrote {len(ordered)} structure(s) to {args.output}"
          + (f" with {len(baselines)} baseline distribution(s)" if baselines else ""))


if __name__ == "__main__":
    main()
