#!/usr/bin/env bash
# Author: Victoria Begley (begley@ebi.ac.uk)
# License: MIT
set -euo pipefail

# Deduplicate XML files staged from multiple replicates into unique_xml/.
# Nextflow stages each file under input<N>/ with the original filename preserved.
# Keeps the first occurrence of each transcript XML when replicates overlap.

mkdir -p unique_xml

for _f in input*/*.xml; do
    [[ -f "$_f" ]] || continue
    _base=$(basename "$_f")
    if [[ ! -e "unique_xml/$_base" ]]; then
        ln -sf "$(readlink -f "$_f")" "unique_xml/$_base"
    fi
done
