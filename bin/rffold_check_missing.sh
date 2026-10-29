#!/usr/bin/env bash
# Author: Victoria Begley (begley@ebi.ac.uk)
# License: MIT
set -euo pipefail

# Check for transcripts expected but absent from rf-fold output.
# Writes a warning log if any are missing; cleans up temp files on success.
# Usage: rffold_check_missing.sh <fold_dir> <xml_dir>

fold_dir="$1"
xml_dir="$2"

missing_list="${fold_dir}/missing_transcripts.txt"
warning_log="${fold_dir}/partial_fold_warning.log"

expected_list=$(mktemp)
folded_list=$(mktemp)

find "$xml_dir" -maxdepth 1 -name '*.xml' \
    | sed 's#.*/##; s#\.xml$##' | sort -u >| "$expected_list"
find "${fold_dir}/structures" -maxdepth 1 -name '*.db' \
    | sed 's#.*/##; s#\.db$##' | sort -u >| "$folded_list"
comm -23 "$expected_list" "$folded_list" >| "$missing_list"

expected_count=$(wc -l < "$expected_list" | tr -d ' ')
folded_count=$(wc -l < "$folded_list" | tr -d ' ')
missing_count=$(wc -l < "$missing_list" | tr -d ' ')

rm -f "$expected_list" "$folded_list"

if [[ "$missing_count" -gt 0 ]]; then
    {
        printf '[RNAFRAMEWORK_RFFOLD] Partial fold output detected.\n'
        printf '[RNAFRAMEWORK_RFFOLD] Expected %s transcript(s), folded %s, missing %s.\n' \
            "$expected_count" "$folded_count" "$missing_count"
        printf '[RNAFRAMEWORK_RFFOLD] Missing transcript IDs:\n'
        cat "$missing_list"
    } | tee "$warning_log" >&2
else
    rm -f "$warning_log" "$missing_list"
fi
