#!/usr/bin/env bash
# Author: Victoria Begley (begley@ebi.ac.uk)
# License: MIT
set -euo pipefail

# Mirror a fold output directory into a publish directory using hard links where possible.
# Usage: rffold_publish.sh <fold_dir> <publish_dir>

fold_dir="$1"
publish_dir="$2"

rm -rf "$publish_dir"
find "$fold_dir" -type f -print | while IFS= read -r file; do
    rel="${file#${fold_dir}/}"
    dest="${publish_dir}/${rel}"
    mkdir -p "$(dirname "$dest")"
    ln "$file" "$dest" 2>/dev/null || cp -p "$file" "$dest"
done
