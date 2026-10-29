process MERGE_BP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(bp_files, stageAs: "inputs/*.bp")

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.bp"), emit: bp, optional: true
    tuple val("${task.process}"), val('python'), eval("python3 --version | cut -d' ' -f2"), topic: versions, emit: versions_python

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 - "${prefix}" <<'PY'
import sys
from pathlib import Path

EXPECTED_DATA_FIELDS = 6


def is_valid_data_line(fields):
    if len(fields) != EXPECTED_DATA_FIELDS:
        return False
    for f in fields:
        if not f or f.lower() == "none":
            return False
    try:
        for i in range(1, 6):
            int(fields[i])
    except ValueError:
        return False
    return True


prefix = sys.argv[1]
bp_files = sorted(Path("inputs").glob("*.bp"))
output_path = Path(f"{prefix}.bp")

color_lines = []
data_lines = []
color_written = False
skipped = 0

for bp_file in bp_files:
    for line in bp_file.read_text().splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("color"):
            if not color_written:
                color_lines.append(stripped)
        else:
            fields = stripped.split("\\t")
            if is_valid_data_line(fields):
                data_lines.append(stripped)
            else:
                skipped += 1
                print(f"WARNING: skipping malformed/null bp line in {bp_file.name}: {stripped!r}", file=sys.stderr)
    color_written = True

if skipped:
    print(f"WARNING: {skipped} line(s) were dropped due to missing or non-numeric fields.", file=sys.stderr)

if data_lines:
    with output_path.open("w") as fh:
        for line in color_lines:
            fh.write(line + "\\n")
        for line in data_lines:
            fh.write(line + "\\n")
PY
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bp
    """
}
