process GTF_SANITIZE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("${meta.id}.sanitized.gtf"), emit: gtf
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | sed 's/^Python //'"), topic: versions, emit: versions_python

    script:
    """
    python3 - "${gtf}" "${meta.id}.sanitized.gtf" <<'PY'
import gzip
import re
import sys
from pathlib import Path

UNSAFE = re.compile(r"[^A-Za-z0-9._-]")
ATTR = re.compile(r'((?:transcript_id|gene_id)\\s+")([^"]*)(")')


def sanitize_id(value):
    return UNSAFE.sub("_", value)


def _open(path, mode):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8")
    return path.open(mode, encoding="utf-8")


def sanitize_line(line):
    return ATTR.sub(lambda m: f"{m.group(1)}{sanitize_id(m.group(2))}{m.group(3)}", line)


input_path = Path(sys.argv[1])
output_path = Path(sys.argv[2])
changed = 0
with _open(input_path, "rt") as fin, output_path.open("wt", encoding="utf-8") as fout:
    for line in fin:
        if line.startswith("#"):
            fout.write(line)
            continue
        new_line = sanitize_line(line)
        if new_line != line:
            changed += 1
        fout.write(new_line)

print(f"[sanitize_gtf_ids] Rewrote IDs on {changed} line(s) of {input_path.name}.", file=sys.stderr)
PY
    """

    stub:
    """
    touch ${meta.id}.sanitized.gtf
    """
}
