process FASTA_SANITIZE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.sanitized.fa"), emit: fasta
    tuple val("${task.process}"), val('python'), eval("python3 --version | cut -d' ' -f2"), topic: versions, emit: versions_python

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 - "${fasta}" "${prefix}.sanitized.fa" "${meta.organism ?: meta.id}" <<'PY'
import gzip
import re
import sys
from pathlib import Path

YEAST_ISOFORM_PATTERN = re.compile(
    r'^(Y[A-P][LR][0-9]{3}[CW])-([A-Z])(_(?:mRNA|ncRNA|snRNA|snoRNA|rRNA|tRNA))\$'
)
UNSAFE_ID_CHARS = re.compile(r'[^A-Za-z0-9._-]')


def open_fasta(path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def sanitize_header(header, organism):
    token, sep, remainder = header.partition(" ")
    if organism == "saccharomyces_cerevisiae":
        token = YEAST_ISOFORM_PATTERN.sub(r'\\1_\\2\\3', token)
    token = UNSAFE_ID_CHARS.sub("_", token)
    return f"{token}{sep}{remainder}" if sep else token


input_path = Path(sys.argv[1])
output_path = Path(sys.argv[2])
organism = sys.argv[3].strip().lower().replace(" ", "_")

with open_fasta(input_path) as handle, output_path.open("wt", encoding="utf-8") as out_handle:
    for line in handle:
        if line.startswith(">"):
            header = sanitize_header(line[1:].rstrip("\\n"), organism)
            out_handle.write(f">{header}\\n")
        else:
            out_handle.write(line)
PY
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sanitized.fa
    """
}
