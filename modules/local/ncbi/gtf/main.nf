process NCBI_GTF {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.annotation.gtf"), emit: gtf
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/^Python //'"), topic: versions, emit: versions_python

    script:
    """
    python3 - "${fasta}" "${meta.id}.annotation.gtf" <<'PY'
import gzip
import sys


def iter_fasta_lengths(fasta_path):
    opener = gzip.open if fasta_path.endswith(".gz") else open
    accession = None
    length = 0
    with opener(fasta_path, "rt", errors="ignore") as fh:
        for line in fh:
            line = line.rstrip("\\n")
            if line.startswith(">"):
                if accession is not None:
                    yield accession, length
                accession = line[1:].split()[0]
                length = 0
            else:
                length += len(line.strip())
    if accession is not None:
        yield accession, length


fasta_path, output_path = sys.argv[1], sys.argv[2]
records = list(iter_fasta_lengths(fasta_path))

if not records:
    print(f"[NCBI_GTF] Error: no sequences found in {fasta_path!r}", file=sys.stderr)
    sys.exit(1)

with open(output_path, "w", encoding="utf-8") as out:
    for accession, length in records:
        attrs = f'gene_id "{accession}"; transcript_id "{accession}";'
        out.write(f"{accession}\\tncbi\\ttranscript\\t1\\t{length}\\t.\\t+\\t.\\t{attrs}\\n")
        out.write(f"{accession}\\tncbi\\texon\\t1\\t{length}\\t.\\t+\\t.\\t{attrs}\\n")

print(f"[NCBI_GTF] Written {len(records)} transcript record(s) to {output_path!r}.", file=sys.stderr)
PY
    """

    stub:
    """
    touch ${meta.id}.annotation.gtf
    """
}
