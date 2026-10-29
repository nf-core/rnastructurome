process NCBI_FASTA {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), val(accessions)

    output:
    tuple val(meta), path("${meta.id}.transcripts.fa.gz"), emit: fasta
    tuple val(meta), path("ncbi_source_accessions.txt"), emit: source_accessions
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/^Python //'"), topic: versions, emit: versions_python

    script:
    // accessions may be null for unseen organisms — the script falls back to esearch
    def acc_arg = accessions ? "--accessions \"${accessions}\"" : ""
    def args = task.ext.args ?: ''
    """
    ncbi_fasta.py \
        ${acc_arg} \
        --organism "${meta.original_organism}" \
        --output "${meta.id}.transcripts.fa.gz" \
        --source-accessions "ncbi_source_accessions.txt" \
        ${args}
    """

    stub:
    """
    touch ${meta.id}.transcripts.fa.gz
    printf '%s\\n' "stub://${meta.id}" > ncbi_source_accessions.txt
    """
}
