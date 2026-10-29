process FASTA_SORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'
        : 'quay.io/biocontainers/seqkit:2.9.0--h9ee0642_0'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.sorted.fa"), emit: fasta
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+'"), topic: versions, emit: versions_seqkit

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    seqkit sort -w 80 "${fasta}" > "${prefix}.sorted.fa"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sorted.fa
    """
}
