process ENSEMBL_GENOME {
    tag "${meta.id}:${ensembl_species}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), val(ensembl_species)
    val  ensembl_config_input

    output:
    tuple val(meta), path("${meta.id}.genome.fa"), emit: fasta, optional: true
    tuple val(meta), path("ensembl_source_url.txt"), emit: source_url, optional: true
    tuple val(meta), path("${meta.id}.not_found"), emit: not_found, optional: true
    tuple val("${task.process}"), val('ensembl'), val("${(defaultEnsemblConfig() + (ensembl_config_input ?: [:])).ensembl_release}"), topic: versions, emit: versions_ensembl

    when:
    task.ext.when == null || task.ext.when

    script:
    def ensembl_config = defaultEnsemblConfig() + (ensembl_config_input ?: [:])
    """
    ensembl_genome.py \
        --species   "${ensembl_species}" \
        --release   "${ensembl_config.ensembl_release}" \
        --base-url  "${ensembl_config.ensembl_base_url}" \
        --output    "${meta.id}.genome.fa" \
        --source-url    "ensembl_source_url.txt" \
        --not-found-file "${meta.id}.not_found"
    """

    stub:
    """
    touch ${meta.id}.genome.fa
    printf '%s\\n' \
        "stub://${meta.id}.genome.fa" \
        > ensembl_source_url.txt
    """
}

def defaultEnsemblConfig() {
    [
        ensembl_release : 'current',
        ensembl_base_url: 'https://ftp.ensembl.org/pub'
    ]
}
