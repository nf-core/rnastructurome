process ENSEMBL_TRANSCRIPTOME {
    tag "${meta.id}:${ensembl_species}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), val(ensembl_species)
    val ensembl_config_input

    output:
    tuple val(meta), path("${meta.id}.transcripts.fa.gz"), emit: fasta, optional: true
    tuple val(meta), path("ensembl_source_url.txt"), emit: source_urls, optional: true
    tuple val(meta), path("${meta.id}.not_found"), emit: not_found, optional: true
    path "ensembl_warnings.log", emit: warnings, optional: true
    tuple val("${task.process}"), val('ensembl'), val("${(defaultEnsemblConfig() + (ensembl_config_input ?: [:])).ensembl_release}"), topic: versions, emit: versions_ensembl

    script:
    def ensembl_config = defaultEnsemblConfig() + (ensembl_config_input ?: [:])
    """
    ensembl_transcriptome.py \
        --species "${ensembl_species}" \
        --release "${ensembl_config.ensembl_release}" \
        --base-url "${ensembl_config.ensembl_base_url}" \
        --output "${meta.id}.transcripts.fa.gz" \
        --source-urls "ensembl_source_url.txt" \
        --warnings-log "ensembl_warnings.log" \
        --not-found-file "${meta.id}.not_found"
    """

    stub:
    """
    touch ${meta.id}.transcripts.fa.gz
    printf '%s\n' \
        "stub://${meta.id}.transcripts.fa.gz" \
        > ensembl_source_url.txt
    """
}

def defaultEnsemblConfig() {
    [
        ensembl_release : 'current',
        ensembl_base_url: 'https://ftp.ensembl.org/pub'
    ]
}
