process ENSEMBL_GTF {
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
    tuple val(meta), path("${meta.id}.annotation.gtf"), emit: gtf, optional: true
    tuple val(meta), path("${meta.id}.annotation.gtf.gz"), emit: gtf_gz, optional: true
    tuple val(meta), path("ensembl_source_url.txt"), emit: source_urls, optional: true
    tuple val(meta), path("${meta.id}.not_found"), emit: not_found, optional: true
    tuple val("${task.process}"), val('ensembl'), val("${(defaultEnsemblConfig() + (ensembl_config_input ?: [:])).ensembl_release}"), topic: versions, emit: versions_ensembl

    script:
    def ensembl_config = defaultEnsemblConfig() + (ensembl_config_input ?: [:])
    """
    ensembl_gtf.py \
        --species "${ensembl_species}" \
        --release "${ensembl_config.ensembl_release}" \
        --base-url "${ensembl_config.ensembl_base_url}" \
        --output "${meta.id}.annotation.gtf.gz" \
        --source-urls "ensembl_source_url.txt" \
        --not-found-file "${meta.id}.not_found"

    [[ -f "${meta.id}.annotation.gtf.gz" ]] && gzip -dc "${meta.id}.annotation.gtf.gz" > "${meta.id}.annotation.gtf" || true
    """

    stub:
    """
    touch ${meta.id}.annotation.gtf.gz
    touch ${meta.id}.annotation.gtf
    printf '%s\n' \
        "stub://${meta.id}.annotation.gtf.gz" \
        > ensembl_source_url.txt
    """
}

def defaultEnsemblConfig() {
    [
        ensembl_release : 'current',
        ensembl_base_url: 'https://ftp.ensembl.org/pub'
    ]
}
