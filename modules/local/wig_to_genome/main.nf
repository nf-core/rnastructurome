process WIG_TO_GENOME {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(wig), path(gtf)

    output:
    tuple val(meta), path("${prefix}.genomic.wig"), emit: wig
    tuple val(meta), path("${prefix}_genomic.chrom.sizes"), emit: chrom_sizes
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/^Python //'"), topic: versions, emit: versions_python

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    remap_wig_to_genome.py \
        --wig "${wig}" \
        --gtf "${gtf}" \
        --output-wig "${prefix}.genomic.wig" \
        --chrom-sizes "${prefix}_genomic.chrom.sizes" \
        --organism "${meta.organism ?: meta.id}" \
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.genomic.wig
    touch ${prefix}_genomic.chrom.sizes
    """
}
