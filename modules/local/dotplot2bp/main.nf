process RNAFRAMEWORK_DOTPLOT2BP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(fold_dir), path(gtf)

    output:
    tuple val(meta), path("${meta.id}_bp/dotplot/*.bp"), emit: bp, optional: true
    tuple val(meta), path("${meta.id}_bp/conversion_warnings.log"), emit: warnings, optional: true
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/^Python //'"), topic: versions, emit: versions_python

    script:
    def args = task.ext.args ?: ''
    // Transcript-coordinate conversion (--transcript-coords) needs no GTF; genome-coordinate does.
    def gtf_arg = gtf ? "--gtf \"${gtf}\"" : ''
    """
    rnaframework_dotplot2bp.py \
        --organism "${meta.organism ?: meta.id}" \
        --prefix "${meta.id}" \
        --fold-dir "${fold_dir}" \
        ${gtf_arg} \
        ${args}
    """

    stub:
    """
    mkdir -p ${meta.id}_bp/dotplot
    touch ${meta.id}_bp/dotplot/stub.bp
    """
}
