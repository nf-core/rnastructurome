process MERGE_SOURCE_URLS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/coreutils:04a9693aaeef79aa'
        : 'community.wave.seqera.io/library/coreutils:9a0aa1447088f229'}"

    input:
    tuple val(meta), path(url_files, stageAs: "inputs/url_??.txt")

    output:
    tuple val(meta), path("ensembl_source_url.txt"), emit: urls
    tuple val("${task.process}"), val('cat'), eval("cat --version 2>&1 | head -n 1 | sed 's/^.*coreutils) //; s/ .*\$//'"), topic: versions, emit: versions_cat

    script:
    """
    cat inputs/* > ensembl_source_url.txt
    """

    stub:
    """
    touch ensembl_source_url.txt
    """
}
