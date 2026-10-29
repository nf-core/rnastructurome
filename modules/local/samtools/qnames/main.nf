process SAMTOOLS_QNAMES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.qnames.txt"), emit: qnames
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed -n '1s/samtools //p'"), topic: versions, emit: versions_samtools

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Bound sort's RAM below the request and spill to the work dir (-T .), not the node's small $TMPDIR.
    def sort_mem = task.memory ? Math.max(1, task.memory.toGiga().intValue() - 2) : 4
    """
    samtools view -@ ${task.cpus} ${bam} | cut -f1 | LC_ALL=C sort -u -T . -S ${sort_mem}G > ${prefix}.qnames.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.qnames.txt
    """
}
