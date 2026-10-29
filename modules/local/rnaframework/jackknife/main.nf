process RNAFRAMEWORK_RFJACKKNIFE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/rnaframework_findutils:a2ed531525463ed6' :
        'community.wave.seqera.io/library/rnaframework_findutils:70fa3c51935922d9' }"

    input:
    tuple val(meta), path(xml, stageAs: "input*/*")
    path reference

    output:
    tuple val(meta), path("${prefix}_jackknife/*.csv"), emit: csv
    tuple val(meta), path("${prefix}_jackknife/*.pdf"), emit: heatmap, optional: true
    tuple val(meta), path("${prefix}_jackknife/rfjackknife.log"), emit: log, optional: true
    tuple val("${task.process}"), val('rnaframework'), eval("rf-jackknife -h 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1 | grep . || echo unknown"), topic: versions, emit: versions_rnaframework

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    export TERM="\${TERM:-xterm}"
    log_tmp="\$(mktemp "${prefix}_jackknife.XXXXXX.log")"

    rf-jackknife \\
        -p ${task.cpus} \\
        -o ${prefix}_jackknife \\
        -ow \\
        -r ${reference} \\
        ${args} \\
        input*/ 2>&1 | tee "\${log_tmp}"

    # Strip CR-based progress animation and ANSI escape codes; keep only results
    perl -pe 's/\\r/\\n/g; s/\\e\\[[0-9;]*[A-Za-z]//g' "\${log_tmp}" \\
        | grep -v 'Jackknifing folding parameters' \\
        | cat -s \\
        > "\${log_tmp}.clean"
    mv "\${log_tmp}.clean" "\${log_tmp}"

    mkdir -p ${prefix}_jackknife
    mv "\${log_tmp}" ${prefix}_jackknife/rfjackknife.log
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_jackknife

    cat <<-END_CSV > ${prefix}_jackknife/FMI.csv
    FMI;-2;-1;0
    1.0;0.72;0.75;0.60
    1.2;0.70;0.78;0.65
    END_CSV

    touch ${prefix}_jackknife/rfjackknife.log
    """
}
