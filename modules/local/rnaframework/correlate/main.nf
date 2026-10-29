process RNAFRAMEWORK_RFCORRELATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/rnaframework_findutils:a2ed531525463ed6' :
        'community.wave.seqera.io/library/rnaframework_findutils:70fa3c51935922d9' }"

    input:
    tuple val(meta), path(xmls, stageAs: "input*/*")

    output:
    tuple val(meta), path("${prefix}_correlate/matrix.csv"), emit: matrix, optional: true
    tuple val(meta), path("${prefix}_correlate/"), emit: results, optional: true
    tuple val(meta), path("${prefix}.rfcorrelate.log"), emit: log, optional: true
    tuple val("${task.process}"), val('rnaframework'), eval("rf-correlate -h 2>&1 | grep -oE 'v[0-9]+\\.[0-9]+\\.[0-9]+' | sed 's/v//' | head -1 | grep . || echo unknown"), topic: versions, emit: versions_rnaframework

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    def sizes  = meta.correlate_replicate_sizes.join(' ')
    def labels = meta.correlate_replicate_labels.join(' ')
    def total  = meta.correlate_replicate_sizes.sum()
    """
    export TERM="\${TERM:-xterm}"

    abs_path() { perl -MCwd=abs_path -e 'print abs_path(shift)' "\$1"; }

    # Rebuild one folder per replicate from the per-file input*/ staging. correlate_replicate_sizes
    # gives how many input dirs belong to each replicate (in order), correlate_replicate_labels the
    # rf-correlate sample label for each. Transcript XML names repeat across replicates, so each
    # replicate must live in its own folder (rf-correlate is called as label:folder/).
    sizes=(${sizes})
    labels=(${labels})
    _input_dirs=()
    while IFS= read -r d; do _input_dirs+=( "\${d}" ); done < <(printf '%s\\n' input*/ | sort -t t -k2,2n)  # printf is a builtin: no ARG_MAX limit at high replicate/transcript counts
    if [[ \${#_input_dirs[@]} -ne ${total} ]]; then
        echo "[RNAFRAMEWORK_RFCORRELATE] staged input dir count (\${#_input_dirs[@]}) != expected (${total})." >&2
        exit 1
    fi
    corr_args=()
    idx=0
    for k in "\${!sizes[@]}"; do
        repdir="corr_input/\${labels[k]}"
        mkdir -p "\${repdir}"
        for ((j = 0; j < sizes[k]; j++)); do
            for x in "\${_input_dirs[idx]}"*.xml; do
                [[ -e "\${x}" ]] && ln -sf "\$(abs_path "\${x}")" "\${repdir}/\$(basename "\${x}")"
            done
            idx=\$((idx + 1))
        done
        corr_args+=( "\${labels[k]}:\${repdir}/" )
    done

    rf-correlate \\
        -p ${task.cpus} \\
        -o ${prefix}_correlate \\
        -ow \\
        ${args} \\
        "\${corr_args[@]}" 2>&1 | tee "${prefix}.rfcorrelate.log"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_correlate/pairwise
    printf 'Sample,repA,repB\\nrepA,1,0.9\\nrepB,0.9,1\\n' > ${prefix}_correlate/matrix.csv
    touch ${prefix}.rfcorrelate.log
    """
}
