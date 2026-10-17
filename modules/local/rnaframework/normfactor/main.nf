process RNAFRAMEWORK_RFNORMFACTOR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/rnaframework_findutils:a2ed531525463ed6' :
        'community.wave.seqera.io/library/rnaframework_findutils:70fa3c51935922d9' }"

    input:
    tuple val(meta), path(treated), path(untreated), path(denatured), path(rci_files)

    output:
    tuple val(meta), path("${prefix}.norm_factors.txt"), emit: factors, optional: true
    tuple val(meta), path("${prefix}.rfnormfactor.log"), emit: log, optional: true
    tuple val("${task.process}"), val('rnaframework'), eval("rf-normfactor -h 2>&1 | grep -oE 'v[0-9]+\\.[0-9]+\\.[0-9]+' | sed 's/v//' | head -1 | grep . || echo unknown"), topic: versions, emit: versions_rnaframework

    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args ?: ''
    prefix             = task.ext.prefix ?: "${meta.id}"
    // Staged file names are unique (a shared control can't be staged twice under the same name), so the
    // positional -t/-u/-d order (which repeats a shared control per treated sample) is carried in meta.
    def treatedNames   = (treated instanceof List ? treated : [treated]).collect { f -> f.name }
    def untreatedNames = untreated ? (untreated instanceof List ? untreated : [untreated]).collect { f -> f.name } : []
    def denaturedNames = denatured ? (denatured instanceof List ? denatured : [denatured]).collect { f -> f.name } : []
    def alignNames     = (treatedNames + untreatedNames + denaturedNames).unique()
    def treatedOrder   = meta.nf_treated_order   ?: treatedNames
    def untreatedOrder = meta.nf_untreated_order ?: untreatedNames
    def denaturedOrder = meta.nf_denatured_order ?: denaturedNames
    """
    export TERM="\${TERM:-xterm}"
    rfnormfactor_log_tmp="${prefix}.rfnormfactor.log"

    # rf-normfactor requires every input RC file to contain the SAME transcript set (it errors with
    # "Provided RC files have unequal sizes" otherwise). On the genome route, the per-sample coverage
    # pre-filter in rf-rctools extract leaves treated/untreated covering different transcripts. Align
    # the UNIQUE staged RC files to their common transcript set by re-extracting them
    # (see bin/rfnormfactor_align_rc.sh); transcript order is identical across outputs.
    align_files=(${alignNames.collect { n -> "\"${n}\"" }.join(' ')})
    rfnormfactor_align_rc.sh "${prefix}" aligned "\${align_files[@]}"

    # Build the positionally-paired -t/-u/-d lists from the meta order (repeats a shared control once
    # per treated sample). Each name points at its aligned copy; repeats reuse the same aligned file.
    aligned_list() { local out=""; for f in "\$@"; do out+="aligned/\${f},"; done; echo "\${out%,}"; }
    treated_order=(${treatedOrder.collect { n -> "\"${n}\"" }.join(' ')})
    untreated_order=(${untreatedOrder ? untreatedOrder.collect { n -> "\"${n}\"" }.join(' ') : ''})
    denatured_order=(${denaturedOrder ? denaturedOrder.collect { n -> "\"${n}\"" }.join(' ') : ''})
    treated_arg="-t \$(aligned_list "\${treated_order[@]}")"
    untreated_arg=""
    [[ \${#untreated_order[@]} -gt 0 ]] && untreated_arg="-u \$(aligned_list "\${untreated_order[@]}")"
    denatured_arg=""
    [[ \${#denatured_order[@]} -gt 0 ]] && denatured_arg="-d \$(aligned_list "\${denatured_order[@]}")"

    # Best-effort: rf-normfactor can legitimately fail to produce factors (e.g. "No bases covered
    # across all samples" when no base meets -mc across every sample). Cross-experiment normalisation
    # is auto-enabled, so a shortfall must NOT break the run — capture the status, warn, and emit no
    # factor file. The subworkflow then leaves rf-norm to normalise this reference per-sample.
    set +e
    rf-normfactor \\
        -p ${task.cpus} \\
        -o ${prefix}.norm_factors.txt \\
        -ow \\
        ${args} \\
        \${treated_arg} \\
        \${untreated_arg} \\
        \${denatured_arg} 2>&1 | tee "\${rfnormfactor_log_tmp}"
    nf_status="\${PIPESTATUS[0]}"
    set -e

    mv "\${rfnormfactor_log_tmp}" ${prefix}.rfnormfactor.log || true

    if [[ "\${nf_status}" -ne 0 || ! -s "${prefix}.norm_factors.txt" ]]; then
        echo "[RNAFRAMEWORK_RFNORMFACTOR] rf-normfactor did not produce normalisation factors for '${prefix}' (exit \${nf_status}); rf-norm will fall back to per-sample normalisation for this reference. If coverage is the cause, lower --rfnorm_normfactor_min_coverage." >&2
        rm -f "${prefix}.norm_factors.txt"
    fi

    # rf-normfactor prints its version banner only under -h and exits non-zero; capture with '|| true'
    # so 'set -o pipefail' does not abort, and fall back to 'unknown' if no version is found.
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "${args}"
    touch ${prefix}.norm_factors.txt
    touch ${prefix}.rfnormfactor.log
    """
}
