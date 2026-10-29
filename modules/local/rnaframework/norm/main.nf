process RNAFRAMEWORK_RFNORM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/rnaframework_findutils:a2ed531525463ed6' :
        'community.wave.seqera.io/library/rnaframework_findutils:70fa3c51935922d9' }"

    input:
    tuple val(meta), path(treated), path(untreated), path(denatured), path(rci_files), path(norm_factor)

    output:
    tuple val(meta), path("${prefix}_norm/*.xml"), emit: xml
    tuple val(meta), path("${prefix}_norm/rfnorm.log"), emit: log, optional: true
    tuple val(meta), path("${prefix}_norm/plots/*.pdf"), emit: plots, optional: true
    tuple val("${task.process}"), val('rnaframework'), eval("rf-norm -h 2>&1 | grep -oE 'v[0-9]+\\.[0-9]+\\.[0-9]+' | sed 's/v//' | head -1 | grep . || echo unknown"), topic: versions, emit: versions_rnaframework

    script:
    def args          = task.ext.args ?: ''
    prefix            = task.ext.prefix ?: "${meta.id}"
    def untreated_arg = untreated ? "-u ${untreated}" : ''
    def denatured_arg = denatured ? "-d ${denatured}" : ''
    // Cross-experiment factor from rf-normfactor (opt-in). IMPORTANT: rf-norm's -nf takes a NUMERIC
    // value, not a path — looked up from rf-normfactor's table by this group's treated RC name.
    // Only attempted for single-treated groups; otherwise falls back to per-sample normalisation.
    def treatedNames  = (treated instanceof List ? treated : [treated]).collect { f -> f.name }
    def nfKeyName     = (norm_factor && treatedNames.size() == 1) ? treatedNames[0] : ''
    def treated_list  = treated instanceof List ? treated.join(' ') : "${treated}"
    """
    export TERM="\${TERM:-xterm}"
    rfnorm_log_tmp="${prefix}.rfnorm.log"

    # Resolve the numeric -nf value for this group from the rf-normfactor table (see note above).
    norm_factor_arg=""
    nf_key_name="${nfKeyName}"
    if [[ -n "\${nf_key_name}" && -s "${norm_factor}" ]]; then
        nf_key="\${nf_key_name%.rc}"
        if [[ "${meta.four_base_dms ? 'true' : 'false'}" == "true" ]]; then
            # rf-normfactor -ni writes one column per base; rf-norm takes "A:x,C:y", skipping NaN bases.
            nf_value="\$(awk -F'\\t' -v k="\${nf_key}" 'NR == 1 { for (i = 2; i <= NF; i++) h[i] = \$i; next } \$1 == k { for (i = 2; i <= NF; i++) if (\$i != "NaN") s = s (s == "" ? "" : ",") h[i] ":" \$i; print s; exit }' "${norm_factor}")"
        else
            nf_value="\$(awk -v k="\${nf_key}" '\$1 == k { print \$2; exit }' "${norm_factor}")"
        fi
        if [[ -n "\${nf_value}" ]]; then
            norm_factor_arg="-nf \${nf_value}"
        else
            echo "[RNAFRAMEWORK_RFNORM] no rf-normfactor entry for '\${nf_key}' in ${norm_factor}; using rf-norm internal per-sample normalisation" >&2
        fi
    fi

    rf-norm \\
        -p ${task.cpus} \\
        -o ${prefix}_norm \\
        -ow \\
        ${args} \\
        \${norm_factor_arg} \\
        -t ${treated_list} \\
        ${untreated_arg} \\
        ${denatured_arg} 2>&1 | tee "\${rfnorm_log_tmp}"

    # Strip ANSI colour codes and collapse \\r-terminated progress-bar updates into real newlines,
    # so parseRfnormLog's line-based regex can find the covered/discarded summary reliably.
    cleaned_log="${prefix}.rfnorm.clean.log"
    sed -E 's/\\x1b\\[[0-9;]*[A-Za-z]//g' "\${rfnorm_log_tmp}" | tr '\\r' '\\n' > "\${cleaned_log}"

    # Publish only the Normalization statistics block (the covered/discarded summary parseRfnormLog reads);
    # the preceding progress output is noise. Keep the full log if the marker is absent (e.g. a failure
    # before stats are printed) so errors stay debuggable.
    if grep -Fq '[+] Normalization statistics:' "\${cleaned_log}"; then
        awk 'index(\$0, "[+] Normalization statistics:") { f=1 } f' "\${cleaned_log}" > "\${cleaned_log}.trim"
        mv "\${cleaned_log}.trim" "\${cleaned_log}"
    fi

    if [[ -d ${prefix}_norm ]]; then
        mv "\${cleaned_log}" ${prefix}_norm/rfnorm.log
    fi
    rm -f "\${rfnorm_log_tmp}"

    if [[ -d "${prefix}_norm/plots" ]]; then
        find "${prefix}_norm/plots" -maxdepth 1 -type f -name '*.pdf' | while IFS= read -r plot; do
            transcript_id="\$(basename "\${plot}" .pdf)"
            if [[ ! -f "${prefix}_norm/\${transcript_id}.xml" ]]; then
                rm -f "\${plot}"
            fi
        done
    fi
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "${args}"
    mkdir -p ${prefix}_norm
    mkdir -p ${prefix}_norm/plots
    touch ${prefix}_norm/${prefix}.xml
    touch ${prefix}_norm/plots/${prefix}.pdf
    touch ${prefix}_norm/rfnorm.log
    """
}
