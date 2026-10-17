process RNAFRAMEWORK_RFCOUNT_GENOME {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/rnaframework_findutils:a2ed531525463ed6' :
        'community.wave.seqera.io/library/rnaframework_findutils:70fa3c51935922d9' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta_ref), path(fasta)

    output:
    tuple val(meta), path("*_rfcount_genome/*.rc"), emit: rc, optional: true
    tuple val(meta), path("*_rfcount_genome/*.rc.rci"), emit: rci, optional: true
    tuple val(meta), path("*_rfcount_genome/index.rci"), emit: index_rci, optional: true
    tuple val(meta), path("*_rfcount_genome/error.out"), emit: error_log, optional: true
    tuple val(meta), path("*_rfcount_genome/samtools.log"), emit: samtools_log, optional: true
    tuple val(meta), path("*_rfcount_genome/*.rfcount_genome_summary.tsv"), emit: summary, optional: true
    tuple val(meta), path("*_rfcount_genome/*.rfcount_genome.log"), emit: log, optional: true
    tuple val(meta), path("*_rfcount_genome/plots/*.pdf"), emit: plots, optional: true
    tuple val("${task.process}"), val('rnaframework'), eval("rf-count-genome -h 2>&1 | grep -oE 'v[0-9]+\\.[0-9]+\\.[0-9]+' | sed 's/v//' | head -1 | grep . || echo unknown"), topic: versions, emit: versions_rnaframework

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outdir = "${prefix}_rfcount_genome"
    // Mirrors the -m (mutation mode) decision in conf/modules.config, so the summary parser knows
    // which rf-count-genome table layout to expect (MaP has a Mutated-alignments column; RT-stop does not).
    def is_map = ((meta.principle ?: '').toLowerCase() == 'map') ? '1' : '0'
    """
    FASTA_PATH="${fasta}"

    export TERM="\${TERM:-xterm}"

    mkdir -p ${outdir}
    rfcount_outdir="${outdir}"
    rfcount_log_tmp="${prefix}.rfcount_genome.log"

    set -o pipefail
    set +e
    rf-count-genome \\
        -p ${task.cpus} \\
        -f "\${FASTA_PATH}" \\
        -o ${outdir} \\
        -ow \\
        ${args} \\
        "${bam}" 2>&1 | tee "\${rfcount_log_tmp}"
    pipeline_statuses=( "\${PIPESTATUS[@]}" )
    rfcount_status="\${pipeline_statuses[0]}"
    tee_status="\${pipeline_statuses[1]}"
    set -e

    if [[ "\${tee_status}" -ne 0 ]]; then
        echo "[RNAFRAMEWORK_RFCOUNT_GENOME] tee failed while writing rf-count-genome log for sample '${prefix}'." >&2
        exit "\${tee_status}"
    fi

    cleaned_log="${prefix}.rfcount_genome.clean.log"
    sed -E 's/\\x1b\\[[0-9;]*[A-Za-z]//g' "\${rfcount_log_tmp}" | tr '\\r' '\\n' > "\${cleaned_log}"

    rfcount_completed_with_nonzero=0
    if [[ "\${rfcount_status}" -ne 0 ]]; then
        if grep -Fq '[+] All done.' "\${cleaned_log}"; then
            rfcount_completed_with_nonzero=1
        else
            echo "[RNAFRAMEWORK_RFCOUNT_GENOME] rf-count-genome failed with exit status \${rfcount_status} for sample '${prefix}'." >&2
            exit "\${rfcount_status}"
        fi
    fi

    summary_tsv="${outdir}/${prefix}.rfcount_genome_summary.tsv"
    # match_mode=prefix: the summary's sample column carries the staged BAM's filename suffix.
    rfcount_parse_summary.awk \\
        -v sample="${prefix}" -v is_map="${is_map}" -v match_mode="prefix" \\
        "\${cleaned_log}" > "\${summary_tsv}"

    # rf-count-genome reports "Covered: 0" when run without a -a annotation file
    # (genome-wide mode); that is expected here — rf-rctools extract handles transcript
    # extraction using the GTF in the next step.  Fail only if no RC files were produced.
    rc_count=\$(find "\${rfcount_outdir}" -type f -name '*.rc' 2>/dev/null | wc -l)
    if [[ "\${rc_count}" -eq 0 ]]; then
        echo "[RNAFRAMEWORK_RFCOUNT_GENOME] rf-count-genome produced no RC files for sample '${prefix}'." >&2
        exit 1
    fi
    if [[ "\${rfcount_completed_with_nonzero}" -eq 1 ]]; then
        echo "[RNAFRAMEWORK_RFCOUNT_GENOME] rf-count-genome exited with status \${rfcount_status} after reporting completion; continuing because RC files were produced." >&2
    fi

    # Publish only the final statistics section (from "[+] Statistics:" to the end); the verbose
    # per-transcript progress above it is noise. Fall back to the full log if the marker is absent,
    # e.g. an early failure, so nothing useful is lost.
    if grep -Fq '[+] Statistics:' "\${cleaned_log}"; then
        awk '/^\\[\\+\\] Statistics:/{p=1} p' "\${cleaned_log}" > "${outdir}/${prefix}.rfcount_genome.log"
    else
        cp "\${cleaned_log}" "${outdir}/${prefix}.rfcount_genome.log"
    fi
    rm -f "\${cleaned_log}" "\${rfcount_log_tmp}"
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outdir = "${prefix}_rfcount_genome"
    """
    echo "${args}"
    mkdir -p ${outdir}
    touch ${outdir}/${prefix}.plus.rc
    touch ${outdir}/${prefix}.plus.rc.rci
    touch ${outdir}/index.rci
    touch ${outdir}/error.out
    touch ${outdir}/samtools.log
    cat <<-END_SUMMARY > ${outdir}/${prefix}.rfcount_genome_summary.tsv
    sample	covered	pct_mutated	pct_a_muts	pct_c_muts	pct_g_muts	pct_u_muts
    ${prefix}	1	25.0	25.0	25.0	25.0
    END_SUMMARY
    """
}
