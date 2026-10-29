process RNAFRAMEWORK_RFEVAL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/rnaframework_findutils:a2ed531525463ed6' :
        'community.wave.seqera.io/library/rnaframework_findutils:70fa3c51935922d9' }"

    input:
    tuple val(meta), path(xml, stageAs: "xml_input*/*")
    path structures
    path windows

    output:
    tuple val(meta), path("${prefix}_rfeval/${prefix}_rfeval.metrics.tsv"), emit: metrics
    tuple val(meta), path("${prefix}_rfeval/plots/**.pdf"), emit: plots, optional: true
    tuple val(meta), path("${prefix}_rfeval/rfeval.log"), emit: log, optional: true
    tuple val(meta), path("${prefix}_rfeval_windows/*.xml"), emit: windows, optional: true
    tuple val("${task.process}"), val('rnaframework'), eval("rf-eval -h 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1 | grep . || echo unknown"), topic: versions, emit: versions_rnaframework
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/^Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    // Plot flags, main call only: the baseline scores hundreds of decoys per structure and
    // plotting them would dominate runtime. Kept separate rather than stripped back out of
    // args, because -R's value can be a $(command -v R) substitution that no regex survives.
    def args3 = task.ext.args3 ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}"
    // With a windows manifest, slice reactivities to each reference region first so a
    // sub-region structure is scored against a matching XML, not the full transcript.
    def window_cmd = windows
        ? "rnaframework_rfeval_window.py --windows \"${windows}\" --xml-glob 'rfeval_xml/*.xml' --outdir ${prefix}_rfeval_windows"
        : ''
    def reactivity_dir = windows ? "${prefix}_rfeval_windows" : 'rfeval_xml'
    def baseline_cmd = args2 ? """
    printf '\\n===== rf-eval: rotation baseline =====\\n' >> "\${log_tmp}"
    printf 'Decoys below are rotations of each profile, not real structures.\\n\\n' >> "\${log_tmp}"

    rnaframework_rfeval_baseline.py \\
        --xml-glob '${reactivity_dir}/*.xml' \\
        --db ${structures} \\
        --outdir ${prefix}_rfeval_baseline/decoys \\
        --decoy-db ${prefix}_rfeval_baseline/decoy.db \\
        ${args2}

    rf-eval \\
        -p ${task.cpus} \\
        -o ${prefix}_rfeval_baseline/scores \\
        -ow \\
        -no \\
        -s ${prefix}_rfeval_baseline/decoy.db \\
        -r ${prefix}_rfeval_baseline/decoys \\
        ${args} >> "\${log_tmp}" 2>&1
""" : ''
    def baseline_metrics_arg = args2 ? "--baseline-metrics ${prefix}_rfeval_baseline/scores/metrics.txt" : ''
    """
    export TERM="\${TERM:-xterm}"
    log_tmp="\$(mktemp "${prefix}_rfeval.XXXXXX.log")"

    # stageAs indexes each XML into its own xml_input<N>/ dir; rf-eval takes a single
    # -r folder, so gather them into one before slicing or scoring.
    mkdir -p rfeval_xml
    for f in xml_input*/*.xml; do
        if [ -e "\$f" ]; then ln -sf "\$PWD/\$f" rfeval_xml/; fi
    done

    ${window_cmd}

    printf '===== rf-eval: reference structures =====\\n' >> "\${log_tmp}"

    # -no is required: the pooled "Overall" stats divide by zero when every reference
    # window is NaN, and rnaframework_rfeval_metrics.py discards that row regardless.
    rf-eval \\
        -p ${task.cpus} \\
        -o ${prefix}_rfeval \\
        -ow \\
        -no \\
        -s ${structures} \\
        -r ${reactivity_dir} \\
        ${args} ${args3} 2>&1 | tee -a "\${log_tmp}"
${baseline_cmd}
    rnaframework_rfeval_metrics.py \\
        --metrics ${prefix}_rfeval/metrics.txt \\
        ${baseline_metrics_arg} \\
        --output ${prefix}_rfeval/${prefix}_rfeval.metrics.tsv

    # Strip ANSI/CR progress noise; drop [+] status and | progress-bar lines
    perl -pe 's/\\r/\\n/g; s/\\e\\[[0-9;]*[A-Za-z]//g' "\${log_tmp}" \\
        | grep -vE '^\\[\\+\\]|^\\|' \\
        | cat -s \\
        > "\${log_tmp}.clean"
    mv "\${log_tmp}.clean" "\${log_tmp}"

    mkdir -p ${prefix}_rfeval
    mv "\${log_tmp}" ${prefix}_rfeval/rfeval.log
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def baseline_on = task.ext.args2 as Boolean
    def stub_windows = windows
        ? "mkdir -p ${prefix}_rfeval_windows && touch ${prefix}_rfeval_windows/stub_element.xml"
        : ''
    def header = baseline_on
        ? 'structure\\tcoeff_unpaired\\tcoeff_unpaired_baseline_mean\\tcoeff_unpaired_baseline_sd\\tdsci\\tdsci_baseline_mean\\tdsci_baseline_sd\\tauroc\\tauroc_baseline_mean\\tauroc_baseline_sd\\tbaseline_num'
        : 'structure\\tcoeff_unpaired\\tdsci\\tauroc'
    def row = baseline_on
        ? 'stub_element\\t0.8200\\t0.7500\\t0.0500\\t0.7900\\t0.3800\\t0.0600\\t0.8500\\t0.5000\\t0.0600\\t40'
        : 'stub_element\\t0.8200\\t0.7900\\t0.8500'
    """
    mkdir -p ${prefix}_rfeval/plots
    ${stub_windows}

    printf '${header}\\n${row}\\n' > ${prefix}_rfeval/${prefix}_rfeval.metrics.tsv

    touch ${prefix}_rfeval/rfeval.log
    """
}
