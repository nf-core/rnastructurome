process RNAFRAMEWORK_RFRCTOOLS_EXTRACT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/rnaframework_findutils:a2ed531525463ed6' :
        'community.wave.seqera.io/library/rnaframework_findutils:70fa3c51935922d9' }"

    input:
    tuple val(meta), path(rc, stageAs: "input/*"), path(rci, stageAs: "input/*"), path(summary)
    tuple val(meta_ref), path(gtf)

    output:
    tuple val(meta), path("${prefix}_rctools_extract/${prefix}.rc"), emit: rc
    tuple val(meta), path("${prefix}_rctools_extract/${prefix}.rc.rci"), emit: rci, optional: true
    tuple val(meta), path("${prefix}_rctools_extract/${prefix}.rfcount_genome_summary.tsv"), emit: summary, optional: true
    tuple val("${task.process}"), val('rnaframework'), eval("rf-rctools -h 2>&1 | sed -nE 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/p' | head -1 | grep . || echo unknown"), topic: versions, emit: versions_rnaframework

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def feature_name  = task.ext.gtf_feature ?: 'exon'
    def attr_name     = task.ext.gtf_attribute ?: 'transcript_id'
    def min_coverage  = task.ext.min_coverage != null ? task.ext.min_coverage : 1
    prefix        = task.ext.prefix ?: "${meta.id}"
    def outdir    = "${prefix}_rctools_extract"
    """
    export TERM="\${TERM:-xterm}"
    mkdir -p ${outdir}

    # Build per-file RCI indexes so rf-rctools can do strand-aware extraction.
    # set -e here so index failures surface rather than silently producing a
    # successful cached task with no .rci output.
    set -e
    for f in input/*.rc; do
        rf-rctools index "\${f}"
    done
    set +e

    # rf-rctools extract requires the BASENAME (no extension) to activate strand-aware
    # extraction: it auto-discovers Sample.plus.rc + Sample.minus.rc and uses the GTF
    # strand field to select the correct file per transcript.  Passing individual *.rc
    # paths disables this logic and produces all-zero counts.
    _rc_first=\$(ls input/*.rc 2>/dev/null | head -1)
    _rc_filename=\$(basename "\${_rc_first}")
    _rc_basename="\${_rc_filename%.plus.rc}"
    _rc_basename="\${_rc_basename%.minus.rc}"
    rc_basename="\${_rc_basename%.rc}"

    rf-rctools extract \\
        -a ${gtf} \\
        -o ${outdir}/${prefix}.rc \\
        -ow \\
        ${args} \\
        input/\${rc_basename}

    set -e
    rf-rctools index ${outdir}/${prefix}.rc
    set +e

    # Coverage filter: the genome route extracts EVERY GTF transcript, the vast majority with zero
    # coverage. Filter the RC down to transcripts with real coverage so the emitted (and published)
    # RC is transcriptome-sized — rf-norm can then run on it whole, using its own -p threading,
    # instead of needing to be split into chunks. The view 'coverage' track (4th line per transcript)
    # is the ground truth (rf-rctools stats does not report usable per-transcript coverage).
    # min_coverage=0 disables the filter and keeps the full annotation RC.
    if [[ ${min_coverage} -gt 0 ]]; then
        set -e
        rf-rctools view ${outdir}/${prefix}.rc 2>/dev/null | awk -v mc=${min_coverage} '
            NF == 0 { line = 0; next }
            { line++ }
            line == 1 { id = \$0; next }
            line == 4 {
                n = split(\$0, cov, ",")
                for (i = 1; i <= n; i++) if (cov[i] + 0 >= mc) { print id; break }
                line = 0
            }
        ' > covered_ids.txt

        if [[ ! -s covered_ids.txt ]]; then
            echo "ERROR: no covered transcripts in ${outdir}/${prefix}.rc (min_coverage=${min_coverage})." >&2
            exit 1
        fi

        # Real per-transcript lengths from the GTF (spliced length = sum of exon lengths). These match
        # the RC exactly because it was built from this same GTF. 4-column BED (id 0 length id) keeps
        # the clean transcript ID instead of renaming the region to <id>_0-<end>.
        rctools_covered_bed.py "${gtf}" "${feature_name}" "${attr_name}"

        if [[ ! -s covered.bed ]]; then
            echo "ERROR: covered transcripts did not match any GTF ${attr_name} for ${prefix}." >&2
            exit 1
        fi

        rf-rctools extract -a covered.bed -o ${outdir}/${prefix}.filtered.rc -ow ${outdir}/${prefix}.rc
        rf-rctools index ${outdir}/${prefix}.filtered.rc
        mv ${outdir}/${prefix}.filtered.rc     ${outdir}/${prefix}.rc
        mv ${outdir}/${prefix}.filtered.rc.rci ${outdir}/${prefix}.rc.rci
        set +e
    fi

    # Rewrite the rf-count-genome summary's 'covered' column (a genome reference/contig count) with the
    # number of transcripts in the emitted RC, so the published metric reflects covered transcripts
    # rather than genome locations. Genome route only — this module is not used for --transcriptome,
    # where rf-count's summary already reports covered transcripts.
    if [[ -s covered_ids.txt ]]; then
        _covered=\$(wc -l < covered_ids.txt | tr -d ' ')
    else
        _covered=\$(rf-rctools view ${outdir}/${prefix}.rc 2>/dev/null | awk 'NF == 0 { line = 0; next } { line++ } line == 1 { c++ } line == 4 { line = 0 } END { print c + 0 }')
    fi
    awk -v cov="\${_covered}" -v sample="${prefix}" 'BEGIN { FS = OFS = "\\t" }
        NR == 1 { print; next }
        { \$2 = cov; print; seen = 1 }
        END { if (!seen) print sample, cov, "", "", "", "", "", "" }' \\
        ${summary} > ${outdir}/${prefix}.rfcount_genome_summary.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def outdir = "${prefix}_rctools_extract"
    """
    mkdir -p ${outdir}
    touch ${outdir}/${prefix}.rc
    touch ${outdir}/${prefix}.rc.rci
    cp ${summary} ${outdir}/${prefix}.rfcount_genome_summary.tsv 2>/dev/null || touch ${outdir}/${prefix}.rfcount_genome_summary.tsv
    """
}
