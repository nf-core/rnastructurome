/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC             } from '../modules/nf-core/multiqc/main'
include { RNAFRAMEWORK_TORDAT } from '../modules/local/tordat/main'

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_rnastructurome_pipeline'
include { FASTQ_QC_TRIM            } from '../subworkflows/local/fastq_qc_trim/main'
include { PREPARE_REFERENCES       } from '../subworkflows/local/prepare_references/main'
include { ALIGN_READS              } from '../subworkflows/local/align_reads/main'
include { QUANTIFY_REACTIVITY      } from '../subworkflows/local/quantify_reactivity/main'
include { NORMALISE_REACTIVITIES   } from '../subworkflows/local/normalise_reactivities/main'
include { FOLD_STRUCTURES          } from '../subworkflows/local/fold_structures/main'
include { CORRELATE_REPLICATES     } from '../subworkflows/local/correlate_replicates/main'
include { VISUALISE_STRUCTURES     } from '../subworkflows/local/visualise/main'
include { BROWSER_TRACKS           } from '../subworkflows/local/browser_tracks/main'

// Pure helper functions (parsers, arg renderers, MultiQC table builders)
include {
    resolveReferenceKey
    parseFlagstatMappedReads
    parseRfcountCoveredTranscripts
    parseRfnormLog
    parseRffoldLog
    parseRfevalMetrics
    countProgressionMultiqc
    rfnormStatsMultiqc
    rffoldStatsMultiqc
    rfevalStatsMultiqc
    filterSummaryParams
    addModuleOptionsSummary
} from '../subworkflows/local/utils_nfcore_rnastructurome_pipeline/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASTRUCTUROME {

    take:
    ch_samplesheet         // channel: samplesheet read in from --input
    transcriptome          // boolean: transcriptome (Bowtie) route, including auto-detection
    main:

    ch_multiqc_files = channel.empty()
    // SUBWORKFLOW: FASTQ_QC_TRIM — cat → FastQC → UMI → cutadapt → FastQC
    FASTQ_QC_TRIM (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_QC_TRIM.out.multiqc_files)

    def ch_samplesheet_for_branching = FASTQ_QC_TRIM.out.reads_branched
    def ch_rtstop_trimmed_for_align = FASTQ_QC_TRIM.out.rtstop_trimmed
    def ch_map_trimmed_for_align    = FASTQ_QC_TRIM.out.map_trimmed

    // SUBWORKFLOW: PREPARE_REFERENCES — resolve/fetch/sort/index reference artifacts
    PREPARE_REFERENCES (
        ch_samplesheet_for_branching,
        transcriptome
    )

    def ch_reference_fasta_map          = PREPARE_REFERENCES.out.fasta_map
    def ch_reference_gtf_map            = PREPARE_REFERENCES.out.gtf_map
    def ch_all_reference_gtf            = PREPARE_REFERENCES.out.all_gtf
    def ch_reference_genome_fasta_keyed = PREPARE_REFERENCES.out.genome_fasta_keyed
    def ch_genome_transcript_fasta_map     = PREPARE_REFERENCES.out.genome_transcript_fasta_map
    def ch_genome_transcript_fasta_fai_map = PREPARE_REFERENCES.out.genome_transcript_fasta_fai_map
    def ch_star_index                   = PREPARE_REFERENCES.out.star_index
    def ch_bowtie_index_map             = PREPARE_REFERENCES.out.bowtie_index_map
    def ch_bowtie2_index_map            = PREPARE_REFERENCES.out.bowtie2_index_map

    // SUBWORKFLOW: ALIGN_READS — align → sort → dedup → stats → strandedness
    ALIGN_READS (
        ch_rtstop_trimmed_for_align,
        ch_map_trimmed_for_align,
        ch_star_index,
        ch_all_reference_gtf,
        ch_bowtie_index_map,
        ch_bowtie2_index_map,
        ch_reference_fasta_map,
        ch_genome_transcript_fasta_fai_map,
        transcriptome
    )
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN_READS.out.multiqc_files)

    def ch_markdup_bam_bai    = ALIGN_READS.out.markdup_bam_bai
    def ch_dedup_bam          = ALIGN_READS.out.dedup_bam
    def ch_strandedness_by_id = ALIGN_READS.out.strandedness_by_id
    def ch_transcript_bam_bai = ALIGN_READS.out.transcript_bam_bai

    // SUBWORKFLOW: QUANTIFY_REACTIVITY — rf-count(-genome) → rf-rctools → RC files
    QUANTIFY_REACTIVITY (
        ch_markdup_bam_bai,
        ch_transcript_bam_bai,
        ch_strandedness_by_id,
        ch_reference_genome_fasta_keyed,
        ch_reference_gtf_map,
        ch_reference_fasta_map,
        ch_genome_transcript_fasta_map,
        transcriptome
    )

    def ch_rfcount_rc      = QUANTIFY_REACTIVITY.out.rc
    def ch_rfcount_rci     = QUANTIFY_REACTIVITY.out.rci
    def ch_rfcount_summary = QUANTIFY_REACTIVITY.out.summary
    def ch_rfcount_plots   = QUANTIFY_REACTIVITY.out.plots

    // Aggregate per-sample rf-count summaries into a single TSV for the whole run
    ch_rfcount_summary
        .map { _meta, tsv -> tsv }
        .collectFile(
            name: 'rfcount_summary_all_samples.tsv',
            keepHeader: true,
            skip: 1,
            storeDir: "${params.outdir}/count"
        )

    def ch_pre_dedup_mapped_reads = ALIGN_READS.out.flagstat_pre
        .map { meta, flagstat -> [ meta.id.toString(), parseFlagstatMappedReads(flagstat) ] }

    def ch_post_dedup_mapped_reads = ALIGN_READS.out.flagstat_post
        .map { meta, flagstat -> [ meta.id.toString(), parseFlagstatMappedReads(flagstat) ] }

    def ch_rfcount_covered_transcripts = ch_rfcount_summary
        .map { meta, summary_tsv -> [ meta.id.toString(), parseRfcountCoveredTranscripts(summary_tsv) ] }

    // Anchor on post-dedup + covered (present for every counted sample); pre-dedup counts exist only for
    // samples that were deduplicated (others skip the redundant flagstat), so fall back pre = post.
    def ch_count_progression_mqc = ch_post_dedup_mapped_reads
        .join(ch_rfcount_covered_transcripts)
        .join(ch_pre_dedup_mapped_reads, remainder: true)
        .filter { entry -> entry[1] != null }
        .map { sample_id, mapped_post, covered_transcripts, mapped_pre ->
            def mappedPostLong = mapped_post as long
            def mappedPreLong = (mapped_pre != null ? mapped_pre : mapped_post) as long
            def removed = Math.max(mappedPreLong - mappedPostLong, 0L)
            def pctRemoved = mappedPreLong ? ((removed as double) / (mappedPreLong as double)) * 100.0d : 0.0d
            [
                sample_id,
                [
                    mapped_reads_pre_dedup     : mappedPreLong,
                    mapped_reads_post_dedup    : mappedPostLong,
                    pct_removed_by_dedup       : pctRemoved,
                    rfcount_covered_transcripts: covered_transcripts as long
                ]
            ]
        }
        .collect()
        .map { rows -> countProgressionMultiqc(rows) }

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_count_progression_mqc.collectFile(
            name: 'count_progression_mqc.yaml',
            sort: true
        )
    )

    // SUBWORKFLOW: NORMALISE_REACTIVITIES — group RC files, pair treated/untreated, run rf-norm
    NORMALISE_REACTIVITIES (
        ch_rfcount_rc,
        ch_rfcount_rci
    )

    // SUBWORKFLOW: FOLD_STRUCTURES — group by sample_group, optional jackknife, rf-fold, dotplot→bp
    FOLD_STRUCTURES (
        NORMALISE_REACTIVITIES.out.xml,
        ch_reference_gtf_map
    )

    if (!params.stop_after_jackknife) {
        // SUBWORKFLOW: VISUALISE_STRUCTURES — R2DT / ViennaRNA 2D structure diagrams
        VISUALISE_STRUCTURES (
            NORMALISE_REACTIVITIES.out.xml,
            FOLD_STRUCTURES.out.structures,
            FOLD_STRUCTURES.out.fold_input,
            ch_reference_fasta_map,
            ch_reference_gtf_map
        )

        // SUBWORKFLOW: BROWSER_TRACKS — rf-wiggle → BigWig genome/transcript tracks (+ Shannon)
        BROWSER_TRACKS (
            NORMALISE_REACTIVITIES.out.xml,
            FOLD_STRUCTURES.out.shannon_wig,
            ch_reference_gtf_map
        )
    }

    if (!params.stop_after_jackknife) {
        // MODULE: tordat — compile rf-norm XML + rf-fold .db structures into RDAT format. Build per-reference
        // source file name channels so the RDAT COMMENT records exact Ensembl filenames/NCBI accessions.
        def ch_ref_fasta_names
        def ch_ref_gtf_names
        if (params.fasta) {
            def local_fasta_name = file(params.fasta.toString()).name
            def local_gtf_name   = params.gtf ? file(params.gtf.toString()).name : ''
            ch_ref_fasta_names = PREPARE_REFERENCES.out.local_fasta_sorted
                .map { meta, _fasta -> [ meta.id.toString(), local_fasta_name ] }
            ch_ref_gtf_names = PREPARE_REFERENCES.out.local_fasta_sorted
                .map { meta, _fasta -> [ meta.id.toString(), local_gtf_name ] }
        } else {
            ch_ref_fasta_names = PREPARE_REFERENCES.out.ensembl_fasta_source_url
                .map { meta, urls_file ->
                    def names = urls_file.readLines().findAll { line -> line.trim() }
                        .collect { line -> line.tokenize('/').last() }.join(' + ')
                    [ meta.id.toString(), names ]
                }
                .mix(PREPARE_REFERENCES.out.ncbi_source_accessions
                    .map { meta, acc_file ->
                        def accs = acc_file.readLines()
                            .findAll { line -> line.trim() && !line.startsWith('stub:') }
                            .collect { line -> line.tokenize('/').last() }.join(', ')
                        [ meta.id.toString(), accs ?: meta.id.toString() ]
                    })
            ch_ref_gtf_names = PREPARE_REFERENCES.out.ensembl_gtf_source_urls
                .map { meta, urls_file ->
                    def name = urls_file.readLines().find { line -> line.trim() }?.tokenize('/')?.last() ?: ''
                    [ meta.id.toString(), name ]
                }
                .mix(PREPARE_REFERENCES.out.gtf_local
                    .map { meta, gtf_file -> [ meta.id.toString(), gtf_file.name ] })
                .mix(PREPARE_REFERENCES.out.ncbi_source_accessions
                    .map { meta, _acc -> [ meta.id.toString(), '' ] })
        }

        def ch_rdat_input = FOLD_STRUCTURES.out.fold_input
            .map { meta, xml -> [ meta.id.toString(), meta, xml ] }
            .combine(
                FOLD_STRUCTURES.out.structures
                    .map { meta, fold_dir -> [ meta.id.toString(), fold_dir ] },
                by: 0
            )
            .map { _key, fold_meta, xml, fold_dir ->
                def reference_key = resolveReferenceKey(fold_meta)
                [ reference_key.toString(), fold_meta, xml, fold_dir ]
            }
            .combine(ch_ref_fasta_names, by: 0)
            .combine(ch_ref_gtf_names, by: 0)
            .map { _ref_key, fold_meta, xml, fold_dir, fasta_name, gtf_name ->
                [ fold_meta + [ fasta_name: fasta_name, gtf_name: gtf_name ], xml, fold_dir ]
            }

        RNAFRAMEWORK_TORDAT (
            ch_rdat_input
        )
    }

    // Add RNAframework outputs to MultiQC input collection.
    ch_multiqc_files = ch_multiqc_files.mix(ch_rfcount_rc.collect { rc_file -> rc_file[1] })
    ch_multiqc_files = ch_multiqc_files.mix(ch_rfcount_plots.collect { plot_file -> plot_file[1] })
    ch_multiqc_files = ch_multiqc_files.mix(NORMALISE_REACTIVITIES.out.xml.collect { xml_file -> xml_file[1] })
    ch_multiqc_files = ch_multiqc_files.mix(NORMALISE_REACTIVITIES.out.plots.collect { plot_file -> plot_file[1] })
    ch_multiqc_files = ch_multiqc_files.mix(FOLD_STRUCTURES.out.structures.collect { fold_dir -> fold_dir[1] })

    // RF-norm summary table: one row per normalisation group (sample_group + replicate).
    // Join RF-count covered transcript counts per group (max across treated replicates).
    def ch_rfcount_covered_by_group = NORMALISE_REACTIVITIES.out.norm_groups
        .filter  { _group, condition, _meta, _rc, _rci -> condition == 'treated' }
        .map     { group, _condition, meta, _rc, _rci -> [ meta.id.toString(), group ] }
        .join(ch_rfcount_covered_transcripts)
        .map     { _sample_id, group, covered -> [ group, covered as long ] }
        .groupTuple()
        .map     { group, covered_list -> [ group, covered_list.max() ] }

    def ch_rfnorm_stats_mqc = NORMALISE_REACTIVITIES.out.rfnorm_log
        .map { meta, log -> [ meta.id.toString(), parseRfnormLog(log) ] }
        .join(ch_rfcount_covered_by_group, remainder: true)
        .map { group_id, rfnorm_stats, rfcount_covered ->
            [ group_id, [ rfcount_covered: (rfcount_covered ?: 0L) as long, covered: rfnorm_stats?.covered ?: 0L ] ]
        }
        .collect()
        .map { rows -> rfnormStatsMultiqc(rows) }

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_rfnorm_stats_mqc.collectFile(name: 'rfnorm_stats_mqc.yaml', sort: true)
    )

    // RF-fold summary table: one row per fold group (sample_group; may span replicates).
    // Skipped when stop_after_jackknife is true (rffold_log is empty).
    if (!params.stop_after_jackknife) {
        // Sum of per-replicate rf-norm covered transcripts, keyed by fold group (sample_group).
        def ch_norm_covered_by_foldgroup = NORMALISE_REACTIVITIES.out.rfnorm_log
            .map { meta, log -> [ meta.sample_group.toString(), parseRfnormLog(log).covered as long ] }
            .groupTuple()
            .map { sample_group, covered_list -> [ sample_group, covered_list.sum() as long ] }

        // "discarded" = transcripts normalised in some replicate but dropped from the consensus fold
        // because they were not present in every replicate: sum(per-rep covered) - n_reps * folded.
        def ch_rffold_stats_mqc = FOLD_STRUCTURES.out.rffold_log
            .map { meta, log -> [ meta.id.toString(), (meta.fold_experiment_count ?: 1) as long, parseRffoldLog(log).folded as long ] }
            .join(ch_norm_covered_by_foldgroup)
            .map { group_id, n_reps, folded, norm_covered_sum ->
                def discarded = Math.max(0L, (norm_covered_sum as long) - (n_reps * folded))
                [ group_id, [ folded: folded, discarded: discarded ] ]
            }
            .collect()
            .map { rows -> rffoldStatsMultiqc(rows) }

        ch_multiqc_files = ch_multiqc_files.mix(
            ch_rffold_stats_mqc.collectFile(name: 'rffold_stats_mqc.yaml', sort: true)
        )
    }

    // RF-eval summary table: one row per rfnorm group. Only present when reference structures
    // were supplied, so the section is absent from the report on a normal run.
    if (params.rfeval_reference) {
        def ch_rfeval_stats_mqc = FOLD_STRUCTURES.out.rfeval_metrics
            .map { meta, tsv -> [ meta.id.toString(), parseRfevalMetrics(tsv) ] }
            .collect()
            .map { rows -> rfevalStatsMultiqc(rows) }

        ch_multiqc_files = ch_multiqc_files.mix(
            ch_rfeval_stats_mqc.collectFile(name: 'rfeval_stats_mqc.yaml', sort: true)
        )
    }

    // rf-correlate — replicate-reproducibility QC (see CORRELATE_REPLICATES subworkflow).
    if (params.correlate_replicates) {
        CORRELATE_REPLICATES(FOLD_STRUCTURES.out.fold_input)
        ch_multiqc_files = ch_multiqc_files.mix(CORRELATE_REPLICATES.out.multiqc)
    }


    // MODULE: multiqc — aggregate pipeline quality control reports
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    summary_params      = filterSummaryParams(summary_params, transcriptome)
    summary_params      = addModuleOptionsSummary(summary_params, transcriptome)
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description, transcriptome))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect().map { f -> [f] }
            .combine(ch_multiqc_config.mix(ch_multiqc_custom_config).collect().map { c -> [c] })
            .combine(ch_multiqc_logo.collect().ifEmpty([]).map { l -> [l] })
            .map { files, config, logo -> [ [ id: 'multiqc' ], files, config, logo, [], [] ] }
    )

    // Collate and save software versions. Every module — nf-core and local alike — reports via the
    // `versions` topic, so there is nothing to mix in per-process; softwareVersionsToYAML is called
    // with an empty channel purely to contribute the workflow/Nextflow version block.
    def topic_versions_string = channel.topic("versions")
        .distinct()
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(channel.empty())
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'rnastructurome_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )


    emit:
    multiqc_report   = MULTIQC.out.report.toList()             // channel: /path/to/multiqc_report.html
    mapped_bam       = ch_dedup_bam                            // channel: [ val(meta), path(bam) ]
    normalized_xml   = NORMALISE_REACTIVITIES.out.xml           // channel: [ val(meta), path(xml) ]
    jackknife_csv    = FOLD_STRUCTURES.out.jackknife_csv        // channel: [ val(meta), path(csv) ] — empty when --jackknife_reference not set
    rfeval_metrics   = FOLD_STRUCTURES.out.rfeval_metrics       // channel: [ val(meta), path(tsv) ] — empty when --rfeval_reference not set
    fold_structures  = FOLD_STRUCTURES.out.structures           // channel: [ val(meta), path(dir) ]
    fold_bp          = FOLD_STRUCTURES.out.bp_dotplot           // channel: [ val(meta), path(bp) ]
    merged_bp        = FOLD_STRUCTURES.out.bp                  // channel: [ val(meta), path(*_merged.bp) ]

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
