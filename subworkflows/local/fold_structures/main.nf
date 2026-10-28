// FOLD_STRUCTURES — group rfnorm XML files by sample_group, optionally run rf-jackknife for slope/intercept
// calibration, run rf-fold, then convert dot-plot outputs to .bp arc files (genome+transcript) and merge.

include { RNAFRAMEWORK_RFJACKKNIFE                                     } from '../../../modules/local/rnaframework/jackknife/main'
include { RNAFRAMEWORK_RFEVAL                                          } from '../../../modules/local/rnaframework/eval/main'
include { RNAFRAMEWORK_RFFOLD                                          } from '../../../modules/local/rnaframework/fold/main'
include { RNAFRAMEWORK_RFSTRUCTEXTRACT                                 } from '../../../modules/local/rnaframework/structextract/main'
include { RNAFRAMEWORK_DOTPLOT2BP                                      } from '../../../modules/local/dotplot2bp/main'
include { RNAFRAMEWORK_DOTPLOT2BP as RNAFRAMEWORK_DOTPLOT2BP_TRANSCRIPT } from '../../../modules/local/dotplot2bp/main'
include { MERGE_BP                                                     } from '../../../modules/local/merge_bp/main'
include { MERGE_BP as MERGE_BP_TRANSCRIPT                              } from '../../../modules/local/merge_bp/main'
include { resolveReferenceKey } from '../utils_nfcore_rnastructurome_pipeline/main'

workflow FOLD_STRUCTURES {

    take:
    ch_rfnorm_xml        // channel: [ val(meta), path(xml) ]  — one entry per rfnorm group
    ch_reference_gtf_map // value:   map of reference_key -> [ val(meta), path(gtf) ]

    main:

    // Group rfnorm XMLs by sample_group so biological replicates are folded together.
    // The group key intentionally excludes replicate.
    def ch_fold_input = ch_rfnorm_xml
        .map { meta, xml ->
            if (!meta.sample_group) {
                error("Missing sample_group for sample '${meta.id}'. rf-fold replicate grouping requires sample_group.")
            }
            def fold_group = meta.sample_group.toString()
            [ fold_group, [ meta, xml ] ]
        }
        .groupTuple()
        .map { fold_group, entries ->
            // rf-fold folds replicates via majority voting, taking one experiment DIRECTORY per replicate
            // (not a merged dir). Sort by replicate, then concatenate XMLs so the module can rebuild
            // per-replicate dirs from input*/ staging using the aligned fold_replicate_sizes counts.
            def sortedEntries  = entries.sort { entry -> (entry[0].replicate ?: 'na').toString() }
            def base           = sortedEntries[0][0]
            def perRepXmls     = sortedEntries.collect { entry -> [ entry[1] ].flatten() }
            def orderedXmls    = perRepXmls.flatten()
            def replicateSizes = perRepXmls.collect { repXmls -> repXmls.size() }
            def replicates     = sortedEntries.collect { entry -> (entry[0].replicate ?: 'na').toString() }
            def sampleIds      = sortedEntries.collect { entry -> (entry[0].id ?: 'na').toString() }
            def foldMeta   = base + [
                id                    : fold_group,
                fold_group            : fold_group,
                fold_replicates       : replicates.join(','),
                fold_source_ids       : sampleIds.join(','),
                fold_experiment_count : sortedEntries.size(),   // replicate experiments; drives rf-fold -oc
                fold_replicate_sizes  : replicateSizes          // XMLs per replicate, aligned to orderedXmls
            ]
            [ foldMeta, orderedXmls ]
        }

    // Optional rf-jackknife — quality assessment against a reference structure set; rf-fold is gated on
    // its completion. With --rfjackknife_pool_all, the optimal slope/intercept is injected into fold meta.
    def ch_fold_for_rffold = ch_fold_input
    def ch_jackknife_csv   = channel.empty()
    def ch_rfeval_metrics  = channel.empty()

    if (params.jackknife_reference) {
        def ch_jackknife_reference = channel.value(file(params.jackknife_reference.toString(), checkIfExists: true))

        // Jackknife runs per rfnorm group (sample_group + replicate), not per fold group, since fold
        // groups flatten identically-named XMLs (e.g. 16S_rRNA.xml) across replicates.
        def ch_jackknife_input
        if (params.rfjackknife_pool_all as Boolean) {
            ch_jackknife_input = ch_rfnorm_xml
                .flatMap { _meta, xmls -> [xmls].flatten() }
                .collect()
                .map { all_xmls -> [ [ id: 'all_groups', fold_group: 'all_groups' ], all_xmls ] }
        } else {
            ch_jackknife_input = ch_rfnorm_xml
        }

        RNAFRAMEWORK_RFJACKKNIFE (
            ch_jackknife_input,
            ch_jackknife_reference
        )
        ch_jackknife_csv = RNAFRAMEWORK_RFJACKKNIFE.out.csv

        if (params.rfjackknife_pool_all as Boolean) {
            // Parse optimal slope/intercept from the pooled jackknife CSV (FMI.csv is a semicolon-delimited
            // matrix: rows = slopes, columns = intercepts, cells = FMI) and inject into fold meta.
            def ch_calibration = RNAFRAMEWORK_RFJACKKNIFE.out.csv
                .map { _meta, csv_files ->
                    def csv_file = [csv_files].flatten()[0]
                    def lines = csv_file.readLines()
                    def intercepts = lines[0].split(';').drop(1)*.trim()
                    def best_slope = null; def best_intercept = null; def best_fmi = -1d
                    lines.drop(1).findAll { line -> line.trim() }.each { line ->
                        def parts = line.split(';')*.trim()
                        def slope = parts[0]
                        parts.drop(1).eachWithIndex { fmi_str, i ->
                            def fmi = fmi_str.toDouble()
                            if (fmi > best_fmi) {
                                best_fmi = fmi; best_slope = slope; best_intercept = intercepts[i]
                            }
                        }
                    }
                    [ best_slope, best_intercept ]
                }
            ch_fold_for_rffold = ch_fold_input
                .combine(ch_calibration)
                .map { meta, xmls, slope, intercept ->
                    [ meta + [ jackknife_slope: slope, jackknife_intercept: intercept ], xmls ]
                }
        } else {
            // Gate each fold group on all jackknife jobs for that sample_group completing.
            ch_fold_for_rffold = ch_fold_input
                .map { meta, xmls -> [ meta.id.toString(), meta, xmls ] }
                .join(
                    RNAFRAMEWORK_RFJACKKNIFE.out.csv
                        .map { meta, _csv -> [ meta.sample_group.toString(), 'done' ] }
                        .groupTuple()
                        .map { sample_group, _dones -> [ sample_group, 'done' ] }
                )
                .map { _id, meta, xmls, _done -> [ meta, xmls ] }
        }
    }

    // Optional rf-eval — evaluate agreement between reactivity data and a reference structure set.
    // Runs per rfnorm group when --rfeval_reference is provided.
    if (params.rfeval_reference) {
        def ch_rfeval_reference = channel.value(file(params.rfeval_reference.toString(), checkIfExists: true))

        // With --rfeval_windows, rf-eval first slices each group's XML reactivities to the manifest
        // regions, so a sub-region reference (e.g. an Rfam element on a whole chromosome) is scored
        // against a matching windowed XML instead of the diluted full transcript.
        def ch_rfeval_windows = params.rfeval_windows
            ? channel.value(file(params.rfeval_windows.toString(), checkIfExists: true))
            : channel.value([])

        RNAFRAMEWORK_RFEVAL (
            ch_rfnorm_xml,
            ch_rfeval_reference,
            ch_rfeval_windows
        )
        ch_rfeval_metrics = RNAFRAMEWORK_RFEVAL.out.metrics
    }

    def ch_fold_structures  = channel.empty()
    def ch_shannon_wig      = channel.empty()
    def ch_rffold_log       = channel.empty()
    def ch_bp_dotplot       = channel.empty()
    def ch_bp               = channel.empty()
    def ch_bp_transcript    = channel.empty()
    def ch_structextract    = channel.empty()

    if (!params.stop_after_jackknife) {
        RNAFRAMEWORK_RFFOLD (
            ch_fold_for_rffold
        )
        ch_fold_structures = RNAFRAMEWORK_RFFOLD.out.structures
        ch_shannon_wig     = RNAFRAMEWORK_RFFOLD.out.shannon_wig
        ch_rffold_log      = RNAFRAMEWORK_RFFOLD.out.log

        // Optional rf-structextract — pull high-confidence, low-reactivity/low-Shannon motifs out of the
        // rf-fold output, pairing each fold dir (-ro) with the group's XML reactivities (-xf), deduped by filename.
        if (params.structextract) {
            def ch_structextract_input = RNAFRAMEWORK_RFFOLD.out.structures
                .map { meta, fold_dir -> [ meta.id.toString(), meta, fold_dir ] }
                .join(
                    ch_fold_for_rffold.map { meta, xmls -> [ meta.id.toString(), xmls ] }
                )
                .map { _id, meta, fold_dir, xmls ->
                    def seen    = [] as Set
                    def deduped = [xmls].flatten().findAll { x -> seen.add(x.name) }
                    [ meta, fold_dir, deduped ]
                }

            RNAFRAMEWORK_RFSTRUCTEXTRACT(
                ch_structextract_input
            )
            ch_structextract = RNAFRAMEWORK_RFSTRUCTEXTRACT.out.motifs
        }

        // Resolve each fold group's GTF (absent for NCBI/viral references). The genome-coordinate
        // conversion needs it; the transcript-coordinate one does not — so only the genome variant is
        // gated on GTF presence, while the transcript variant runs for every reference (gtf = []).
        def ch_dotplot_bp_resolved = RNAFRAMEWORK_RFFOLD.out.structures
            .combine(ch_reference_gtf_map)
            .map { meta, fold_dir, gtf_map ->
                def reference_key = resolveReferenceKey(meta)
                def gtf_tuple = gtf_map[reference_key]
                if (!gtf_tuple) {
                    log.warn("No GTF for '${reference_key}': skipping genome-coordinate bp (transcript-coordinate bp still produced).")
                }
                [ meta, fold_dir, gtf_tuple ? gtf_tuple[1] : [] ]
            }

        def ch_dotplot_bp_genome = ch_dotplot_bp_resolved.filter { _meta, _fold_dir, gtf -> gtf }

        RNAFRAMEWORK_DOTPLOT2BP (
            ch_dotplot_bp_genome
        )

        RNAFRAMEWORK_DOTPLOT2BP_TRANSCRIPT (
            ch_dotplot_bp_resolved
        )

        MERGE_BP (
            RNAFRAMEWORK_DOTPLOT2BP.out.bp
        )

        MERGE_BP_TRANSCRIPT (
            RNAFRAMEWORK_DOTPLOT2BP_TRANSCRIPT.out.bp
        )


        ch_bp_dotplot    = RNAFRAMEWORK_DOTPLOT2BP.out.bp
        ch_bp            = MERGE_BP.out.bp
        ch_bp_transcript = MERGE_BP_TRANSCRIPT.out.bp
    }

    emit:
    fold_input    = ch_fold_input       // channel: [ val(meta), list(path(xml)) ] — grouped by sample_group
    structures    = ch_fold_structures  // channel: [ val(meta), path(fold_dir) ]  — empty when stop_after_jackknife
    shannon_wig   = ch_shannon_wig      // channel: [ val(meta), path(wig) ]        — empty when stop_after_jackknife
    rffold_log    = ch_rffold_log       // channel: [ val(meta), path(log) ]        — empty when stop_after_jackknife
    jackknife_csv = ch_jackknife_csv    // channel: [ val(meta), path(csv) ]        — empty when --jackknife_reference not set
    rfeval_metrics = ch_rfeval_metrics  // channel: [ val(meta), path(tsv) ]        — empty when --rfeval_reference not set
    bp_dotplot    = ch_bp_dotplot       // channel: [ val(meta), path(bp) ]         — empty when stop_after_jackknife
    bp            = ch_bp              // channel: [ val(meta), path(bp) ]          — empty when stop_after_jackknife
    bp_transcript = ch_bp_transcript   // channel: [ val(meta), path(bp) ]          — empty when stop_after_jackknife
    structextract = ch_structextract   // channel: [ val(meta), path(dir) ]         — empty unless --structextract
}
