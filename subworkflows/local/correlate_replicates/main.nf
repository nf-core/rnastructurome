// CORRELATE_REPLICATES — rf-correlate replicate-reproducibility QC. Reuses FOLD_STRUCTURES' per-sample_group
// XML grouping, runs only for groups with >1 replicate, and reports both Pearson (reactivity-capped) and
// Spearman correlations in a single combined MultiQC table.

include { RNAFRAMEWORK_RFCORRELATE as RNAFRAMEWORK_RFCORRELATE_PEARSON  } from '../../../modules/local/rnaframework/correlate/main'
include { RNAFRAMEWORK_RFCORRELATE as RNAFRAMEWORK_RFCORRELATE_SPEARMAN } from '../../../modules/local/rnaframework/correlate/main'
include {
    parseRfcorrelateMatrix
    rfCorrelateMultiqc
} from '../utils_nfcore_rnastructurome_pipeline/main'

workflow CORRELATE_REPLICATES {

    take:
    ch_fold_input  // channel: [ val(meta), list(path(xml)) ] — grouped by sample_group (FOLD_STRUCTURES.out.fold_input)

    main:

    def ch_correlate_input = ch_fold_input
        .filter { meta, _xmls -> (meta.fold_replicate_sizes?.size() ?: 0) >= 2 }
        .map { meta, xmls ->
            // rf-correlate sample labels must be unique and free of ':' / whitespace; derive from
            // the per-replicate labels, sanitised, with an index suffix to guarantee uniqueness.
            def rawLabels = (meta.fold_replicates ?: '').toString().tokenize(',')
            def labels    = rawLabels.withIndex().collect { lbl, i -> "rep${i + 1}_${lbl.replaceAll(/[^A-Za-z0-9_.-]/, '_')}" }
            def cmeta = meta + [
                correlate_replicate_sizes : meta.fold_replicate_sizes,
                correlate_replicate_labels: labels
            ]
            [ cmeta, xmls ]
        }

    RNAFRAMEWORK_RFCORRELATE_PEARSON(ch_correlate_input)
    RNAFRAMEWORK_RFCORRELATE_SPEARMAN(ch_correlate_input)

    // Summarise both methods' overall pairwise correlations into one combined MultiQC table.
    // Join per sample_group so each row carries mean Pearson + mean Spearman side by side.
    def ch_pearson_by_id = RNAFRAMEWORK_RFCORRELATE_PEARSON.out.matrix
        .map { meta, matrix -> [ meta.id.toString(), parseRfcorrelateMatrix(matrix) ] }
    def ch_spearman_by_id = RNAFRAMEWORK_RFCORRELATE_SPEARMAN.out.matrix
        .map { meta, matrix -> [ meta.id.toString(), parseRfcorrelateMatrix(matrix) ] }

    def ch_multiqc = ch_pearson_by_id
        .join(ch_spearman_by_id)
        .map { id, pearson, spearman ->
            [ id, [
                replicates    : pearson.replicates ?: spearman.replicates ?: 0,
                mean_pearson  : pearson.mean_corr,
                mean_spearman : spearman.mean_corr
            ] ]
        }
        .collect()
        .map { rows -> rfCorrelateMultiqc(rows) }
        .collectFile(name: 'rfcorrelate_mqc.yaml', sort: true)

    emit:
    matrix_pearson  = RNAFRAMEWORK_RFCORRELATE_PEARSON.out.matrix  // channel: [ val(meta), path(matrix.csv) ]
    matrix_spearman = RNAFRAMEWORK_RFCORRELATE_SPEARMAN.out.matrix // channel: [ val(meta), path(matrix.csv) ]
    multiqc  = ch_multiqc                                          // channel: path(rfcorrelate_mqc.yaml)
}
