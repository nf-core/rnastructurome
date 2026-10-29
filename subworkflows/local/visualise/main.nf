// VISUALISE_STRUCTURES — render 2D structure diagrams from rf-fold structures. ViennaRNA (RNAplot) always
// draws every structure; R2DT (container-only) additionally draws template-based diagrams in parallel.

include { R2DT                } from '../../../modules/local/r2dt/main'
include { VIENNARNA           } from '../../../modules/local/viennarna/main'
include { resolveReferenceKey } from '../utils_nfcore_rnastructurome_pipeline/main'

workflow VISUALISE_STRUCTURES {

    take:
    ch_rfnorm_xml          // channel: [ val(meta), path(xml) ]
    ch_fold_structures     // channel: [ val(meta), path(dir) ]
    ch_fold_input          // channel: [ val(meta), path(xmls) ]
    ch_reference_fasta_map // value:   map of reference_key -> [meta, fasta]
    ch_reference_gtf_map   // value:   map of reference_key -> [meta, gtf]

    main:

    // Always draw every structure with ViennaRNA, independent of R2DT (empty drawn-id list means
    // "draw all"), so the two renderers run in parallel rather than ViennaRNA filling R2DT's gaps.
    // The empty list is a real asset file, NOT file('/dev/null'): staging /dev/null makes Nextflow
    // bind-mount the host /dev into the container, which breaks /dev/null writability under Singularity.
    def ch_empty_drawn_ids = file("${projectDir}/assets/NO_DRAWN_IDS.txt", checkIfExists: true)
    def ch_rnaplot_input = ch_fold_structures
        .map { meta, dir -> [ meta.id.toString(), meta, dir ] }
        .join(ch_fold_input.map { meta, xmls -> [ meta.id.toString(), xmls ] })
        .map { _id, meta, dir, xmls -> [ meta, dir, xmls, ch_empty_drawn_ids ] }

    VIENNARNA(
        ch_rnaplot_input
    )

    // R2DT (container-only): additionally draw template-based diagrams for every
    // structure, published alongside the ViennaRNA renderings for comparison.
    if (params.r2dt) {
        def ch_r2dt_xml = ch_rfnorm_xml
            .map { meta, xml ->
                def fold_group = meta.sample_group?.toString()
                if (!fold_group) error("Missing sample_group for '${meta.id}' — required for R2DT grouping.")
                [ fold_group, xml instanceof List ? xml : [xml] ]
            }
            .groupTuple()
            .map { fold_group, xml_lists -> [ fold_group, xml_lists.flatten() ] }

        def ch_r2dt_input = ch_fold_structures
            .map { meta, dir -> [ meta.id.toString(), meta, dir ] }
            .join(ch_r2dt_xml)
            .map { _fg, meta, dir, xmls -> [ meta, dir, xmls ] }
            .combine(ch_reference_fasta_map)
            .combine(ch_reference_gtf_map)
            .flatMap { combined ->
                def meta      = combined[0]
                def fold_dir  = combined[1]
                def xmls      = combined[2]
                def fasta_map = combined[3]
                def gtf_map   = combined[4]
                def ref_key   = resolveReferenceKey(meta)
                def fasta_t   = fasta_map[ref_key]
                if (!fasta_t) {
                    log.warn("Skipping R2DT for '${meta.id}': no FASTA for '${ref_key}'")
                    return []
                }
                // No GTF (e.g. viral synthetic reference): draw all — the module skips biotype filtering.
                def gtf_t = gtf_map[ref_key]
                return [ [ meta, fold_dir, xmls, fasta_t[1], gtf_t ? gtf_t[1] : [] ] ]
            }

        R2DT(
            ch_r2dt_input
        )
    }
}
