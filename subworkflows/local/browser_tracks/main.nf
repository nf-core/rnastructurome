// BROWSER_TRACKS — turn rf-norm reactivities (and optional rf-fold Shannon entropy) into genome- and
// transcript-coordinate BigWig tracks for IGV/UCSC. rf-wiggle → merge/average replicates → wigToBigWig.

include { RNAFRAMEWORK_RFWIGGLE } from '../../../modules/local/rnaframework/wiggle/main'
include { MERGE_WIG             } from '../../../modules/local/merge_wig/main'
include { AVERAGE_WIG           } from '../../../modules/local/average_wig/main'
include { MERGE_SHANNON_WIG     } from '../../../modules/local/merge_shannon_wig/main'
include { WIG_TO_GENOME as WIG_TO_GENOME_REACTIVITY } from '../../../modules/local/wig_to_genome/main'
include { WIG_TO_GENOME as WIG_TO_GENOME_SHANNON    } from '../../../modules/local/wig_to_genome/main'
include { WIG_CHROM_SIZES as WIG_CHROM_SIZES_REACTIVITY } from '../../../modules/local/wig_chrom_sizes/main'
include { WIG_CHROM_SIZES as WIG_CHROM_SIZES_SHANNON    } from '../../../modules/local/wig_chrom_sizes/main'
include { UCSC_WIGTOBIGWIG as UCSC_WIGTOBIGWIG_REACTIVITY            } from '../../../modules/nf-core/ucsc/wigtobigwig/main'
include { UCSC_WIGTOBIGWIG as UCSC_WIGTOBIGWIG_SHANNON               } from '../../../modules/nf-core/ucsc/wigtobigwig/main'
include { UCSC_WIGTOBIGWIG as UCSC_WIGTOBIGWIG_TRANSCRIPT_REACTIVITY } from '../../../modules/nf-core/ucsc/wigtobigwig/main'
include { UCSC_WIGTOBIGWIG as UCSC_WIGTOBIGWIG_TRANSCRIPT_SHANNON    } from '../../../modules/nf-core/ucsc/wigtobigwig/main'
include { resolveReferenceKey } from '../utils_nfcore_rnastructurome_pipeline/main'

workflow BROWSER_TRACKS {

    take:
    ch_rfnorm_xml          // channel: [ val(meta), path(xml) ]
    ch_fold_shannon_wig    // channel: [ val(meta), path(shannon_wig) ]
    ch_reference_gtf_map   // value:   map of reference_key -> [meta, gtf]

    main:

    // MODULE: rf-wiggle — convert rf-norm XML reactivities to WIG + chrom.sizes
    RNAFRAMEWORK_RFWIGGLE (
        ch_rfnorm_xml
    )

    MERGE_WIG (
        RNAFRAMEWORK_RFWIGGLE.out.wig
    )

    // Group per-replicate merged WIGs by sample_group. Single replicate bypasses AVERAGE_WIG
    // (meta.id stays e.g. "HEK293T_1"); multiple replicates run it with meta.id = sample_group.
    def ch_reactivity_grouped = MERGE_WIG.out.merged_wig
        .map { meta, wig -> [ meta.sample_group.toString(), meta, wig ] }
        .groupTuple(by: 0)
        .map { sample_group, metas, wigs ->
            def base_meta = wigs.size() == 1
                ? metas[0]
                : metas[0] + [ id: sample_group ]
            [ base_meta, wigs.flatten() ]
        }

    def ch_reactivity_branches = ch_reactivity_grouped.branch { _meta, wigs ->
        single: wigs.size() == 1
        multi:  true
    }

    AVERAGE_WIG (
        ch_reactivity_branches.multi
    )

    // MODULE: wig_to_genome + wigToBigWig — convert merged/averaged transcript WIG to genomic BigWig for IGV
    def ch_reactivity_for_bigwig = ch_reactivity_branches.single
        .map { meta, wigs -> [ meta, wigs[0] ] }
        .mix(AVERAGE_WIG.out.merged_wig)
    def ch_reactivity_genomic_wig_input = ch_reactivity_for_bigwig
        .combine(ch_reference_gtf_map)
        .flatMap { combined ->
            def meta = combined[0]
            def wig = combined[1]
            def gtf_map = combined[2]
            def reference_key = resolveReferenceKey(meta)
            def gtf_tuple = gtf_map[reference_key]
            if (!gtf_tuple) {
                log.warn("Skipping genomic reactivity BigWig for '${reference_key}': no GTF available.")
                return []
            }
            return [ [ meta, wig, gtf_tuple[1] ] ]
        }

    WIG_TO_GENOME_REACTIVITY (
        ch_reactivity_genomic_wig_input
    )

    def ch_reactivity_genome_split = WIG_TO_GENOME_REACTIVITY.out.wig
        .map { meta, wig -> [ meta.id.toString(), meta, wig ] }
        .join(WIG_TO_GENOME_REACTIVITY.out.chrom_sizes.map { meta, sizes -> [ meta.id.toString(), sizes ] })
        .filter { _id, _meta, _wig, sizes -> sizes.size() > 0 }
        .multiMap { _id, meta, wig, sizes ->
            wig:   [ meta, wig ]
            sizes: sizes
        }

    UCSC_WIGTOBIGWIG_REACTIVITY (
        ch_reactivity_genome_split.wig,
        ch_reactivity_genome_split.sizes
    )

    // MODULE: WIG_CHROM_SIZES + wigToBigWig — transcript-coordinate reactivity BigWig. Independent
    // .map{} subscriptions on the same sources ensure items aren't consumed by the genome operators above.
    def ch_reactivity_wig_for_sizes = ch_reactivity_branches.single
        .map { meta, wigs -> [ meta, wigs[0] ] }
        .mix(AVERAGE_WIG.out.merged_wig)
    def ch_reactivity_wig_for_transcript_bw = ch_reactivity_branches.single
        .map { meta, wigs -> [ meta, wigs[0] ] }
        .mix(AVERAGE_WIG.out.merged_wig)

    WIG_CHROM_SIZES_REACTIVITY(ch_reactivity_wig_for_sizes)

    def ch_reactivity_transcript_bw_split = ch_reactivity_wig_for_transcript_bw
        .map { meta, wig -> [ meta.id.toString(), meta, wig ] }
        .join(WIG_CHROM_SIZES_REACTIVITY.out.sizes.map { meta, sizes -> [ meta.id.toString(), sizes ] })
        .multiMap { _id, meta, wig, sizes ->
            wig:   [ meta, wig ]
            sizes: sizes
        }

    UCSC_WIGTOBIGWIG_TRANSCRIPT_REACTIVITY(
        ch_reactivity_transcript_bw_split.wig,
        ch_reactivity_transcript_bw_split.sizes
    )


    // MODULE: merge_shannon_wig + wigToBigWig — merge per-transcript Shannon entropy WIG files and
    // convert to genome- and transcript-coordinate BigWigs.
    if (params.rffold_shannon_entropy) {
        MERGE_SHANNON_WIG (
            ch_fold_shannon_wig
        )

        def ch_shannon_genomic_wig_input = MERGE_SHANNON_WIG.out.merged_wig
            .combine(ch_reference_gtf_map)
            .flatMap { combined ->
                def meta = combined[0]
                def wig = combined[1]
                def gtf_map = combined[2]
                def reference_key = resolveReferenceKey(meta)
                def gtf_tuple = gtf_map[reference_key]
                if (!gtf_tuple) {
                    log.warn("Skipping genomic Shannon BigWig for '${reference_key}': no GTF available.")
                    return []
                }
                return [ [ meta, wig, gtf_tuple[1] ] ]
            }

        WIG_TO_GENOME_SHANNON (
            ch_shannon_genomic_wig_input
        )

        def ch_shannon_genome_split = WIG_TO_GENOME_SHANNON.out.wig
            .map { meta, wig -> [ meta.id.toString(), meta, wig ] }
            .join(WIG_TO_GENOME_SHANNON.out.chrom_sizes.map { meta, sizes -> [ meta.id.toString(), sizes ] })
            .filter { _id, _meta, _wig, sizes -> sizes.size() > 0 }
            .multiMap { _id, meta, wig, sizes ->
                wig:   [ meta, wig ]
                sizes: sizes
            }

        UCSC_WIGTOBIGWIG_SHANNON (
            ch_shannon_genome_split.wig,
            ch_shannon_genome_split.sizes
        )

        // Transcript-coordinate Shannon BigWig — independent subscriptions from the genome path above.
        def ch_shannon_wig_for_sizes = MERGE_SHANNON_WIG.out.merged_wig
        def ch_shannon_wig_for_transcript_bw = MERGE_SHANNON_WIG.out.merged_wig

        WIG_CHROM_SIZES_SHANNON(ch_shannon_wig_for_sizes)

        def ch_shannon_transcript_bw_split = ch_shannon_wig_for_transcript_bw
            .map { meta, wig -> [ meta.id.toString(), meta, wig ] }
            .join(WIG_CHROM_SIZES_SHANNON.out.sizes.map { meta, sizes -> [ meta.id.toString(), sizes ] })
            .multiMap { _id, meta, wig, sizes ->
                wig:   [ meta, wig ]
                sizes: sizes
            }

        UCSC_WIGTOBIGWIG_TRANSCRIPT_SHANNON(
            ch_shannon_transcript_bw_split.wig,
            ch_shannon_transcript_bw_split.sizes
        )
    }
}
