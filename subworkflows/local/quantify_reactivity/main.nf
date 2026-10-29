// QUANTIFY_REACTIVITY — per-base RT-stop/mutation counts as RNAFramework RC files.
// Genome route: rf-count on STAR's TranscriptomeSAM output (default) or rf-count-genome + rf-rctools
// extract (count_genome=true). Transcriptome route: rf-count directly on the transcript BAM.

include { RNAFRAMEWORK_RFCOUNT           } from '../../../modules/local/rnaframework/count/main'
include { RNAFRAMEWORK_RFCOUNT as RNAFRAMEWORK_RFCOUNT_STAR } from '../../../modules/local/rnaframework/count/main'
include { RNAFRAMEWORK_RFCOUNT_GENOME    } from '../../../modules/local/rnaframework/count_genome/main'
include { RNAFRAMEWORK_RFRCTOOLS_EXTRACT } from '../../../modules/local/rnaframework/rctools/extract/main'

include { collectToMap        } from '../utils_nfcore_rnastructurome_pipeline/main'
include { isFourBaseDms       } from '../utils_nfcore_rnastructurome_pipeline/main'
include { resolveReferenceKey } from '../utils_nfcore_rnastructurome_pipeline/main'

workflow QUANTIFY_REACTIVITY {

    take:
    ch_markdup_bam_bai             // channel: [ val(meta), path(bam), path(bai) ]
    ch_transcript_bam_bai          // channel: [ val(meta), path(bam), path(bai) ] (genome route, count_genome=false)
    ch_strandedness_by_id          // channel: [ id, strandedness ]
    ch_reference_genome_fasta_keyed // channel: [ key, [meta, fasta] ]
    ch_reference_gtf_map           // value:   map ref_key -> [meta, gtf]
    ch_reference_fasta_map         // value:   map ref_key -> [meta, fasta]
    ch_genome_transcript_fasta_map // value:   map ref_key -> [meta, fasta] (genome route, count_genome=false)
    transcriptome                // boolean: transcriptome (Bowtie) route, including auto-detection

    main:

    def ch_rfcount_rc       = channel.empty()
    def ch_rfcount_rci       = channel.empty()
    def ch_rfcount_summary   = channel.empty()
    def ch_rfcount_plots     = channel.empty()

    // Annotate each BAM with the per-sample strandedness inferred by RSeQC. remainder: true keeps
    // samples with no BED (e.g. viral) at strandedness = null, falling back to params.rfcount_strandedness.
    // Used by both genome-route branches below (rf-count-genome and rf-count-direct).
    def ch_stranded_meta_by_id = ch_strandedness_by_id

    if (!transcriptome && params.count_genome) {
        def ch_genome_fasta_map = collectToMap(ch_reference_genome_fasta_keyed)

        def ch_bam_stranded = ch_markdup_bam_bai
            .map { meta, bam, bai -> [ meta.id.toString(), meta, bam, bai ] }
            .join(ch_stranded_meta_by_id, remainder: true)
            .map { _id, meta, bam, bai, strandedness ->
                [ meta + [strandedness: strandedness], bam, bai ]
            }

        def ch_rfcount_genome_inputs = ch_bam_stranded
            .combine(ch_genome_fasta_map)
            .map { combined ->
                def meta    = combined[0]
                def bam     = combined[1]
                def bai     = combined[2]
                def ref_map = combined[3]
                def ref_key = resolveReferenceKey(meta)
                def fasta_t = ref_map[ref_key]
                if (!fasta_t) error("No genome FASTA resolved for reference '${ref_key}' for rf-count-genome.")
                [ [meta + [four_base_dms: isFourBaseDms(meta)], bam, bai], fasta_t ]
            }
        def ch_rg_split = ch_rfcount_genome_inputs.multiMap { entry ->
            bam:   entry[0]
            fasta: entry[1]
        }
        RNAFRAMEWORK_RFCOUNT_GENOME(ch_rg_split.bam, ch_rg_split.fasta)
        ch_rfcount_plots   = RNAFRAMEWORK_RFCOUNT_GENOME.out.plots

        // rf-rctools extract: genome RC → transcript-level RC using the reference GTF, calling extract
        // with the BASENAME for strand-aware extraction. Also rewrites the summary's genome-level
        // 'covered' count to the covered-transcript count, so pair the RC with its summary (strict join).
        def ch_rc_with_summary = RNAFRAMEWORK_RFCOUNT_GENOME.out.rc
            .map { meta, rc -> [ meta.id.toString(), meta, rc ] }
            .join(RNAFRAMEWORK_RFCOUNT_GENOME.out.summary.map { meta, summary -> [ meta.id.toString(), summary ] })
            .map { _id, meta, rc, summary -> [ meta, rc, summary ] }

        def ch_rctools_inputs = ch_rc_with_summary
            .combine(ch_reference_gtf_map)
            .map { combined ->
                def meta    = combined[0]
                def rc      = combined[1]
                def summary = combined[2]
                def gtf_map = combined[3]
                def ref_key = resolveReferenceKey(meta)
                def gtf_t   = gtf_map[ref_key]
                if (!gtf_t) error("No GTF resolved for reference '${ref_key}' for rf-rctools extract.")
                [ [meta, rc, [], summary], gtf_t ]
            }
        def ch_rct_split = ch_rctools_inputs.multiMap { entry ->
            rc:  entry[0]
            gtf: entry[1]
        }
        RNAFRAMEWORK_RFRCTOOLS_EXTRACT(
            ch_rct_split.rc,
            ch_rct_split.gtf
        )
        ch_rfcount_rc      = RNAFRAMEWORK_RFRCTOOLS_EXTRACT.out.rc
        ch_rfcount_rci     = RNAFRAMEWORK_RFRCTOOLS_EXTRACT.out.rci
        ch_rfcount_summary = RNAFRAMEWORK_RFRCTOOLS_EXTRACT.out.summary
    } else if (!transcriptome) {
        // Default genome route: rf-count directly on the dedup-reconciled, calmd-corrected
        // transcript-coordinate BAM from ALIGN_READS — no rf-rctools extract needed.
        def ch_transcript_bam_stranded = ch_transcript_bam_bai
            .map { meta, bam, bai -> [ meta.id.toString(), meta, bam, bai ] }
            .join(ch_stranded_meta_by_id, remainder: true)
            .map { _id, meta, bam, bai, strandedness ->
                [ meta + [strandedness: strandedness], bam, bai ]
            }

        def ch_rfcount_direct_inputs = ch_transcript_bam_stranded
            .combine(ch_genome_transcript_fasta_map)
            .map { combined ->
                def meta    = combined[0]
                def bam     = combined[1]
                def bai     = combined[2]
                def ref_map = combined[3]
                def ref_key = resolveReferenceKey(meta)
                def fasta_t = ref_map[ref_key]
                if (!fasta_t) error("No transcript FASTA resolved for reference '${ref_key}' for rf-count.")
                [ [meta + [four_base_dms: isFourBaseDms(meta)], bam, bai], fasta_t ]
            }
        def ch_rd_split = ch_rfcount_direct_inputs.multiMap { entry ->
            bam:   entry[0]
            fasta: entry[1]
        }
        RNAFRAMEWORK_RFCOUNT_STAR(ch_rd_split.bam, ch_rd_split.fasta)
        ch_rfcount_rc      = RNAFRAMEWORK_RFCOUNT_STAR.out.rc
        ch_rfcount_rci     = RNAFRAMEWORK_RFCOUNT_STAR.out.rci
        ch_rfcount_summary = RNAFRAMEWORK_RFCOUNT_STAR.out.summary
        ch_rfcount_plots   = RNAFRAMEWORK_RFCOUNT_STAR.out.plots
    } else {
        def ch_rfcount_inputs = ch_markdup_bam_bai
            .combine(ch_reference_fasta_map)
            .map { combined ->
                def meta    = combined[0]
                def bam     = combined[1]
                def bai     = combined[2]
                def ref_map = combined[3]
                def ref_key = resolveReferenceKey(meta)
                def fasta_t = ref_map[ref_key]
                if (!fasta_t) error("No transcript FASTA resolved for reference '${ref_key}' for rf-count.")
                [ [meta + [four_base_dms: isFourBaseDms(meta)], bam, bai], fasta_t ]
            }
        def ch_rc_split = ch_rfcount_inputs.multiMap { entry ->
            bam:   entry[0]
            fasta: entry[1]
        }
        RNAFRAMEWORK_RFCOUNT(ch_rc_split.bam, ch_rc_split.fasta)
        ch_rfcount_rc      = RNAFRAMEWORK_RFCOUNT.out.rc
        ch_rfcount_rci     = RNAFRAMEWORK_RFCOUNT.out.rci
        ch_rfcount_summary = RNAFRAMEWORK_RFCOUNT.out.summary
        ch_rfcount_plots   = RNAFRAMEWORK_RFCOUNT.out.plots
    }

    emit:
    rc       = ch_rfcount_rc       // channel: [ val(meta), path(rc) ]
    rci      = ch_rfcount_rci      // channel: [ val(meta), path(rci) ]
    summary  = ch_rfcount_summary  // channel: [ val(meta), path(summary_tsv) ]
    plots    = ch_rfcount_plots    // channel: [ val(meta), path(plots) ]
}
