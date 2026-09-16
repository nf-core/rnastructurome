// FASTQ_QC_TRIM — per-sample read intake, QC and adapter/UMI trimming: cat -> pre-trim FastQC -> optional
// UMI extraction -> principle-specific cutadapt -> post-trim FastQC. Emits post-cat samplesheet, trimmed
// reads per principle, and MultiQC contributions (versions come via the `versions` topic).

include { FASTQC as FASTQC_PRE        } from '../../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_POST       } from '../../../modules/nf-core/fastqc/main'
include { CAT_FASTQ                   } from '../../../modules/nf-core/cat/fastq/main'
include { CUTADAPT as CUTADAPT_RTSTOP } from '../../../modules/nf-core/cutadapt/main'
include { CUTADAPT as CUTADAPT_MAP    } from '../../../modules/nf-core/cutadapt/main'
include { UMITOOLS_EXTRACT            } from '../../../modules/nf-core/umitools/extract/main'

include { parseCutadaptCommandArg } from '../utils_nfcore_rnastructurome_pipeline/main'
include { cutadaptAdaptersMultiqc } from '../utils_nfcore_rnastructurome_pipeline/main'

workflow FASTQ_QC_TRIM {

    take:
    ch_samplesheet  // channel: [ val(meta), [ reads ] ]

    main:
    ch_multiqc_files = channel.empty()

    def ch_samplesheet_checked = ch_samplesheet.map { meta, reads ->
        def principle = (meta.principle ?: '').toLowerCase()
        if (!(principle in ['rt-stop', 'map'])) {
            error("Unsupported principle '${meta.principle}' for sample '${meta.id}'. Expected one of: RT-stop, MaP.")
        }
        [meta, reads]
    }

    // MODULE: cat/fastq — merge resequenced FASTQ files per sample before QC
    CAT_FASTQ (
        ch_samplesheet_checked
    )

    def ch_pretrim_reads_split = CAT_FASTQ.out.reads.multiMap { meta, reads ->
        fastqc: [ meta, reads ]
        branching: [ meta, reads ]
    }
    def ch_pretrim_fastqc_input = ch_pretrim_reads_split.fastqc.map { meta, reads -> [ meta, reads ] }
    def ch_samplesheet_for_branching = ch_pretrim_reads_split.branching.map { meta, reads -> [ meta, reads ] }

    // MODULE: fastqc (pre-trim) — quality control on raw reads
    FASTQC_PRE (
        ch_pretrim_fastqc_input
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_PRE.out.zip.collect { fastqc_zip -> fastqc_zip[1] })

    // Branch by probing principle so RT-stop and MaP can use different default cutadapt args
    def principle_branches = ch_samplesheet_for_branching.branch { meta, _reads ->
        rtstop: (meta.principle ?: '').toLowerCase() == 'rt-stop'
        map:    (meta.principle ?: '').toLowerCase() == 'map'
    }

    def ch_reads_for_umi = principle_branches.rtstop.mix(principle_branches.map)
    def ch_reads_with_umi = ch_reads_for_umi.filter { meta, _reads ->
        (meta.umi_pattern ?: '').toString().trim()
    }
    def ch_reads_without_umi = ch_reads_for_umi.filter { meta, _reads ->
        !((meta.umi_pattern ?: '').toString().trim())
    }

    // MODULE: umi_tools extract — extract UMIs from reads
    UMITOOLS_EXTRACT (
        ch_reads_with_umi
    )

    def ch_reads_after_umi = ch_reads_without_umi.mix(UMITOOLS_EXTRACT.out.reads)
    def reads_for_cutadapt = ch_reads_after_umi.branch { meta, _reads ->
        rtstop: (meta.principle ?: '').toLowerCase() == 'rt-stop'
        map:    (meta.principle ?: '').toLowerCase() == 'map'
    }
    def ch_rtstop_reads_for_cutadapt = reads_for_cutadapt.rtstop
    def ch_map_reads_for_cutadapt = reads_for_cutadapt.map
    ch_multiqc_files = ch_multiqc_files.mix(UMITOOLS_EXTRACT.out.log.collect { umi_log -> umi_log[1] })

    // MODULE: cutadapt (RT-stop) — trim RT-stop reads
    CUTADAPT_RTSTOP (
        ch_rtstop_reads_for_cutadapt
    )

    // MODULE: cutadapt (MaP) — trim MaP reads
    CUTADAPT_MAP (
        ch_map_reads_for_cutadapt
    )

    def ch_rtstop_trimmed_split = CUTADAPT_RTSTOP.out.reads.multiMap { meta, reads ->
        align: [ meta, reads ]
        fastqc: [ meta, reads ]
    }
    def ch_map_trimmed_split = CUTADAPT_MAP.out.reads.multiMap { meta, reads ->
        align: [ meta, reads ]
        fastqc: [ meta, reads ]
    }
    def ch_rtstop_trimmed_for_align = ch_rtstop_trimmed_split.align.map { meta, reads -> [ meta, reads ] }
    def ch_rtstop_trimmed_for_fastqc = ch_rtstop_trimmed_split.fastqc.map { meta, reads -> [ meta, reads ] }
    def ch_map_trimmed_for_align = ch_map_trimmed_split.align.map { meta, reads -> [ meta, reads ] }
    def ch_map_trimmed_for_fastqc = ch_map_trimmed_split.fastqc.map { meta, reads -> [ meta, reads ] }
    def ch_trimmed_reads = ch_rtstop_trimmed_for_fastqc.mix(ch_map_trimmed_for_fastqc)
    ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT_RTSTOP.out.log.collect { cutadapt_log -> cutadapt_log[1] })
    ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT_MAP.out.log.collect { cutadapt_log -> cutadapt_log[1] })
    def ch_cutadapt_adapter_mqc = CUTADAPT_RTSTOP.out.log
        .map { meta, cutadapt_log ->
            [
                meta.id.toString(),
                [
                    cutadapt_mode: 'RT-stop',
                    adapter_5p  : parseCutadaptCommandArg(cutadapt_log, '-g'),
                    adapter_3p  : parseCutadaptCommandArg(cutadapt_log, '-a')
                ]
            ]
        }
        .mix(
            CUTADAPT_MAP.out.log.map { meta, cutadapt_log ->
                [
                    meta.id.toString(),
                    [
                        cutadapt_mode: 'MaP',
                        adapter_5p  : parseCutadaptCommandArg(cutadapt_log, '-g'),
                        adapter_3p  : parseCutadaptCommandArg(cutadapt_log, '-a')
                    ]
                ]
            }
        )
        .collect()
        .map { rows -> cutadaptAdaptersMultiqc(rows) }
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_cutadapt_adapter_mqc.collectFile(
            name: 'cutadapt_adapters_mqc.yaml',
            sort: true
        )
    )

    // MODULE: fastqc (post-trim) — quality control on trimmed reads
    FASTQC_POST (
        ch_trimmed_reads.map { meta, reads -> [ meta + [id: "${meta.id}_trimmed"], reads ] }
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_POST.out.zip.collect { fastqc_zip -> fastqc_zip[1] })

    emit:
    reads_branched = ch_samplesheet_for_branching // channel: [ val(meta), [ reads ] ] (post-cat)
    rtstop_trimmed = ch_rtstop_trimmed_for_align   // channel: [ val(meta), [ reads ] ]
    map_trimmed    = ch_map_trimmed_for_align      // channel: [ val(meta), [ reads ] ]
    multiqc_files  = ch_multiqc_files              // channel: queue of mqc files
}
