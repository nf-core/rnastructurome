#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/rnastructurome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/rnastructurome
    Website: https://nf-co.re/rnastructurome
    Slack  : https://nfcore.slack.com/channels/rnastructurome
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RNASTRUCTUROME  } from './workflows/rnastructurome'
include { PIPELINE_INITIALISATION   } from './subworkflows/local/utils_nfcore_rnastructurome_pipeline'
include { PIPELINE_COMPLETION       } from './subworkflows/local/utils_nfcore_rnastructurome_pipeline'
include { allReferencesUseNcbiRoute } from './subworkflows/local/utils_nfcore_rnastructurome_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// WORKFLOW: Run main analysis pipeline depending on type of input
workflow NFCORE_RNASTRUCTUROME {

    take:
    samplesheet      // channel: samplesheet read in from --input
    transcriptome    // boolean: transcriptome (Bowtie) route, including auto-detection

    main:

    // WORKFLOW: Run pipeline
    RNASTRUCTUROME (
        samplesheet,
        transcriptome
    )
    emit:
    multiqc_report  = RNASTRUCTUROME.out.multiqc_report  // channel: /path/to/multiqc_report.html
    normalized_xml  = RNASTRUCTUROME.out.normalized_xml  // channel: [ val(meta), path(xml) ]
    fold_structures = RNASTRUCTUROME.out.fold_structures // channel: [ val(meta), path(dir) ]
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    // The only pipeline setting not readable straight from params: the transcriptome route can be
    // auto-detected below, and params is immutable (assigning to it is silently ignored). Everything
    // downstream reads params directly; this one boolean is threaded through instead.
    def transcriptome = params.transcriptome as Boolean

    // NCBI (bacteria/viral) references have no introns, so STAR offers nothing over Bowtie — auto-enable
    // the transcriptome route when every reference resolves to NCBI, so users don't need --transcriptome.
    if (!transcriptome
            && allReferencesUseNcbiRoute(params.input, "${projectDir}/assets/schema_input.json")) {
        log.info('[rnastructurome] All references resolve to the NCBI route (no-intron organisms) — enabling the transcriptome (Bowtie) route automatically. Pass --transcriptome to set it explicitly.')
        transcriptome = true
    }

    // A user-supplied --fasta without --gtf can only be used on the transcriptome route: the genome route
    // needs the GTF to extract and count transcripts. Assume the FASTA is a transcriptome and enable the
    // route automatically rather than failing deep in STAR — pass --gtf to use the genome route instead.
    if (!transcriptome && params.fasta && !params.gtf) {
        log.warn('[rnastructurome] --fasta was supplied without --gtf; enabling the transcriptome route.')
        transcriptome = true
    }
    // SUBWORKFLOW: Run initialisation tasks
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    // WORKFLOW: Run main workflow
    NFCORE_RNASTRUCTUROME (
        PIPELINE_INITIALISATION.out.samplesheet,
        transcriptome
    )
    // SUBWORKFLOW: Run completion tasks
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_RNASTRUCTUROME.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
