// Subworkflow with functionality specific to the nf-core/rnastructurome pipeline

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    // Print version and exit if required and dump pipeline parameters to JSON file
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    // Validate parameters and generate parameter summary to stdout
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/rnastructurome ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/rnastructurome/blob/master/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        false
    )

    // Check config provided to the pipeline
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    // Create channel from input file provided through `input`

    def samplesheet_rows = samplesheetToList(input, "${projectDir}/assets/schema_input.json")
    warnOnSamplesheetOverrides(samplesheet_rows)

    channel
        .fromList(samplesheet_rows)
        .map {
            meta, fastq_1, fastq_2 ->
                def resolved_meta = meta + [
                    sample_id     : meta.sample_id ?: params.sample_id,
                    library_layout: fastq_2 ? 'PAIRED' : 'SINGLE',
                    method        : meta.method ?: params.method,
                    principle     : meta.principle ?: params.principle,
                    chemical      : meta.chemical ?: params.chemical,
                    RT_enzyme     : meta.RT_enzyme ?: params.RT_enzyme,
                    pH            : hasMetadataValue(meta.pH) ? meta.pH : params.pH,
                    organism      : meta.organism ?: params.organism,
                    adapter_3p    : meta.adapter_3p,
                    adapter_5p    : meta.adapter_5p,
                    umi_pattern   : meta.umi_pattern ?: params.umi_pattern,
                    sample_group  : meta.sample_group,
                    replicate     : meta.replicate
                ]

                if (!fastq_2) {
                    return [ resolved_meta.id, resolved_meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ resolved_meta.id, resolved_meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .map {
            meta, fastqs ->
                return [ meta, fastqs.flatten() ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    // Completion email and summary
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def normaliseMetadataValue(value) {
    return value.toString().trim().toLowerCase()
}

// The per-sample params are fallbacks for rows leaving the column blank — a row that sets it wins.
// Warn once per field so the precedence is visible rather than silent. Compared case-insensitively,
// as principle/method are read downstream: warning about a value the user did not actually
// contradict is worse than staying quiet.
def warnOnSamplesheetOverrides(rows) {
    ['sample_id', 'method', 'principle', 'chemical', 'RT_enzyme', 'pH', 'organism', 'umi_pattern'].each { field ->
        def param_value = params[field]
        if (!hasMetadataValue(param_value)) {
            return
        }
        def conflicting = rows
            .collect { row -> row[0] }
            .findAll { meta ->
                hasMetadataValue(meta[field]) &&
                    normaliseMetadataValue(meta[field]) != normaliseMetadataValue(param_value)
            }
        if (conflicting) {
            log.warn(
                "--${field}=${param_value} is overridden by the samplesheet for ${conflicting.size()} sample(s) " +
                "(e.g. '${conflicting[0].id}' = ${conflicting[0][field]}); the samplesheet takes precedence."
            )
        }
    }
}

def hasMetadataValue(value) {
    if (value == null) {
        return false
    }
    if (value instanceof Collection) {
        return !value.isEmpty()
    }
    return value.toString().trim()
}

// Validate channels from input samplesheet
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
// Tools credited in the MultiQC methods section, in pipeline order. Citation and reference live in one
// entry so the two lists cannot drift; `when` mirrors the conditions the workflow actually branches on.
def pipelineToolReferences(transcriptome) {
    def counts_on_genome = !transcriptome && params.count_genome
    [
        [ when: true,
          cite: 'FastQC (Andrews 2010)',
          ref : '<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/.</li>' ],
        [ when: params.umi_pattern as Boolean,
          cite: 'UMI-tools (Smith et al. 2017)',
          ref : '<li>Smith T, et al. (2017). UMI-tools: modelling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. Genome Research, 27(3), 491–499. doi: 10.1101/gr.209601.116</li>' ],
        [ when: true,
          cite: 'Cutadapt (Martin 2011)',
          ref : '<li>Martin M (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10–12. doi: 10.14806/ej.17.1.200</li>' ],
        [ when: !transcriptome && !params.count_genome,
          cite: 'GffRead (Pertea & Pertea 2020)',
          ref : '<li>Pertea G & Pertea M (2020). GFF Utilities: GffRead and GffCompare. F1000Research, 9, 304. doi: 10.12688/f1000research.23297.2</li>' ],
        [ when: !transcriptome,
          cite: 'STAR (Dobin et al. 2013)',
          ref : '<li>Dobin A, et al. (2013). STAR: ultrafast universal RNA-seq aligner. Bioinformatics, 29(1), 15–21. doi: 10.1093/bioinformatics/bts635</li>' ],
        [ when: transcriptome as Boolean,
          cite: 'Bowtie (Langmead et al. 2009)',
          ref : '<li>Langmead B, et al. (2009). Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biology, 10(3), R25. doi: 10.1186/gb-2009-10-3-r25</li>' ],
        [ when: transcriptome as Boolean,
          cite: 'Bowtie2 (Langmead & Salzberg 2012)',
          ref : '<li>Langmead B & Salzberg SL (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9(4), 357–359. doi: 10.1038/nmeth.1923</li>' ],
        [ when: true,
          cite: 'SAMtools (Danecek et al. 2021)',
          ref : '<li>Danecek P, et al. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. doi: 10.1093/gigascience/giab008</li>' ],
        [ when: counts_on_genome,
          cite: 'BEDOPS (Neph et al. 2012)',
          ref : '<li>Neph S, et al. (2012). BEDOPS: high-performance genomic feature operations. Bioinformatics, 28(14), 1919–1920. doi: 10.1093/bioinformatics/bts277</li>' ],
        [ when: counts_on_genome,
          cite: 'RSeQC (Wang et al. 2012)',
          ref : '<li>Wang L, Wang S & Li W (2012). RSeQC: quality control of RNA-seq experiments. Bioinformatics, 28(16), 2184–2185. doi: 10.1093/bioinformatics/bts356</li>' ],
        [ when: true,
          cite: 'RNAFramework (Incarnato et al. 2018)',
          ref : '<li>Incarnato D, et al. (2018). RNA Framework: an all-in-one toolkit for the analysis of RNA structures and post-transcriptional modifications. Nucleic Acids Research, 46(W1), W121–W127. doi: 10.1093/nar/gky486</li>' ],
        [ when: !params.stop_after_jackknife,
          cite: 'ViennaRNA (Lorenz et al. 2011)',
          ref : '<li>Lorenz R, et al. (2011). ViennaRNA Package 2.0. Algorithms for Molecular Biology, 6, 26. doi: 10.1186/1748-7188-6-26</li>' ],
        [ when: !params.stop_after_jackknife && params.r2dt,
          cite: 'R2DT (Sweeney et al. 2021)',
          ref : '<li>Sweeney BA, et al. (2021). R2DT is a framework for predicting and visualising RNA secondary structure using templates. Nature Communications, 12(1), 3494. doi: 10.1038/s41467-021-23555-5</li>' ],
        [ when: !params.stop_after_jackknife,
          cite: 'UCSC wigToBigWig (Kent et al. 2010)',
          ref : '<li>Kent WJ, et al. (2010). BigWig and BigBed: enabling browsing of large distributed datasets. Bioinformatics, 26(17), 2204–2207. doi: 10.1093/bioinformatics/btq351</li>' ],
        [ when: true,
          cite: 'MultiQC (Ewels et al. 2016)',
          ref : '<li>Ewels P, et al. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047–3048. doi: 10.1093/bioinformatics/btw354</li>' ]
    ].findAll { tool -> tool.when }
}

def toolCitationText(transcriptome) {
    def tools = pipelineToolReferences(transcriptome).collect { tool -> tool.cite }
    return "Tools used in the workflow included: ${tools.join(', ')}."
}

def toolBibliographyText(transcriptome) {
    return pipelineToolReferences(transcriptome).collect { tool -> tool.ref }.join(' ')
}

def methodsDescriptionText(mqc_methods_yaml, transcriptome) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Loop to handle multiple DOIs, stripping `https://doi.org/` (DOIs vs resolvers) and spaces
        // (manifest.doi is a string, not a proper list).
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    meta["tool_citations"] = toolCitationText(transcriptome)
    meta["tool_bibliography"] = toolBibliographyText(transcriptome)


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PIPELINE HELPER FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Value transforms, file parsers, arg renderers and MultiQC table builders.
// No channel or process logic.

def normaliseEnsemblSpecies(value) {
    value
        ?.toString()
        ?.trim()
        ?.toLowerCase()
        ?.replaceAll(/[^a-z0-9_]+/, '_')
        ?.replaceAll(/^_+|_+$/, '')
}

def resolveReferenceKey(meta) {
    def sampleId = meta?.id ?: 'unknown'
    def rawReference = (meta?.organism ?: params.organism)?.toString()?.trim()
    if (!rawReference) {
        error("Missing organism for sample '${sampleId}'. Set organism in the samplesheet or provide --organism.")
    }
    if (rawReference.contains(' ')) {
        return normaliseEnsemblSpecies(rawReference)
    }
    if (rawReference ==~ /[a-z]+_[a-z0-9_]+/) {
        return rawReference.toLowerCase()
    }
    rawReference
}

// Resolve how one sample's reference artifact (kind = 'fasta' | 'gtf') is obtained, returning
// [ reference_key, "<scheme>::<value>", original_organism ] (scheme: path:: | ensembl:: | ncbi:: | none::).
def resolveReferenceResolution(meta, kind) {
    def reference_key     = resolveReferenceKey(meta)
    def original_organism = meta.organism?.toString() ?: reference_key
    def genome_entry      = params.genomes?.containsKey(reference_key) ? params.genomes[reference_key] : null

    // 1. Explicit local path (user-supplied or from params.genomes)
    if (kind == 'fasta') {
        def local_fasta = params.fasta
        if (local_fasta) {
            return [ reference_key, "path::${local_fasta.toString()}", original_organism ]
        }
        def transcript_fasta = genome_entry?.transcript_fasta ?: genome_entry?.transcriptome ?: genome_entry?.cdna
        if (transcript_fasta) {
            return [ reference_key, "path::${transcript_fasta.toString()}", original_organism ]
        }
    } else {
        if (params.gtf) {
            return [ reference_key, "path::${params.gtf.toString()}", original_organism ]
        }
        def gtf_path = genome_entry?.gtf
        if (gtf_path) {
            return [ reference_key, "path::${gtf_path.toString()}", original_organism ]
        }
    }

    // 2. Explicit Ensembl species (per-genome override or species map)
    def explicit_ensembl = genome_entry?.ensembl_species ?: params.ensembl_species_map?.get(reference_key)
    if (explicit_ensembl) {
        return [ reference_key, "ensembl::${explicit_ensembl.toLowerCase()}", original_organism ]
    }

    // 3. Pre-configured NCBI accessions. For FASTA this skips Ensembl (avoids spurious 404
    // log lines); for GTF the annotation is the synthetic NCBI_GTF, signalled by none::.
    def ncbi_accessions = genome_entry?.ncbi_accessions ?: params.ncbi_accessions_map?.get(reference_key)
    if (ncbi_accessions) {
        if (kind == 'fasta') {
            def acc_str = (ncbi_accessions instanceof List) ? ncbi_accessions.join(',') : ncbi_accessions.toString()
            return [ reference_key, "ncbi::${acc_str}", original_organism ]
        }
        return [ reference_key, "none::", original_organism ]
    }

    // 4. Default: Ensembl first. A 404 there triggers the NCBI fallback downstream.
    if (kind == 'fasta' && !reference_key) {
        error("No organism specified for sample '${meta.id}'. Provide --fasta, --organism, or set params.genomes.")
    }
    return [ reference_key, "ensembl::${reference_key}", original_organism ]
}

// Decide whether a run should default to the transcriptome (Bowtie) route: NCBI (bacteria/viral)
// references have no introns, so STAR offers nothing — auto-enable when every reference resolves to
// NCBI (samplesheets are single-organism-class by contract). Returns false on any parse error so the
// normal validation path in PIPELINE_INITIALISATION still runs and reports.
def allReferencesUseNcbiRoute(samplesheetPath, schemaPath) {
    if (!samplesheetPath) return false
    try {
        def rows = samplesheetToList(samplesheetPath.toString(), schemaPath.toString())
        if (!rows) return false
        def organisms = rows.collect { row ->
            def meta = (row instanceof List) ? row[0] : row
            (meta?.organism ?: params.organism)?.toString()?.trim()
        }
        if (organisms.any { org -> !org }) return false
        return organisms.every { org ->
            resolveReferenceResolution([ id: 'route-probe', organism: org ], 'fasta')[1].startsWith('ncbi::')
        }
    } catch (Exception _ignored) {
        return false
    }
}

// Collect a queue channel of [key, value] pairs into a single value channel holding
// a [key: value] map (last write wins on duplicate keys). Empty input yields [:].
def collectToMap(ch_keyed) {
    // .collect() on an empty channel emits nothing; ifEmpty keeps the result a real
    // (reusable value) map so downstream .combine() isn't silently emptied (e.g. no --gtf).
    ch_keyed
        .map { key, value -> [ (key): value ] }
        .collect()
        .map { entries -> entries.inject([:]) { acc, entry -> acc + entry } }
        .ifEmpty([:])
}

// Build STAR_ALIGN inputs for a set of trimmed reads: pair each sample with its reference's STAR index +
// GTF (combine by:0) and emit [ [meta,reads], [idx_meta,index], [gtf_meta,gtf], ignore_gtf ].
// ignore_gtf skips --sjdbGTFfile at align time only when no GTF was resolved for that reference.
def buildStarAlignInputs(ch_trimmed, ch_star_index, ch_gtf) {
    def ch_keyed = ch_trimmed
        .map { meta, reads -> [ resolveReferenceKey(meta), meta, reads ] }
    def ch_idx_keyed = ch_star_index
        .map { meta, index -> [ meta.id.toString(), meta, index ] }
    def ch_gtf_keyed = ch_gtf
        .map { meta, gtf -> [ meta.id.toString(), meta, gtf ] }
    def ch_star_ref = ch_idx_keyed
        .join(ch_gtf_keyed, remainder: true)
        .map { combined ->
            def ref_key  = combined[0]
            def idx_meta = combined[1]
            def index    = combined[2]
            def gtf_meta = combined.size() > 3 ? combined[3] : null
            def gtf      = combined.size() > 4 ? combined[4] : null
            [ ref_key, idx_meta, index, gtf_meta ?: [id: 'no_gtf'], gtf ?: [], gtf_meta != null && gtf != null ]
        }
    ch_keyed
        .combine(ch_star_ref, by: 0)
        .map { _ref_key, sample_meta, reads, idx_meta, index, gtf_meta, gtf, has_gtf ->
            [ [sample_meta, reads], [idx_meta, index], [gtf_meta, gtf], !has_gtf ]
        }
}

// Collapse per-sample reference resolutions to one resolution per reference_key,
// erroring if a single reference resolved inconsistently across its samples.
def uniqueReferenceResolution(ch_requests, label) {
    ch_requests
        .groupTuple()
        .map { reference_key, resolutions, organisms ->
            def unique_resolutions = resolutions.unique()
            if (unique_resolutions.size() != 1) {
                error("Multiple ${label} resolutions were detected for reference '${reference_key}': ${unique_resolutions.join(', ')}")
            }
            [ reference_key, unique_resolutions[0], organisms[0] ]
        }
}

def parseFlagstatMappedReads(flagstatFile) {
    def mappedLine = flagstatFile.readLines().find { line ->
        line ==~ /^\d+\s+\+\s+\d+\s+mapped\s+\(.*/
    }
    if (!mappedLine) {
        error("Could not parse mapped read count from flagstat file: ${flagstatFile}")
    }
    (mappedLine.tokenize()[0]) as long
}

def parseRfcountCoveredTranscripts(summaryFile) {
    def summaryLines = summaryFile.readLines().findAll { line -> line?.trim() }
    if (summaryLines.size() < 2) {
        error("Could not parse rf-count summary TSV: ${summaryFile}")
    }
    def fields = summaryLines[1].split('\t')
    if (fields.size() < 2) {
        error("rf-count summary TSV is missing the covered transcript column: ${summaryFile}")
    }
    (fields[1]) as long
}

def parseInferExperiment(txtFile) {
    def forward = 0.0
    def reverse = 0.0
    txtFile.readLines().each { line ->
        def m = line =~ /Fraction of reads explained by "(?:1\+\+,1--,2\+-,2-\+|\+\+,--)": (.+)/
        if (m) forward = m[0][1].trim() as double
        m = line =~ /Fraction of reads explained by "(?:1\+-,1-\+,2\+\+,2--|\+-,-\+)": (.+)/
        if (m) reverse = m[0][1].trim() as double
    }
    // rf-count-genome's "second-strand" assigns transcript strand from read1's own mapped orientation
    // (read1=sense); "first-strand" uses read2's orientation (read2=sense, e.g. dUTP/TruSeq-directional).
    // RSeQC's "forward" fraction ("1++,1--,2+-,2-+") is read1=sense, so it maps to 'second' here; its
    // "reverse" fraction ("1+-,1-+,2++,2--", the common dUTP pattern) is read2=sense, mapping to 'first'.
    if (forward > 0.7) return 'second'
    if (reverse > 0.7) return 'first'
    return 'unstranded'
}

def filterSummaryParams(summaryParams, transcriptome) {
    def hiddenKeys = [
        'ensembl_species_map',
        'ncbi_accessions_map',
        'genomes',
        'container',
        'configFiles',
        'launchDir',
        'projectDir',
        'userName',
        'workDir'
    ] as Set

    // Hide aligner-specific params for the route not in use.
    // Transcriptome route uses Bowtie (RT-stop) / Bowtie2 (MaP); genome route uses STAR.
    if (transcriptome) {
        hiddenKeys.addAll(['star_map_sjdb_overhang', 'star_multimap_nmax'])
    } else {
        hiddenKeys.addAll(summaryParams.values()
            .findAll { section -> section instanceof Map }
            .collectMany { section -> section.keySet() as List }
            .findAll { key -> key.startsWith('bowtie_') || key.startsWith('bowtie2_') })
    }

    // Hide genome-route-specific params for transcriptome runs, and vice versa
    if (transcriptome) {
        hiddenKeys.addAll(summaryParams.values()
            .findAll { section -> section instanceof Map }
            .collectMany { section -> section.keySet() as List }
            .findAll { key -> key.startsWith('rfcount_genome_') })
    } else {
        hiddenKeys.addAll(summaryParams.values()
            .findAll { section -> section instanceof Map }
            .collectMany { section -> section.keySet() as List }
            .findAll { key -> key.startsWith('rfcount_map_') })
    }

    // Hide jackknife params when jackknife is not configured
    if (!params.jackknife_reference) {
        hiddenKeys.addAll(summaryParams.values()
            .findAll { section -> section instanceof Map }
            .collectMany { section -> section.keySet() as List }
            .findAll { key -> key.startsWith('rfjackknife_') || key == 'jackknife_reference' })
    }

    // Hide rfeval params when rfeval is not configured
    if (!params.rfeval_reference) {
        hiddenKeys.addAll(summaryParams.values()
            .findAll { section -> section instanceof Map }
            .collectMany { section -> section.keySet() as List }
            .findAll { key -> key.startsWith('rfeval_') || key == 'rfeval_reference' })
    }

    summaryParams.collectEntries { sectionName, sectionParams ->
        if (!(sectionParams instanceof Map)) {
            return [(sectionName): sectionParams]
        }

        def filteredSection = sectionParams.findAll { key, _value ->
            if (hiddenKeys.contains(key)) {
                return false
            }
            true
        }

        [(sectionName): filteredSection]
    }.findAll { _sectionName, sectionParams ->
        !(sectionParams instanceof Map) || !sectionParams.isEmpty()
    }
}

def addModuleOptionsSummary(summaryParams, transcriptome) {
    def sampleMetadata = parseInputSamplesheetMetadata(params.input)
    def moduleOptions = buildModuleOptionsSummary(sampleMetadata, transcriptome)
    if (moduleOptions.isEmpty()) {
        return summaryParams
    }
    summaryParams + ['Module options': moduleOptions]
}

def parseInputSamplesheetMetadata(inputPath) {
    if (!inputPath) {
        return [principles: [], conditions: [], methods: [], adapter_5p: [], adapter_3p: []]
    }

    def inputFile = file(inputPath.toString())
    if (!inputFile.exists()) {
        return [principles: [], conditions: [], methods: [], adapter_5p: [], adapter_3p: []]
    }

    def lines = inputFile.readLines().findAll { line -> line?.trim() }
    if (lines.size() < 2) {
        return [principles: [], conditions: [], methods: [], adapter_5p: [], adapter_3p: []]
    }

    def header = lines[0].split(',', -1)*.trim()
    def principleIdx = header.indexOf('principle')
    def conditionIdx = header.indexOf('condition')
    def methodIdx = header.indexOf('method')
    def adapter5pIdx = header.indexOf('adapter_5p')
    def adapter3pIdx = header.indexOf('adapter_3p')

    def principles = []
    def conditions = []
    def methods = []
    def adapter5p = []
    def adapter3p = []

    lines.drop(1).each { line ->
        def fields = line.split(',', -1)
        if (principleIdx >= 0 && principleIdx < fields.size()) {
            def value = fields[principleIdx]?.trim()
            if (value) principles << value
        }
        if (conditionIdx >= 0 && conditionIdx < fields.size()) {
            def value = fields[conditionIdx]?.trim()
            if (value) conditions << value
        }
        if (methodIdx >= 0 && methodIdx < fields.size()) {
            def value = fields[methodIdx]?.trim()
            if (value) methods << value
        }
        if (adapter5pIdx >= 0 && adapter5pIdx < fields.size()) {
            def value = fields[adapter5pIdx]?.trim()
            if (value) adapter5p << value
        }
        if (adapter3pIdx >= 0 && adapter3pIdx < fields.size()) {
            def value = fields[adapter3pIdx]?.trim()
            if (value) adapter3p << value
        }
    }

    [
        principles: principles.unique(),
        conditions: conditions.collect { condition -> condition.toLowerCase() }.unique(),
        methods   : methods.unique(),
        adapter_5p: adapter5p.unique(),
        adapter_3p: adapter3p.unique()
    ]
}

def buildModuleOptionsSummary(sampleMetadata, transcriptome) {
    def moduleOptions = [:]
    def principles = (sampleMetadata.principles ?: []).collect { principle -> principle.toLowerCase() }
    def adapter5p = (sampleMetadata.adapter_5p ?: []).findAll { adapter -> adapter?.trim() }
    def adapter3p = (sampleMetadata.adapter_3p ?: []).findAll { adapter -> adapter?.trim() }

    if (!adapter5p.isEmpty()) {
        moduleOptions['cutadapt_adapter_5p'] = adapter5p.join(', ')
    }
    if (!adapter3p.isEmpty()) {
        moduleOptions['cutadapt_adapter_3p'] = adapter3p.join(', ')
    }

    if (transcriptome) {
        moduleOptions['aligner'] = 'transcriptome'
        if (!principles || principles.contains('rt-stop')) {
            moduleOptions['rtstop_aligner'] = 'bowtie'
            moduleOptions['rtstop_aligner_args'] = renderBowtie1Args()
        }
        if (principles.contains('map')) {
            moduleOptions['map_aligner'] = 'bowtie2'
            moduleOptions['map_aligner_args'] = renderBowtie2Args()
        }
    } else {
        moduleOptions['aligner'] = 'star'
    }

    def rfnormSummary = renderRfNormSummary(sampleMetadata)
    moduleOptions.putAll(rfnormSummary)

    moduleOptions.findAll { _k, v -> v != null && v.toString().trim() }
}

def renderBowtie1Args() {
    def args = []
    if (params.bowtie_all as Boolean) {
        args << '-a'
    } else if (params.bowtie_k != null) {
        args << "-k ${params.bowtie_k as Integer}"
    }
    if ((params.bowtie_trim5 as Integer) > 0) {
        args << "--trim5 ${params.bowtie_trim5 as Integer}"
    }
    if ((params.bowtie_trim3 as Integer) > 0) {
        args << "--trim3 ${params.bowtie_trim3 as Integer}"
    }
    args << '-l 28'
    if (params.bowtie_v != null) {
        args << "-v ${params.bowtie_v as Integer}"
    } else {
        args << "-n ${params.bowtie_n as Integer}"
    }
    if (!(params.bowtie_all as Boolean) && params.bowtie_k == null && params.bowtie_max != null) {
        args << "-m ${params.bowtie_max as Integer}"
        if ((params.bowtie_max as Integer) > 1) {
            args << '-a'
        }
    }
    args << '--best'
    args << '--strata'
    args << "--chunkmbs ${params.bowtie_chunkmbs as Integer}"
    args.join(' ').trim()
}

def renderBowtie2Args() {
    def args = []
    def preset = params.bowtie2_preset?.toString()?.trim() ?: ''
    if (params.bowtie_all as Boolean) {
        args << '-a'
    } else if (!preset && params.bowtie_k != null) {
        args << "-k ${params.bowtie_k as Integer}"
    }
    if ((params.bowtie_trim5 as Integer) > 0) {
        args << "--trim5 ${params.bowtie_trim5 as Integer}"
    }
    if ((params.bowtie_trim3 as Integer) > 0) {
        args << "--trim3 ${params.bowtie_trim3 as Integer}"
    }
    if (preset) {
        args << preset
        if (params.bowtie2_softclip as Boolean) {
            args << "--ma ${params.bowtie2_ma as Integer}"
        }
    } else {
        args << '-L 22'
        if (params.bowtie2_softclip as Boolean) {
            args << '--local'
            args << "--ma ${params.bowtie2_ma as Integer}"
        }
    }
    args << "--mp ${params.bowtie2_mp}"
    args << "--dpad ${params.bowtie2_dpad as Integer}"
    args << "--rdg ${params.bowtie2_rdg}"
    args << "--rfg ${params.bowtie2_rfg}"
    if (params.bowtie2_dovetail as Boolean) {
        args << '--dovetail'
    }
    args.join(' ').trim()
}

def renderRfNormSummary(sampleMetadata) {
    def principles = (sampleMetadata.principles ?: []).collect { principle -> principle.toLowerCase() }.unique()
    def conditions = (sampleMetadata.conditions ?: []).collect { condition -> condition.toLowerCase() }.unique()
    if (principles.size() != 1) {
        return [
            rfnorm_mode: 'dynamic (mixed principles across samples)'
        ]
    }

    def principle = principles[0]
    def hasUntreated = conditions.contains('untreated')
    def hasDenatured = conditions.contains('denatured')
    def scoringMethod = resolveRfNormScoreMethod(principle, hasUntreated)
    def normMethod = resolveRfNormNormMethod(scoringMethod)
    def isDmsOnly = ((sampleMetadata.methods ?: []).collect { method -> method.toLowerCase() }.unique()) == ['dms']
    def isDmsBroad = isDmsOnly && sampleMetadata.pH != null && (sampleMetadata.pH as Double) >= 8.0
    def reactiveBases = params.rfnorm_reactive_bases ?: (isDmsOnly ? (isDmsBroad ? 'ACGU' : 'AC') : null)
    def dynamicWindow = params.rfnorm_dynamic_window != null ? (params.rfnorm_dynamic_window as Integer) : (isDmsOnly && !isDmsBroad ? 50 : null)

    def args = [
        "-sm ${scoringMethod}",
        "-nm ${normMethod}"
    ]
    if (params.rfnorm_remap_reactivities as Boolean) args << '--remap-reactivities'
    if (reactiveBases) args << "--reactive-bases ${reactiveBases}"
    if (params.rfnorm_norm_window != null) args << "--norm-window ${params.rfnorm_norm_window as Integer}"
    if (params.rfnorm_window_offset != null) args << "--window-offset ${params.rfnorm_window_offset as Integer}"
    if (dynamicWindow != null) args << "--dynamic-window ${dynamicWindow}"
    if (params.rfnorm_norm_independent as Boolean) args << '--norm-independent'
    if (params.rfnorm_raw as Boolean) args << '--raw'
    if (params.rfnorm_pseudocount != null) args << "--pseudocount ${params.rfnorm_pseudocount}"
    if (params.rfnorm_ignore_lower_than_untreated as Boolean) args << '--ignore-lower-than-untreated'
    def meanCoverage = params.rfnorm_mean_coverage != null ? params.rfnorm_mean_coverage as BigDecimal : 0
    if (meanCoverage > 0) args << "--mean-coverage ${params.rfnorm_mean_coverage}"
    def medianCoverage = params.rfnorm_median_coverage != null ? params.rfnorm_median_coverage as BigDecimal : 0
    if (medianCoverage > 0) args << "--median-coverage ${params.rfnorm_median_coverage}"
    def nanThreshold = params.rfnorm_nan != null ? params.rfnorm_nan as Integer : 10
    if (nanThreshold != 10) args << "--nan ${nanThreshold}"
    args << '--img'
    args << "-R ${params.rnaframework_r_path}"

    [
        rfnorm_mode         : "${principle.toUpperCase()} ${hasUntreated ? 'with untreated' : 'treated-only'}${hasDenatured ? ' + denatured' : ''}",
        rfnorm_scoring      : "${rfNormScoringLabel(scoringMethod)} (sm=${scoringMethod})",
        rfnorm_normalisation: "${rfNormNormLabel(normMethod)} (nm=${normMethod})",
        rfnorm_args         : args.join(' ').trim()
    ]
}

// Returns the base sample_group identifier (portion before the first underscore), e.g.
// "MDA-MB-231_MTX" → "MDA-MB-231". Used by fuzzy untreated-pairing to match a shared root.
def sampleGroupBaseToken(String sample_group) {
    sample_group.tokenize('_')[0]
}

def resolveRfNormScoreMethod(principle, hasUntreated) {
    def defaultMethod = principle == 'map' ? (hasUntreated ? 3 : 4) : (hasUntreated ? 1 : 2)
    if (params.rfnorm_score_method == null) {
        return defaultMethod
    }

    def requestedMethod = params.rfnorm_score_method as Integer
    if (!(requestedMethod in [1, 2, 3, 4])) {
        error("Unsupported rf-norm scoring method '${params.rfnorm_score_method}'. Expected one of: 1, 2, 3, 4.")
    }

    requestedMethod
}

def resolveRfNormNormMethod(scoringMethod) {
    def defaultMethod = (scoringMethod as Integer) == 2 ? 2 : 3
    if (params.rfnorm_norm_method == null) {
        return defaultMethod
    }

    def requestedMethod = params.rfnorm_norm_method as Integer
    if (!(requestedMethod in [2, 3, 4])) {
        error("Unsupported rf-norm normalization method '${params.rfnorm_norm_method}'. Expected one of: 2, 3, 4.")
    }

    requestedMethod
}

def rfNormScoringLabel(code) {
    def labels = [
        1: 'Ding',
        2: 'Rouskin',
        3: 'Siegfried',
        4: 'Zubradt'
    ]
    labels[code as Integer] ?: 'unknown'
}

def rfNormNormLabel(code) {
    def labels = [
        2: '90% Winsorizing',
        3: 'Box-plot normalisation',
        4: 'Mitchell normalisation'
    ]
    labels[code as Integer] ?: 'unknown'
}

def parseCutadaptCommandArg(logFile, optionName) {
    def commandLine = logFile.readLines().find { line -> line.startsWith('Command line parameters:') }
    if (!commandLine) {
        return 'none'
    }
    def pattern = java.util.regex.Pattern.compile("(?:^|\\s)${java.util.regex.Pattern.quote(optionName)}\\s+(\\S+)")
    def matcher = pattern.matcher(commandLine)
    matcher.find() ? matcher.group(1) : 'none'
}

// Normalise the assorted shapes Nextflow collect() can produce (single, nested, flattened, or Map) into
// a uniform list of [id, map] entries. With a non-null `label`, an unrecognised shape errors; else [].
def normaliseMqcRows(rows, label = null) {
    if (rows instanceof Map) {
        return rows.entrySet().collect { entry -> [entry.key, entry.value] }
    }
    if (rows instanceof List && rows.size() == 2 && rows[1] instanceof Map && !(rows[0] instanceof List)) {
        return [rows]
    }
    if (rows instanceof List && rows.every { row -> row instanceof List && row.size() == 2 && row[1] instanceof Map }) {
        return rows
    }
    if (rows instanceof List && rows.size() % 2 == 0 && rows.collate(2).every { pair -> pair.size() == 2 && pair[1] instanceof Map }) {
        return rows.collate(2)
    }
    if (label) {
        error("Unexpected ${label} row structure: ${rows?.getClass()?.name} -> ${rows}")
    }
    []
}

def countProgressionMultiqc(rows) {
    def orderedRows = normaliseMqcRows(rows, 'count progression').sort { a, b -> a[0] <=> b[0] }
    def dataBlock = orderedRows.collect { row ->
        def sample_id = row[0]
        def metrics = row[1]
        def metricLines = metrics.collect { key, value ->
            def rendered = value instanceof BigDecimal ? String.format(java.util.Locale.ROOT, '%.2f', value) : value.toString()
            "    ${key}: ${rendered}"
        }.join('\n')
        "  ${sample_id}:\n${metricLines}"
    }.join('\n')

    """id: 'nf-core-rnastructurome-count-progression'
section_name: 'nf-core/rnastructurome Count Progression'
description: 'Read counts at each stage of the alignment → deduplication → RF-count pipeline per sample.'
plot_type: 'table'
pconfig:
  id: 'nf-core-rnastructurome-count-progression'
  title: 'nf-core/rnastructurome Count Progression'
  show_table_by_default: true
headers:
  mapped_reads_pre_dedup:
    title: 'Mapped (pre-dedup)'
    description: 'Reads mapped to the reference before deduplication'
    scale: 'Blues'
    format: '{:,.0f}'
  mapped_reads_post_dedup:
    title: 'Mapped (post-dedup)'
    description: 'Reads retained after UMI/positional deduplication'
    scale: 'Blues'
    format: '{:,.0f}'
  pct_removed_by_dedup:
    title: 'Removed by Dedup'
    description: 'Percentage of mapped reads removed as duplicates'
    scale: 'Oranges'
    format: '{:,.1f}'
    suffix: '%'
  rfcount_covered_transcripts:
    title: 'RF-count: Covered Transcripts'
    description: 'Number of transcripts with sufficient coverage in RF-count'
    scale: 'Greens'
    format: '{:,.0f}'
data:
${dataBlock}
"""
}

def parseRfnormLog(logFile) {
    def covered   = 0L
    def discarded = 0L
    logFile.readLines().each { line ->
        def covM = (line =~ /\[\*\]\s+Covered transcripts:\s+(\d+)/)
        if (covM.find()) covered = covM.group(1) as long
        def disM = (line =~ /\[\*\]\s+Discarded transcripts:\s+(\d+)\s+total/)
        if (disM.find()) discarded = disM.group(1) as long
    }
    [covered: covered, discarded: discarded]
}

def parseRffoldLog(logFile) {
    def folded    = 0L
    def discarded = 0L
    logFile.readLines().each { line ->
        def foldM = (line =~ /\[\*\]\s+Folded transcripts:\s+(\d+)/)
        if (foldM.find()) folded = foldM.group(1) as long
        def disM  = (line =~ /\[\*\]\s+Discarded transcripts:\s+(\d+)\s+total/)
        if (disM.find()) discarded = disM.group(1) as long
    }
    [folded: folded, discarded: discarded]
}

def buildSimpleMultiqcTable(rows, id, sectionName, description, headers) {
    def orderedRows = normaliseMqcRows(rows).sort { a, b -> a[0] <=> b[0] }
    def dataBlock = orderedRows.collect { row ->
        def sampleId = row[0]
        def metrics  = row[1]
        def metricLines = metrics.collect { key, value -> "    ${key}: ${value}" }.join('\n')
        "  ${sampleId}:\n${metricLines}"
    }.join('\n')
    def headerBlock = headers.collect { col, cfg ->
        def lines = ["  ${col}:"]
        cfg.each { k, v -> lines << "    ${k}: '${v}'" }
        lines.join('\n')
    }.join('\n')

    """id: '${id}'
section_name: '${sectionName}'
description: '${description}'
plot_type: 'table'
pconfig:
  id: '${id}'
  title: '${sectionName}'
  show_table_by_default: true
headers:
${headerBlock}
data:
${dataBlock}
"""
}

def rfnormStatsMultiqc(rows) {
    buildSimpleMultiqcTable(
        rows,
        'nf-core-rnastructurome-rfnorm-stats',
        'nf-core/rnastructurome RF-norm Statistics',
        'Transcript coverage statistics from RF-count and RF-norm (per normalisation group).',
        [
            rfcount_covered: [title: 'RF-count Covered', description: 'Transcripts covered by RF-count (input to RF-norm)', scale: 'Blues',  format: '{:,.0f}'],
            covered        : [title: 'RF-norm Good Coverage',  description: 'Transcripts passing RF-norm normalisation (sufficient coverage)', scale: 'Greens', format: '{:,.0f}']
        ]
    )
}

def rffoldStatsMultiqc(rows) {
    buildSimpleMultiqcTable(
        rows,
        'nf-core-rnastructurome-rffold-stats',
        'nf-core/rnastructurome RF-fold Statistics',
        'Folded transcripts per fold group, and transcripts dropped from the consensus fold because they were not covered in every replicate.',
        [
            folded   : [title: 'Folded Transcripts',   description: 'Transcripts successfully folded by rf-fold', scale: 'Purples', format: '{:,.0f}'],
            discarded: [title: 'Discarded Transcripts', description: 'Transcripts normalised in a replicate but excluded from the consensus fold (not present in all replicates)', scale: 'Reds', format: '{:,.0f}']
        ]
    )
}

// Summarise one group's rf-eval metrics TSV for the MultiQC table. The per-structure rows live in the
// published TSV; here we reduce to a group-level view: how many structures were scored, and the median
// observed score for each metric against the median of its rotation baseline.
def parseRfevalMetrics(metricsFile) {
    def lines = metricsFile.readLines().findAll { line -> line.trim() }
    if (lines.size() < 2) {
        return [ structures: 0, median_coeff: 0, median_coeff_baseline: 0, median_dsci: 0, median_dsci_baseline: 0, median_auroc: 0, median_auroc_baseline: 0 ]
    }
    def header = lines[0].split('\t').toList()
    def columns = [
        coeff                : header.indexOf('coeff_unpaired'),
        coeff_baseline       : header.indexOf('coeff_unpaired_baseline_mean'),
        dsci                 : header.indexOf('dsci'),
        dsci_baseline        : header.indexOf('dsci_baseline_mean'),
        auroc                : header.indexOf('auroc'),
        auroc_baseline       : header.indexOf('auroc_baseline_mean')
    ]
    def values = [:]
    columns.each { name, _index -> values[name] = [] }

    def structures = 0
    lines.drop(1).each { line ->
        def parts = line.split('\t').toList()
        // rf-eval's own "Overall" row would double-count against the per-structure rows
        if (parts[0] == 'Overall') {
            return
        }
        structures += 1
        columns.each { name, index ->
            if (index >= 0 && index < parts.size() && parts[index]) {
                // rf-eval writes "nan" for structures with no usable reactivity; drop those from the medians
                def value = parts[index].isDouble() ? (parts[index] as Double) : Double.NaN
                if (!value.isNaN()) {
                    values[name] << value
                }
            }
        }
    }
    [
        structures           : structures,
        median_coeff         : roundedMedian(values.coeff),
        median_coeff_baseline: roundedMedian(values.coeff_baseline),
        median_dsci          : roundedMedian(values.dsci),
        median_dsci_baseline : roundedMedian(values.dsci_baseline),
        median_auroc         : roundedMedian(values.auroc),
        median_auroc_baseline: roundedMedian(values.auroc_baseline)
    ]
}

def roundedMedian(values) {
    if (!values) {
        return 0
    }
    def sorted = values.sort(false)
    def mid = sorted.size().intdiv(2)
    def value = sorted.size() % 2 == 1 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2
    Math.round(value * 1000) / 1000.0
}

def rfevalStatsMultiqc(rows) {
    def normalised = normaliseMqcRows(rows)
    if (!normalised) {
        return ''
    }
    buildSimpleMultiqcTable(
        normalised,
        'nf-core-rnastructurome-rfeval',
        'nf-core/rnastructurome RF-eval Structure Agreement',
        'Agreement between measured reactivities and the reference structures given by --rfeval_reference, ' +
        'per normalisation group. Per-structure scores are in the published <code>eval/*_rfeval.metrics.tsv</code>. ' +
        'Each score is paired with its baseline: the same structure scored against rotations of its own ' +
        'reactivity profile, which is what it would score by chance. Read the gap, not the score — only AUROC ' +
        'has a fixed chance level, so the coefficient and DSCI baselines differ from one structure to the next.',
        [
            structures           : [title: 'Structures',        description: 'Reference structures scored in this group', scale: 'Blues', format: '{:,.0f}'],
            median_coeff         : [title: 'Median Coefficient', description: 'Median per-structure unpaired coefficient — the fraction of highly-reactive bases that are unpaired in the reference', scale: 'RdYlGn', min: 0, max: 1, format: '{:,.3f}'],
            median_coeff_baseline: [title: 'Coefficient Baseline', description: 'Median unpaired coefficient across rotation decoys; the raw score saturates near 1, so the gap to this baseline is what carries the signal', scale: 'Greys', min: 0, max: 1, format: '{:,.3f}'],
            median_dsci          : [title: 'Median DSCI',       description: 'Median per-structure DSCI', scale: 'RdYlGn', min: 0, max: 1, format: '{:,.3f}'],
            median_dsci_baseline : [title: 'DSCI Baseline',     description: 'Median DSCI across rotation decoys — what these structures score by chance', scale: 'Greys', min: 0, max: 1, format: '{:,.3f}'],
            median_auroc         : [title: 'Median AUROC',      description: 'Median per-structure AUROC', scale: 'RdYlGn', min: 0, max: 1, format: '{:,.3f}'],
            median_auroc_baseline: [title: 'AUROC Baseline',    description: 'Median AUROC across rotation decoys; sits at ~0.5 by construction, so it doubles as a sanity check', scale: 'Greys', min: 0, max: 1, format: '{:,.3f}']
        ]
    )
}

// Parse an rf-correlate matrix.csv into the off-diagonal summary for the MultiQC reproducibility table
// (replicate count, mean/min pairwise correlation). Header: Sample,<label0>,...; rows: <label_i>,<corr_i0>,...
def parseRfcorrelateMatrix(matrixFile) {
    def lines = matrixFile.readLines().findAll { line -> line.trim() }
    if (lines.size() < 2) {
        return [ replicates: 0, mean_corr: 0, min_corr: 0 ]
    }
    def labels = lines[0].split(',').drop(1)
    def values = []
    lines.drop(1).eachWithIndex { line, i ->
        def parts = line.split(',')
        // Upper-triangle off-diagonal only: column j > row i. parts[0] is the row label.
        ((i + 1)..<labels.size()).each { j ->
            def raw = (j + 1) < parts.size() ? parts[j + 1].toString().trim() : ''
            if (raw && raw.toLowerCase() != 'nan') {
                values << (raw as Double)
            }
        }
    }
    if (!values) {
        return [ replicates: labels.size(), mean_corr: 0, min_corr: 0 ]
    }
    def mean = values.sum() / values.size()
    [
        replicates: labels.size(),
        mean_corr : (Math.round(mean * 1000) / 1000.0),
        min_corr  : (Math.round(values.min() * 1000) / 1000.0)
    ]
}

// rows: per sample_group [id, [replicates, mean_pearson, mean_spearman]], already merged across
// both correlation methods in the CORRELATE_REPLICATES subworkflow (flattened by .collect()).
def rfCorrelateMultiqc(rows) {
    def normalised = normaliseMqcRows(rows)
    if (!normalised) {
        return ''
    }
    buildSimpleMultiqcTable(
        normalised,
        'nf-core-rnastructurome-rfcorrelate',
        'nf-core/rnastructurome Replicate Correlation',
        'Pairwise reactivity-profile correlation between replicates from rf-correlate (per sample group). Higher is more reproducible.',
        [
            replicates    : [title: 'Replicates',      description: 'Number of replicates compared', scale: 'Blues',  format: '{:,.0f}'],
            mean_pearson  : [title: 'Mean Pearson',     description: 'Mean pairwise replicate Pearson correlation (reactivity-capped, overall, transcriptome-wide)', scale: 'RdYlGn', min: 0, max: 1, format: '{:,.3f}'],
            mean_spearman : [title: 'Mean Spearman',    description: 'Mean pairwise replicate Spearman correlation (overall, transcriptome-wide)', scale: 'RdYlGn', min: 0, max: 1, format: '{:,.3f}']
        ]
    )
}

def cutadaptAdaptersMultiqc(rows) {
    def rowEntries
    if (rows instanceof Map) {
        rowEntries = rows.entrySet().collect { entry -> [entry.key, entry.value] }
    } else if (rows instanceof List && rows.size() == 2 && rows[1] instanceof Map && !(rows[0] instanceof List)) {
        rowEntries = [rows]
    } else if (rows instanceof List) {
        rowEntries = rows.collectMany { row ->
            if (row instanceof Map) {
                return row.entrySet().collect { entry -> [entry.key, entry.value] }
            }
            if (row instanceof List && row.size() == 2 && row[1] instanceof Map) {
                return [row]
            }
            return []
        }
    } else {
        rowEntries = []
    }
    def orderedRows = rowEntries
        .findResults { row ->
            if (row instanceof List && row.size() == 2 && row[1] instanceof Map) {
                return [row[0].toString(), row[1]]
            }
            null
        }
        .sort { a, b -> a[0] <=> b[0] }
    def dataBlock = orderedRows.collect { row ->
        def sample_id = row[0]
        def metrics = row[1]
        def metricLines = metrics.collect { key, value ->
            def rendered = value.toString().replace("'", "''")
            "    ${key}: '${rendered}'"
        }.join('\n')
        "  ${sample_id}:\n${metricLines}"
    }.join('\n')

    """id: 'nf-core-rnastructurome-cutadapt-adapters'
section_name: 'Cutadapt: Adapter Sequences Used'
description: 'Adapter sequences used for trimming in each sample.'
parent_id: 'cutadapt'
parent_name: 'Cutadapt'
plot_type: 'table'
pconfig:
  id: 'nf-core-rnastructurome-cutadapt-adapters'
  title: 'Cutadapt: Adapter Sequences Used'
  show_table_by_default: true
headers:
  cutadapt_mode:
    title: 'Trim Mode'
  adapter_5p:
    title: \"5' Adapter\"
  adapter_3p:
    title: \"3' Adapter\"
data:
${dataBlock}
"""
}
