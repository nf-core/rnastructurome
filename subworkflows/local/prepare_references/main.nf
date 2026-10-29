// PREPARE_REFERENCES — resolve each sample's reference to a local path, Ensembl species, or NCBI
// accessions (with Ensembl-404 fallback), fetch/sort FASTA+GTF into lookup maps, then build the aligner
// index (STAR or Bowtie/Bowtie2). Source-provenance channels are emitted for RDAT naming downstream.

include { ENSEMBL_TRANSCRIPTOME } from '../../../modules/local/ensembl/transcriptome/main'
include { ENSEMBL_GENOME        } from '../../../modules/local/ensembl/genome/main'
include { ENSEMBL_GTF           } from '../../../modules/local/ensembl/gtf/main'
include { NCBI_FASTA            } from '../../../modules/local/ncbi/fasta/main'
include { NCBI_GTF              } from '../../../modules/local/ncbi/gtf/main'
include { FASTA_SANITIZE as FASTA_SANITIZE_LOCAL   } from '../../../modules/local/fasta/sanitize/main'
include { FASTA_SANITIZE as FASTA_SANITIZE_ENSEMBL } from '../../../modules/local/fasta/sanitize/main'
include { FASTA_SANITIZE as FASTA_SANITIZE_NCBI    } from '../../../modules/local/fasta/sanitize/main'
include { GTF_SANITIZE          } from '../../../modules/local/gtf/sanitize/main'
include { FASTA_SORT as FASTA_SORT_LOCAL   } from '../../../modules/local/fasta/sort/main'
include { FASTA_SORT as FASTA_SORT_ENSEMBL } from '../../../modules/local/fasta/sort/main'
include { FASTA_SORT as FASTA_SORT_NCBI    } from '../../../modules/local/fasta/sort/main'
include { STAR_GENOMEGENERATE   } from '../../../modules/nf-core/star/genomegenerate/main'
include { BOWTIE_BUILD          } from '../../../modules/nf-core/bowtie/build/main'
include { BOWTIE2_BUILD         } from '../../../modules/nf-core/bowtie2/build/main'
include { SAMTOOLS_FAIDX        } from '../../../modules/nf-core/samtools/faidx/main'
include { GFFREAD               } from '../../../modules/nf-core/gffread/main'

include { resolveReferenceResolution } from '../utils_nfcore_rnastructurome_pipeline/main'
include { uniqueReferenceResolution  } from '../utils_nfcore_rnastructurome_pipeline/main'
include { collectToMap               } from '../utils_nfcore_rnastructurome_pipeline/main'
include { resolveReferenceKey        } from '../utils_nfcore_rnastructurome_pipeline/main'

// Organisms known to have transcript/gene IDs (e.g. parenthesised yeast tRNA IDs like tK(UUU)K,
// or yeast systematic isoform ids like YAL016C-A) that hang RNAframework's XML parser or need
// normalizing before FASTA_SORT — see GTF_SANITIZE and FASTA_SANITIZE below. Extend as new
// offenders are found; other organisms skip these processes entirely rather than pay for a
// no-op sanitize pass.
def idSanitizeOrganisms() {
    ['saccharomyces_cerevisiae']
}

workflow PREPARE_REFERENCES {

    take:
    ch_samplesheet_for_branching // channel: [ val(meta), [ reads ] ] (post-cat)
    transcriptome              // boolean: transcriptome (Bowtie) route, including auto-detection

    main:

    // Branch by probing principle (for transcriptome-route reference-FASTA selection)
    def principle_branches = ch_samplesheet_for_branching.branch { meta, _reads ->
        rtstop: (meta.principle ?: '').toLowerCase() == 'rt-stop'
        map:    (meta.principle ?: '').toLowerCase() == 'map'
    }

    def ch_reference_requests = uniqueReferenceResolution(
        ch_samplesheet_for_branching.map { meta, _reads -> resolveReferenceResolution(meta, 'fasta') },
        'transcript reference'
    )

    def ch_reference_local = ch_reference_requests
        .filter { _reference_key, resolution, _org -> resolution.startsWith('path::') }
        .map { reference_key, resolution, _original_organism ->
            def fasta_path = resolution - 'path::'
            [ [ id: reference_key, organism: reference_key ], file(fasta_path, checkIfExists: true) ]
        }

    def ch_reference_ensembl_input = ch_reference_requests
        .filter { _reference_key, resolution, _org -> resolution.startsWith('ensembl::') }
        .map { reference_key, resolution, original_organism ->
            def ensembl_species = resolution - 'ensembl::'
            [ [ id: reference_key, organism: reference_key, ensembl_species: ensembl_species,
                original_organism: original_organism ], ensembl_species ]
        }

    // Ensembl downloads are mutually exclusive: genome mode gets the soft-masked genome FASTA (STAR),
    // transcriptome mode gets the cDNA FASTA (Bowtie). ch_ensembl_not_found triggers NCBI fallback.
    def ch_ensembl_not_found        = channel.empty()
    def ch_ensembl_fasta_source_url = channel.empty()
    def ch_reference_genome_fasta_keyed = channel.empty()

    if (transcriptome) {
        ENSEMBL_TRANSCRIPTOME (
            ch_reference_ensembl_input,
            [
                ensembl_release : params.ensembl_release,
                ensembl_base_url: params.ensembl_base_url
            ]
        )
        ch_ensembl_not_found        = ENSEMBL_TRANSCRIPTOME.out.not_found
        ch_ensembl_fasta_source_url = ENSEMBL_TRANSCRIPTOME.out.source_urls
    }

    // MODULE: ENSEMBL_GENOME — download soft-masked genome FASTA for STAR alignment.
    // ch_reference_genome_fasta_keyed is populated after sorting (see FASTA_SORT_ENSEMBL block below).
    if (!transcriptome) {
        ENSEMBL_GENOME(
            ch_reference_ensembl_input,
            [
                ensembl_release : params.ensembl_release,
                ensembl_base_url: params.ensembl_base_url
            ]
        )
        ch_ensembl_not_found        = ENSEMBL_GENOME.out.not_found
        ch_ensembl_fasta_source_url = ENSEMBL_GENOME.out.source_url
    }

    // Organisms not found on Ensembl FTP are routed to NCBI_FASTA for automatic accession search;
    // meta.original_organism carries the esearch query when no accessions are pre-configured.
    def ch_reference_ncbi_from_ensembl = ch_ensembl_not_found
        .map { meta, _not_found_file -> [ meta, "" ] }

    def ch_reference_ncbi_explicit = ch_reference_requests
        .filter { _reference_key, resolution, _org -> resolution.startsWith('ncbi::') }
        .map { reference_key, resolution, original_organism ->
            def accessions = resolution - 'ncbi::'
            [ [ id: reference_key, organism: reference_key, original_organism: original_organism ], accessions ]
        }

    def ch_reference_ncbi_input = ch_reference_ncbi_explicit.mix(ch_reference_ncbi_from_ensembl)

    NCBI_FASTA (
        ch_reference_ncbi_input
    )

    NCBI_GTF (
        NCBI_FASTA.out.fasta
    )

    // GTF resolution is skipped entirely when stop_after_jackknife + transcriptome, since every
    // GTF-consuming step is gated behind one of those, avoiding wasted Ensembl lookups.
    def ch_all_reference_gtf        = channel.empty()
    def ch_reference_gtf_local      = channel.empty()
    def ch_ensembl_gtf_source_urls  = channel.empty()

    if (!(params.stop_after_jackknife && transcriptome)) {
        def ch_reference_gtf_requests = uniqueReferenceResolution(
            ch_samplesheet_for_branching.map { meta, _reads -> resolveReferenceResolution(meta, 'gtf') },
            'GTF reference'
        )

        // Local GTF (already in GTF format) — used directly
        ch_reference_gtf_local = ch_reference_gtf_requests
            .filter { _k, resolution, _o -> resolution.startsWith('path::') }
            .map { reference_key, resolution, _original_organism ->
                [ [ id: reference_key, organism: reference_key ], file(resolution - 'path::', checkIfExists: true) ]
            }

        def ch_reference_gtf_ensembl_input = ch_reference_gtf_requests
            .filter { _k, resolution, _o -> resolution.startsWith('ensembl::') }
            .map { reference_key, resolution, original_organism ->
                def ensembl_species = resolution - 'ensembl::'
                [ [ id: reference_key, organism: reference_key, ensembl_species: ensembl_species,
                    original_organism: original_organism ], ensembl_species ]
            }

        ENSEMBL_GTF (
            ch_reference_gtf_ensembl_input,
            [
                ensembl_release : params.ensembl_release,
                ensembl_base_url: params.ensembl_base_url
            ]
        )
        ch_ensembl_gtf_source_urls = ENSEMBL_GTF.out.source_urls

        // NCBI references get their annotation from NCBI_GTF's synthetic GTF (already run above);
        // ENSEMBL_GTF.out.not_found is silently dropped since those organisms have a GTF via NCBI_GTF.
        ch_all_reference_gtf = ch_reference_gtf_local
            .mix(ENSEMBL_GTF.out.gtf)
            .mix(NCBI_GTF.out.gtf)
    }

    // Strip regex/shell-unsafe characters from transcript_id/gene_id (e.g. yeast tRNA IDs like
    // tK(UUU)K) that hang RNAframework's XML parser; single source feeding ch_reference_gtf_map.
    // Only organisms in idSanitizeOrganisms() actually need this — skip it for everything else.
    def gtf_sanitize_branches = ch_all_reference_gtf.branch { meta, _gtf ->
        needs_sanitize: meta.organism in idSanitizeOrganisms()
        clean:          true
    }

    GTF_SANITIZE (
        gtf_sanitize_branches.needs_sanitize
    )
    ch_all_reference_gtf = GTF_SANITIZE.out.gtf.mix(gtf_sanitize_branches.clean)

    // Same organism-gated normalization as GTF_SANITIZE above, applied to FASTA header ids
    // (yeast isoform hyphens, unsafe characters) before FASTA_SORT — one branch per source since
    // each keeps its own aliased FASTA_SORT/downstream channel. FASTA_SORT always decompresses
    // regardless of input, so the "clean" branch needs no separate gzip handling.
    def ch_reference_local_branches = ch_reference_local.branch { meta, _fasta ->
        needs_sanitize: meta.organism in idSanitizeOrganisms()
        clean:          true
    }
    FASTA_SANITIZE_LOCAL (
        ch_reference_local_branches.needs_sanitize
    )
    FASTA_SORT_LOCAL (
        FASTA_SANITIZE_LOCAL.out.fasta.mix(ch_reference_local_branches.clean)
    )

    def ch_fasta_sort_ensembl_out = channel.empty()
    if (transcriptome) {
        def ch_ensembl_fasta_branches = ENSEMBL_TRANSCRIPTOME.out.fasta.branch { meta, _fasta ->
            needs_sanitize: meta.organism in idSanitizeOrganisms()
            clean:          true
        }
        FASTA_SANITIZE_ENSEMBL (
            ch_ensembl_fasta_branches.needs_sanitize
        )
        FASTA_SORT_ENSEMBL (
            FASTA_SANITIZE_ENSEMBL.out.fasta.mix(ch_ensembl_fasta_branches.clean)
        )
        ch_fasta_sort_ensembl_out = FASTA_SORT_ENSEMBL.out.fasta
    } else {
        // Genome route: sort the Ensembl genome FASTA so it is published to reference/
        // and so STAR_GENOMEGENERATE receives a chromosome-sorted FASTA.
        def ch_ensembl_fasta_branches = ENSEMBL_GENOME.out.fasta.branch { meta, _fasta ->
            needs_sanitize: meta.organism in idSanitizeOrganisms()
            clean:          true
        }
        FASTA_SANITIZE_ENSEMBL (
            ch_ensembl_fasta_branches.needs_sanitize
        )
        FASTA_SORT_ENSEMBL (
            FASTA_SANITIZE_ENSEMBL.out.fasta.mix(ch_ensembl_fasta_branches.clean)
        )
        ch_fasta_sort_ensembl_out = FASTA_SORT_ENSEMBL.out.fasta
    }

    def ch_reference_ncbi_fasta_branches = NCBI_FASTA.out.fasta.branch { meta, _fasta ->
        needs_sanitize: meta.organism in idSanitizeOrganisms()
        clean:          true
    }
    FASTA_SANITIZE_NCBI (
        ch_reference_ncbi_fasta_branches.needs_sanitize
    )
    FASTA_SORT_NCBI (
        FASTA_SANITIZE_NCBI.out.fasta.mix(ch_reference_ncbi_fasta_branches.clean)
    )

    def ch_reference_fasta_keyed = FASTA_SORT_LOCAL.out.fasta
        .mix(ch_fasta_sort_ensembl_out)
        .mix(FASTA_SORT_NCBI.out.fasta)
        .map { meta, fasta -> [ meta.id.toString(), [meta, fasta] ] }
    def ch_reference_gtf_keyed   = ch_all_reference_gtf.map { meta, gtf -> [ meta.id.toString(), [meta, gtf] ] }
    def ch_reference_fasta_map   = collectToMap(ch_reference_fasta_keyed)
    def ch_reference_gtf_map     = collectToMap(ch_reference_gtf_keyed)

    // Genome FASTA map for STAR index building: Ensembl species use ENSEMBL_GENOME's soft-masked genome;
    // NCBI species (no introns) use the sorted NCBI FASTA (must be sorted/uncompressed — STAR/samtools
    // faidx reject the raw plain-gzip download); local FASTA references use the user-supplied FASTA as-is.
    def ch_reference_genome_ncbi_keyed = FASTA_SORT_NCBI.out.fasta
        .map { meta, fasta -> [ meta.id.toString(), [meta, fasta] ] }
    if (!transcriptome) {
        // Genome route: use the sorted Ensembl genome FASTA for STAR index building.
        ch_reference_genome_fasta_keyed = FASTA_SORT_ENSEMBL.out.fasta
            .map { meta, fasta -> [ meta.id.toString(), [meta, fasta] ] }
            .mix(ch_reference_genome_ncbi_keyed)
            .mix(FASTA_SORT_LOCAL.out.fasta.map { meta, fasta -> [ meta.id.toString(), [meta, fasta] ] })
    } else {
        ch_reference_genome_fasta_keyed = ch_reference_genome_fasta_keyed
            .mix(ch_reference_genome_ncbi_keyed)
    }

    // MODULE: gffread + samtools faidx — extract a spliced transcript FASTA from the genome FASTA
    // + GTF (genome route, default count_genome=false only). ALIGN_READS runs SAMTOOLS_CALMD against
    // this transcript FASTA on STAR's --quantMode TranscriptomeSAM output, so rf-count can run
    // directly on transcript-coordinate BAMs instead of rf-count-genome + rf-rctools extract.
    def ch_genome_transcript_fasta_fai_map = channel.value([:])
    def ch_genome_transcript_fasta_map     = channel.value([:])
    if (!transcriptome && !params.count_genome) {
        def ch_gffread_split = ch_reference_genome_fasta_keyed
            .join(ch_reference_gtf_keyed)
            .map { _key, fasta_entry, gtf_entry -> [ gtf_entry, fasta_entry[1] ] }
            .multiMap { gtf_entry, fasta ->
                gtf:   gtf_entry
                fasta: fasta
            }
        GFFREAD(ch_gffread_split.gtf, ch_gffread_split.fasta)

        SAMTOOLS_FAIDX(
            GFFREAD.out.gffread_fasta.map { meta, fasta -> [ meta, fasta, [] ] },
            false
        )

        def ch_genome_transcript_fasta_fai_keyed = GFFREAD.out.gffread_fasta
            .map { meta, fasta -> [ meta.id.toString(), meta, fasta ] }
            .join(SAMTOOLS_FAIDX.out.fai.map { meta, fai -> [ meta.id.toString(), fai ] })
            .map { key, meta, fasta, fai -> [ key, [ meta, fasta, fai ] ] }
        ch_genome_transcript_fasta_fai_map = collectToMap(ch_genome_transcript_fasta_fai_keyed)

        ch_genome_transcript_fasta_map = collectToMap(
            GFFREAD.out.gffread_fasta.map { meta, fasta -> [ meta.id.toString(), [meta, fasta] ] }
        )
    }

    def ch_rtstop_reference_fasta = principle_branches.rtstop
        .combine(ch_reference_fasta_map)
        .map { combined ->
            def meta = combined[0]
            def ref_map = combined[2]
            def reference_key = resolveReferenceKey(meta)
            def ref_tuple = ref_map[reference_key]
            if (!ref_tuple) {
                error("No transcript FASTA resolved for reference '${reference_key}' in RT-stop branch.")
            }
            ref_tuple
        }
        .map { meta, fasta -> [ meta.id.toString(), [meta, fasta] ] }
        .groupTuple()
        .map { _reference_key, entries -> entries[0] }

    def ch_map_reference_fasta = principle_branches.map
        .combine(ch_reference_fasta_map)
        .map { combined ->
            def meta = combined[0]
            def ref_map = combined[2]
            def reference_key = resolveReferenceKey(meta)
            def ref_tuple = ref_map[reference_key]
            if (!ref_tuple) {
                error("No transcript FASTA resolved for reference '${reference_key}' in MaP branch.")
            }
            ref_tuple
        }
        .map { meta, fasta -> [ meta.id.toString(), [meta, fasta] ] }
        .groupTuple()
        .map { _reference_key, entries -> entries[0] }

    // INDEX BUILDING — conditional on chosen aligner per principle
    def ch_star_index        = channel.empty()
    def ch_bowtie_index_map  = channel.value([:])
    def ch_bowtie2_index_map = channel.value([:])

    if (!transcriptome) {
        // Build one STAR index per reference using the genome FASTA + GTF (Ensembl soft-masked
        // toplevel assembly, or the NCBI FASTA for bacteria/viruses).
        def ch_star_build = ch_reference_genome_fasta_keyed
            .join(ch_reference_gtf_keyed, remainder: true)
            .map { key, fasta_tuple, gtf_tuple ->
                def fasta_meta = fasta_tuple[0]
                def fasta      = fasta_tuple[1]
                def gtf_meta   = gtf_tuple ? gtf_tuple[0] : [id: "${key}_gtf"]
                def gtf        = gtf_tuple ? gtf_tuple[1] : []
                [ [fasta_meta, fasta], [gtf_meta, gtf] ]
            }

        def ch_star_build_split = ch_star_build.multiMap { fasta_entry, gtf_entry ->
            fasta: fasta_entry
            gtf:   gtf_entry
        }
        STAR_GENOMEGENERATE(ch_star_build_split.fasta, ch_star_build_split.gtf)
        ch_star_index = STAR_GENOMEGENERATE.out.index
    }

    if (transcriptome) {
        BOWTIE_BUILD(ch_rtstop_reference_fasta)
        ch_bowtie_index_map = collectToMap(
            BOWTIE_BUILD.out.index.map { meta, index -> [ meta.id.toString(), [meta, index] ] }
        )
    }

    if (transcriptome) {
        BOWTIE2_BUILD(ch_map_reference_fasta)
        ch_bowtie2_index_map = collectToMap(
            BOWTIE2_BUILD.out.index.map { meta, index -> [ meta.id.toString(), [meta, index] ] }
        )
    }

    emit:
    fasta_map               = ch_reference_fasta_map          // value: map ref_key -> [meta, fasta]
    gtf_map                 = ch_reference_gtf_map            // value: map ref_key -> [meta, gtf]
    all_gtf                 = ch_all_reference_gtf            // channel: [ val(meta), path(gtf) ]
    genome_fasta_keyed      = ch_reference_genome_fasta_keyed // channel: [ key, [meta, fasta] ]
    genome_transcript_fasta_map     = ch_genome_transcript_fasta_map     // value: map ref_key -> [meta, fasta] (genome route, count_genome=false)
    genome_transcript_fasta_fai_map = ch_genome_transcript_fasta_fai_map // value: map ref_key -> [meta, fasta, fai] (genome route, count_genome=false)
    star_index              = ch_star_index                  // channel: [ val(meta), path(index) ] (genome route)
    bowtie_index_map        = ch_bowtie_index_map            // value: map ref_key -> [meta, index]
    bowtie2_index_map       = ch_bowtie2_index_map           // value: map ref_key -> [meta, index]
    ensembl_fasta_source_url = ch_ensembl_fasta_source_url   // channel: [ val(meta), path(urls) ]
    gtf_local               = ch_reference_gtf_local         // channel: [ val(meta), path(gtf) ]
    local_fasta_sorted      = FASTA_SORT_LOCAL.out.fasta     // channel: [ val(meta), path(fasta) ] (--fasta route naming)
    ncbi_source_accessions  = NCBI_FASTA.out.source_accessions // channel: [ val(meta), path(acc) ]
    ensembl_gtf_source_urls = ch_ensembl_gtf_source_urls     // channel: [ val(meta), path(urls) ] — empty when stop_after_jackknife + transcriptome
}
