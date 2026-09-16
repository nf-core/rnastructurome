# nf-core/rnastructurome: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/rnastructurome/usage](https://nf-co.re/rnastructurome/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Pipeline parameters

Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files supplied with `-c` can be used for infrastructure settings such as resources, executors, containers, or module arguments, but should not be used for ordinary pipeline parameters.

## Samplesheet input

The easiest way to run this pipeline is to create a full samplesheet that contains all of the information about each sample in a comma-separated file, including the header row shown below and pass it as `--input '[path to samplesheet file]'`:

```csv title="full_samplesheet.csv"
sample,sample_id,fastq_1,fastq_2,method,principle,chemical,RT_enzyme,sample_group,condition,replicate,organism,pH,adapter_3p,adapter_5p,umi_pattern
HEK293T_treated_r1,GSM000001,/data/treated_r1.fastq.gz,,SHAPE,RT-stop,NAI,M-MLV,HEK293T,treated,1,Homo sapiens,7.5,,,
HEK293T_untreated_r1,GSM000002,/data/untreated_r1.fastq.gz,,SHAPE,RT-stop,NAI,M-MLV,HEK293T,untreated,1,Homo sapiens,7.5,,,
```

An [example samplesheet](../assets/samplesheet.csv) is provided.

However, you can also provide a more minimal version if for example you don't need to specify some of the options, like in the example above you could decide to not provide the columns with information about the adapters or umi pattern if you are happy to use the default options.
You can also provide a very minimal samplesheet with just the information required about each invidivual sample and pass the uniform values across all samples as parameters. For example:

```csv title="minimal_samplesheet.csv"
sample,fastq_1,sample_group,condition,replicate
HEK293T_treated_r1,/data/treated_r1.fastq.gz,HEK293T,treated,1
HEK293T_untreated_r1,/data/untreated_r1.fastq.gz,HEK293T,untreated,1
```

then pass essential but uniform options like this:

```bash
--input /path/to/samplesheet.csv --method SHAPE --principle RT-stop --organism "Homo sapiens"
```

### Column reference

| Column         | Required | Required in samplesheet | Description                                                                                                             |
| -------------- | -------- | ----------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| `sample`       | yes      | yes                     | Sample identifier. Re-use the same name to concatenate re-sequenced samples (multi-lane support).                       |
| `sample_id`    | no       | no                      | Unique identifier for an individual sequencing run.                                                                     |
| `fastq_1`      | yes      | yes                     | Read 1 FASTQ path (`.fastq.gz` / `.fq.gz`).                                                                             |
| `fastq_2`      | no       | no                      | Read 2 FASTQ path for paired-end data.                                                                                  |
| `sample_group` | yes      | yes                     | Group key used for control pairing in `rf-norm`.                                                                        |
| `condition`    | yes      | yes                     | One of `treated`, `untreated`, `denatured`.                                                                             |
| `replicate`    | yes      | yes                     | Replicate key used for control pairing in `rf-norm`.                                                                    |
| `method`       | yes      | no                      | Probing chemistry: `SHAPE` or `DMS`. Controls chemistry-specific defaults. Falls back to `--method`.                    |
| `principle`    | yes      | no                      | Readout principle: `rt-stop` or `map`. Controls alignment, rf-count, and rf-norm defaults. Falls back to `--principle`. |
| `chemical`     | no       | no                      | Probing reagent. Falls back to `--chemical`; if unset, rf-fold uses generic slope/intercept defaults.                   |
| `RT_enzyme`    | no       | no                      | Reverse transcriptase enzyme. Falls back to `--RT_enzyme`; if unset, rf-count MaP settings use M-MLV defaults.          |
| `organism`     | yes      | no                      | Used for automatic genome/transcriptome download from Ensembl/NCBI (e.g. `Homo sapiens`). Falls back to `--organism`.   |
| `pH`           | no       | no                      | DMS reaction pH. Falls back to `--pH`; `pH >= 8.0` sets reactive bases to `ACGU`, otherwise DMS defaults to `AC`.       |
| `adapter_3p`   | no       | no                      | 3′ adapter sequence passed to Cutadapt. Falls back to `--cutadapt_adapter_3p`, then `AGATCGGAAGAGC`.                    |
| `adapter_5p`   | no       | no                      | 5′ adapter sequence passed to Cutadapt. Falls back to `--cutadapt_adapter_5p`, then `AGATCGGAAGAGC`.                    |
| `umi_pattern`  | no       | no                      | UMI pattern passed to umi_tools. Falls back to `--umi_pattern`; if unset, UMI extraction is skipped.                    |

## Other considerations and parameters

### Using your own genome/transcriptome files

This pipeline expects an `organism` parameter to be passed either in the samplesheet or as a global parameter, which will trigger an automatic download of FASTA and GTF files from Ensembl. In the cases of bacteria and viruses where Ensembl does not have a reference genome/transcriptome, it will look this up in NCBI and download the files from there instead. However, you can pass your own files if preferred. For this you need to use the `--fasta` and `--gtf` flags. If these flags are used, the pipeline will prioritise them over the Ensembl download even when `organism` is passed. Example of usage:

```bash
--input /path/to/samplesheet.csv --fasta /path/to/genome_or_transcriptome.fa.gz --gtf /path/to/annotation.gtf.gz
```

The `--gtf` flag is optional. If you supply only `--fasta`, the pipeline enables the transcriptome route automatically with a warning, since the default genome route needs the GTF to extract and count transcripts. Without a GTF the run still produces all transcript-coordinate outputs (counts, reactivities, structures, transcript browser tracks and 2D diagrams) and skips only the annotation-dependent steps, each with a warning.

### Using genome versus transcriptome

By default the pipeline aligns to the genome with STAR and extracts transcript-level counts using the GTF. To align directly to the transcriptome instead, pass `--transcriptome` and the rest happens automatically. See the [Alignment routes](#alignment-routes) section for more detail on when each route is appropriate.

### Using optional modules

Three optional RNAframework modules are available (each documented in more detail under [More details for key steps](#more-details-for-key-steps)):

1. **rf-jackknife**: tunes folding parameters against reference RNA structures. Pass `--jackknife_reference` with a path to a `.db` file of known structures. When provided, rf-jackknife runs between rf-norm and rf-fold and calibrates the slope and intercept parameters passed to rf-fold. Add `--stop_after_jackknife` to run calibration only and skip rf-fold and downstream structure outputs.
2. **rf-eval**: evaluates the agreement between your reactivity data and a set of reference RNA structures, computing metrics such as AUROC and DSCI that can be used as quality control. Pass `--rfeval_reference` with a path to a `.db` file of known structures to enable it.
3. **rf-structextract**: pulls high-confidence structural elements out of the folded structures. Enable it with `--structextract true`.

## More details for key steps

### Adapter trimming

Trimming is done with Cutadapt. Defaults differ by principle because RT-stop experiments encode signal at the read 3′ end and disabling 5′ quality trimming preserves those positions.

| Principle | `--cutadapt_5quality` | `--cutadapt_3quality` |
| --------- | --------------------- | --------------------- |
| `RT-stop` | forced to `0`         | default `20`          |
| `MaP`     | default `20`          | default `20`          |

Adapter precedence: per-sample columns (`adapter_5p`, `adapter_3p`) → global flags (`--cutadapt_adapter_5p`, `--cutadapt_adapter_3p`) → fallback `AGATCGGAAGAGC`.

Set `--cutadapt_quality_only` to skip adapter trimming and perform quality/length filtering only.

### Deduplication and UMI handling (optional)

By default the pipeline does **not** deduplicate reads (`--skip_markdup true`). Position-based deduplication is not valid for chemical-probing data without UMIs; reads that start at the same coordinate are independent molecules, not PCR duplicates. Set `--skip_markdup false` to enable SAMtools markdup if you know position-based dedup is appropriate for your library.

UMI-tagged samples are deduplicated with `umi_tools`, enabled via `umi_pattern` (per-sample column or `--umi_pattern` if the pattern is the same for all samples).

- Patterns containing only `N/C/X` use `umi_tools --bc-pattern` mode directly.
- Patterns that contain IUPAC (e.g. containing `D`) are converted to `--extract-method=regex` mode automatically. Each IUPAC character is expanded to a regex character class (`D` → `[AGT]`, `N` → `[ACGT]`, etc.) and the pattern is anchored to the read start as a named capture group. For example, `DNNN` becomes `^(?P<umi_1>[AGT][ACGT][ACGT][ACGT])`.

On the genome (STAR) route, STAR's transcriptome-coordinate BAM is produced at alignment time, before dedup. When dedup removes reads (UMI samples, or `--skip_markdup false`), the pipeline reconciles the transcript BAM against the surviving read names so transcript-level counts match the deduplicated genome BAM. Samples that undergo no read removal (non-UMI with the default `--skip_markdup true`) skip this reconciliation entirely.

### Reference resolution

The `organism` value (per-sample column or `--organism`) is used to resolve and download the reference automatically. Latin binomials (`Homo sapiens`) are preferred, but a small set of shorthands are also accepted: `human`, `mouse`, `rat`, and `yeast`.

For organisms not in Ensembl, the pipeline falls back to NCBI using pre-configured RefSeq accessions. The following are already built-in: Influenza A, SARS-CoV-2, Dengue virus, Zika virus, HIV-1, Rotavirus A, and Escherichia coli (K-12 MG1655). For any other organism, the pipeline will attempt an NCBI auto-search.

The pipeline always downloads the most recent available assembly - the current Ensembl release for eukaryotes, or the pinned RefSeq accession for NCBI organisms. To keep the reference version fixed across runs, run the pipeline once and then reuse the local files published under `reference/` by passing them with `--fasta` and `--gtf`.

Resolution order for each file type:

| File  | 1st       | 2nd              | 3rd                  |
| ----- | --------- | ---------------- | -------------------- |
| FASTA | `--fasta` | Ensembl download | NCBI fallback        |
| GTF   | `--gtf`   | Ensembl download | Synthetic GTF (NCBI) |

### Alignment routes

**Default: genome alignment with STAR**

The pipeline downloads a soft-masked genome FASTA (`dna_sm.toplevel.fa.gz`) and GTF from Ensembl, builds a STAR index, and aligns with splice-junction awareness. By default, STAR also emits transcript-coordinate alignments with `--quantMode TranscriptomeSAM`; these BAM files are reconciled with the deduplicated genome BAMs, corrected with `samtools calmd`, and counted directly with `rf-count` before passing to `rf-norm`. This is the recommended route for most experiments: STAR handles spliced reads correctly, the soft-masked genome reduces spurious multi-mappers from repetitive elements, and counting in transcript coordinates avoids known `rf-count-genome` edge cases at exon/intron boundaries.

If you prefer to use the genome-coordinate counting path, set `--count_genome true`. This runs `rf-count-genome` on the genome BAM and then converts the genome-coordinate RC files to transcript-level RC files with `rf-rctools extract` using the GTF.

`--star_multimap_nmax` (default `10`) sets the maximum number of genome locations a read is allowed to map to; reads exceeding this are discarded. The default of 10 is appropriate for most protein-coding genes. If you are probing highly repetitive RNAs such as rRNA or snRNA you may want to increase this (e.g. `--star_multimap_nmax 50`), though be aware that allowing too many multimappers can introduce noise if reads are assigned ambiguously.

The STAR preset also differs by readout principle. RT-stop uses `--alignEndsType Extend5pOfRead1`, `--outFilterMismatchNoverReadLmax 0.04`, and `--outSAMattributes NH HI AS NM MD` to preserve stop positions while retaining the tags needed downstream. MaP uses `--alignEndsType Local`, `--outFilterMismatchNoverReadLmax 0.15`, `--outSAMprimaryFlag AllBestScore`, `--outSAMattributes All`, and `--sjdbOverhang` (default `200`) for mutation-tolerant alignment. For additional available STAR flags check out the parameters page.

**Optional: transcriptome alignment with Bowtie**

Add `--transcriptome` to align directly to a transcript FASTA instead of the genome. In this option, the pipeline downloads Ensembl cDNA + ncRNA FASTA files or for bacteria and viruses not in Ensembl it falls back to NCBI and builds a transcript FASTA from the genome assembly automatically. The aligner used depends on the probing principle: Bowtie for RT-stop and Bowtie2 for MaP.

This route works well for bacteria and viruses since they lack introns and often their genome annotation is sparse or absent, or when you prefer to map directly to a curated set of transcripts. For bacteria and viruses, the pipeline enables the transcriptome (Bowtie) route automatically, so you do not need to pass `--transcriptome` explicitly.

The transcriptome aligners are configured to retain isoform-compatible mappings by default, but Bowtie v1 and Bowtie2 handle this slightly differently. For RT-stop, Bowtie v1 uses `-a --best --strata`, so it reports all alignments in the best alignment stratum rather than lower-scoring valid hits. For MaP, Bowtie2 uses `-a` with the `--very-sensitive-local` preset by default, reporting all alignments it finds under that local, mutation-tolerant search mode.

Set `--bowtie_all false` if you want more restrictive reporting. For RT-stop/Bowtie v1, `--bowtie_k 1` reports at most one best-stratum alignment per read; increasing `--bowtie_k` reports up to that many. For MaP/Bowtie2, setting `--bowtie_all false` with the default `--bowtie2_preset '--very-sensitive-local'` uses Bowtie2's single-alignment reporting. To apply `--bowtie_k` to Bowtie2 or control lower-level Bowtie2 search settings, first clear or replace `--bowtie2_preset` (for example, set it to an empty string in a params file). Other exposed Bowtie/Bowtie2 parameters are listed on the parameter page.

### rf-count

`rf-count` quantifies chemical probing signal from the aligned transcript-coordinate BAM files produced either by the transcriptome route or by STAR `--quantMode TranscriptomeSAM` on the default genome route. If `--count_genome true` is set, the genome route instead uses `rf-count-genome` followed by `rf-rctools extract`. What exactly is counted depends on the principle: for RT-stop experiments it tallies read 3′ ends that accumulate at modified bases; for MaP it measures per-position mutation rates.

When `--count_genome true` is used, strandedness is handled automatically for `rf-count-genome`: RT-stop libraries are always treated as second-strand (this is fixed by experimental design), while for MaP the pipeline infers strandedness per sample using `RSeQC infer_experiment`. If inference is ambiguous you can override it with `--rfcount_strandedness first|second|unstranded`. After genome-coordinate counting, `rf-rctools extract` converts the genome-coordinate RC files to transcript-level RC files using the GTF before passing to `rf-norm`.

For MaP samples, `RT_enzyme` controls several mutation-cleanup defaults. Set it per sample with the `RT_enzyme` samplesheet column, or globally with `--RT_enzyme` when the same enzyme was used for all rows. The pipeline currently distinguishes Group II Intron RTs, including TGIRT, from all other or unset enzymes; non-Group-II values use M-MLV-like defaults.

With M-MLV-like defaults, consecutive mutations/indels are collapsed (`--rfcount_map_collapse_consecutive true`, with `--rfcount_map_max_collapse_distance 2`) and deletion calls are assigned to the right-most deleted base (`--rfcount_map_right_deletion true`). With Group II Intron RTs, consecutive events are not collapsed by default; instead, insertion calls and ambiguously aligned deletions are ignored (`--rfcount_map_no_insertions true`, `--rfcount_map_no_ambiguous true`) and nearby mutations are filtered with `--rfcount_map_discard_consecutive 3`. Each of these automatic choices can be overridden explicitly with its corresponding `--rfcount_map_*` parameter.

Other parameters worth knowing about:

`--rfcount_img` generates per-transcript coverage plots. Disabled by default (`--rfcount_img false`) because they are slow to produce at transcriptome scale; enable with `--rfcount_img true` if you want them.

`--rfcount_trim_5prime` trims a fixed number of bases from the 5′ end of each read before counting. This is useful for RT-stop experiments where the first few bases after the adapter can carry sequence-context bias that inflates apparent stop rates.

`--rfcount_mask_file` accepts a BED file of regions to exclude from counting entirely; handy for masking rRNA or other highly-expressed contaminating transcripts that would otherwise dominate the output.

`--rfcount_primary_only` restricts counting to alignments marked as primary. This is most relevant when the aligner emits multiple records per read, such as STAR multimappers or Bowtie2 `--bowtie_all`/`-k` output; otherwise a read can contribute to more than one locus or transcript. For Bowtie v1 `--bowtie_all`, use `--bowtie_all false` and `--bowtie_k` if you need stricter one-alignment-per-read reporting.

For the full list of available options see the [rf-count documentation](https://rnaframework-docs.readthedocs.io/en/latest/rf-count/). Any flag not exposed as a pipeline parameter can be passed directly via `ext.args` in a custom config.

### rf-norm

`rf-norm` normalises per-position counts into reactivity scores. Samples are grouped by `sample_group + replicate`, and treated samples are normalised against their matched controls within the same group. The scoring method is selected automatically based on the probing principle and which conditions are present. The normalisation method defaults to box-plot normalisation in most cases.

You can override the scoring method with `--rfnorm_score_method`:

- `1` is Ding scoring (RT-stop, treated-vs-untreated)
- `2` is Rouskin scoring (RT-stop without an untreated control)
- `3` is Siegfried scoring (MaP, treated-vs-untreated)
- `4` is Zubradt scoring (MaP without an untreated control)

Similarly, you can override the normalisation method with `--rfnorm_norm_method`:

- `1` is 2-8% normalisation (takes the top 10% of reactivities, discards the very highest 2%, and uses the mean of the remaining 8% as the scaling factor)
- `2` is 90% Winsorizing (clips any value above the 95th or below the 5th percentile to those boundaries, then divides all values by the 95th percentile). By default samples scored with Rouskin method will use this normalisation method as it is the only one compatible with it.
- `3` is box-plot normalisation (removes outliers beyond 1.5× IQR then divides by the mean of the next top 10%), which is the default for most conditions
- `4` is Mitchell normalisation (MaP only, uses the higher of the mean 90th-95th percentile reactivity or the 75th percentile of non-zero reactivities as the scaling factor)

Treated samples are usually paired with untreated by matching `sample_group + replicate` exactly. However, in cases where the authors did not create an exact matching untreated sample for some specific treatments but have one for other samples in the same dataset, the pipeline falls back to an untreated sample sharing the same `sample_group` base token and replicate; for example, `MDA-MB-231_DMSO_treated_r1` will match exactly to `MDA-MB-231_DMSO_untreated_r1`, but `MDA-MB-231_MTX_treated_r1` will also pair with `MDA-MB-231_DMSO_untreated_r1` if no exact matching untreated sample exists. A warning is emitted when a fallback is used; if more than one candidate matches, the pipeline errors.

If a treated group still has no untreated after that (e.g. a dataset with only one untreated control shared across several differently-named treated replicates) and exactly one untreated control exists anywhere else for the same reference, the pipeline reuses that single control for it, again with a warning. This only fires when the choice is unambiguous: if two or more distinct untreated controls exist for the reference and some treated groups remain unmatched, the pipeline leaves them without a control (or errors downstream if `rf-normfactor` requires one; see below).

Disable both fallbacks with `--fuzzy_untreated_pairing false`, in which case unmatched groups proceed without a negative control.

For DMS experiments, a few defaults change automatically. `--rfnorm_reactive_bases` is set to `AC` (or `ACGU` when `pH ≥ 8`). `--rfnorm_dynamic_window` turns on when `pH < 8`, and `--rfnorm_norm_window` defaults to `50` for DMS or RT-stop samples; it controls the size of the sliding window used to compute local normalisation factors along the transcript. `--rfnorm_nan` (minimum reads at a position before reactivity is reported as NaN) defaults to `1000` for MaP or `50` for RT-stop. These can all be overridden explicitly if needed.

Other parameters worth knowing about:

`--rfnorm_img` generates per-transcript reactivity plots. Disabled by default (`--rfnorm_img false`) because rendering them via R for thousands of transcripts is the main source of rf-norm runtime on large datasets; enable with `--rfnorm_img true` if you want them.

`--rfnorm_mean_coverage` and `--rfnorm_median_coverage` discard transcripts whose mean or median coverage falls below the given threshold (both default `0`, i.e. no filtering).

`--rfnorm_prefilter_min_coverage` (genome route only) shrinks the full-annotation RC before rf-norm by keeping only transcripts with at least one covered base (default `1`; set `0` to keep the full annotation). This is what keeps the genome-route RC transcriptome-sized.

For the full list of available options see the [rf-norm documentation](https://rnaframework-docs.readthedocs.io/en/latest/rf-norm/). Any flag not exposed as a pipeline parameter can be passed directly via `ext.args` in a custom config.

### rf-normfactor (cross-experiment normalisation)

By default `rf-norm` normalises each `sample_group + replicate` group on its own, so reactivities are not directly comparable across samples. To fix that, whenever a reference has more than one treated sample (replicates or different conditions) the pipeline runs `rf-normfactor` automatically. It computes one set of transcriptome-wide normalisation factors from all of that reference's samples together, then applies them to every group (via `-nf`) so the whole reference sits on a single, comparable scale. Force this on or off with `--rfnorm_use_normfactor true|false`.

For the full list of available options see the [rf-normfactor documentation](https://rnaframework-docs.readthedocs.io/en/latest/rf-normfactor/).

### rf-correlate (replicate reproducibility)

When a sample group has more than one replicate, `rf-correlate` measures how reproducible those replicates are by computing pairwise correlations between their reactivity profiles, both transcriptome-wide and per-transcript. It works from the normalised reactivity XMLs produced by `rf-norm`, so the replicates are already on a comparable scale and the correlation reflects genuine reactivity agreement rather than differences in sequencing depth. It runs by default (and is a no-op for single-replicate groups); disable it with `--correlate_replicates false`.

The per-group correlations are summarised in a MultiQC table (one row per sample group, giving the replicate count and the mean pairwise Pearson and Spearman correlation) where a low value flags a discordant replicate before it dilutes the folded consensus. Both methods are computed for every comparison: **Pearson** on reactivities capped at `--rfcorrelate_cap_react` (default `1.5`, so a few extreme positions cannot dominate the score) and **Spearman** on ranks, to which the cap does not apply.

Set `--correlate_img true` to additionally write the rf-correlate correlation heatmap PDF (off by default).

For the full list of available options see the [rf-correlate documentation](https://rnaframework-docs.readthedocs.io/en/latest/rf-correlate/).

### rf-fold

`rf-fold` predicts RNA secondary structures from normalised reactivity profiles using ViennaRNA. The pipeline groups XMLs by `sample_group`, merging replicates, and folds them together. Dot-bracket output is the default; enable CT format with `--rffold_ct`.

`--rffold_slope` and `--rffold_intercept` control how reactivity values are converted into folding constraints. Left unset, they are chosen automatically from the sample's `chemical`: DMS uses slope `4.6`/intercept `-2` ([Borovská et al. 2026](https://doi.org/10.1038/s41587-025-02739-0)); NAI uses `2.2`/`-0.8` and 2A3 uses `1`/`-0.4` ([Marinus et al. 2021](https://doi.org/10.1093/nar/gkaa1255)); any other or unspecified chemical falls back to `1.8`/`-0.6`. Set either flag explicitly to override, or calibrate optimal values from your own data with `rf-jackknife` (see below); jackknife values take priority over both the flags and the chemical-based defaults.

For long transcripts, rf-fold can fold in a sliding window rather than the whole sequence at once. `--rffold_window` (default `1000`) sets the MFE folding window size and enables windowed folding; unset it (`null`) to fold the whole transcript. `--rffold_partition_window` (default `1000`) sets the partition-function window used for the Shannon entropy and dot-plot calculation, independently of MFE windowing, and `--rffold_vienna_max_bp_span` (default `600`) caps the maximum distance ViennaRNA allows between paired bases.

`--rffold_shannon_entropy` (default `true`) computes per-position Shannon entropy alongside the predicted structure, giving a measure of folding confidence. `--rffold_only_common` keeps only transcripts covered in at least N XML experiments; when not set explicitly, the pipeline enables it automatically for fold groups with more than one replicate, setting N to the group's replicate count so only transcripts present in every replicate are retained.

Once structures are folded, the pipeline draws reactivity-coloured 2D structure diagrams using [R2DT](https://github.com/RNAcentral/R2DT) where a template exists (rRNA, snRNA, tRNA, and other well-characterised RNA families). All folded RNAs will also be folded with ViennaRNA's `RNAplot`. R2DT is only available under a container profile (`docker`, `singularity`, `apptainer`); under `conda` or `mamba` all diagrams fall back to `RNAplot`.

The partition-function dot-plot (`--rffold_dotplot`, default `true`) records base-pairing probabilities across the whole structural ensemble. The pipeline converts these into base-pair (`.bp`) arc files and, together with the reactivity and Shannon-entropy profiles, exports them as genome-browser tracks; see [Genome-browser tracks](#genome-browser-tracks) for more information on this process. Each sample group's folded structures are also compiled with their rf-norm reactivities into RMDB-compatible **RDAT** files, recording the sequence, dot-bracket structure, reactivity values, and the parameters used to produce them.

For the full list of available options see the [rf-fold documentation](https://rnaframework-docs.readthedocs.io/en/latest/rf-fold/). Any flag not exposed as a pipeline parameter can be passed directly via `ext.args` in a custom config.

### Genome-browser tracks

For visualisation in IGV or the UCSC browser, the pipeline exports three kinds of track: per-base **reactivity**, folding-confidence (**Shannon entropy**), and base-pairing (**arc**) tracks, in both genome and transcript coordinates.

Reactivity tracks start from the rf-norm reactivities: `rf-wiggle` converts each sample's XML to WIG, replicates in a `sample_group` are merged (and averaged position-by-position when there is more than one), and the result is written both in transcript coordinates and, using the GTF to remap transcript positions onto the genome, in genomic coordinates, each converted to BigWig. Shannon-entropy tracks follow the same path from the rf-fold entropy profiles, and are produced only when `--rffold_shannon_entropy` is enabled. Base-pairing arcs come from the rf-fold dot-plots: the pairing-probability matrix is converted to base-pair (`.bp`) files that a browser renders as arcs between paired positions.

The genome-coordinate tracks need a GTF for the reference in order to remap transcript positions; if none is available for a reference, its genome track is skipped and only the transcript-coordinate track is written. All of these outputs are skipped when `--stop_after_jackknife` is set.

### rf-jackknife (optional)

`rf-jackknife` calibrates the reactivity-to-folding-constraint conversion by grid-searching slope/intercept values against a known reference structure, choosing the pair that best reproduces it. It runs between rf-norm and rf-fold, and the resulting optimal slope/intercept take priority over both the `--rffold_slope`/`--rffold_intercept` flags and the chemical-based defaults. Enable it by supplying a reference `.db` structure file with `--jackknife_reference`; when omitted, jackknife is skipped and rf-fold uses whatever slope/intercept are otherwise set.

The search space is defined by `--rfjackknife_slope` (default `0,5`) and `--rfjackknife_intercept` (default `-3,0`), stepped by `--rfjackknife_slope_step` and `--rfjackknife_intercept_step` (both default `0.2`). The fit is scored with the Fowlkes-Mallows index; `--rfjackknife_mfmi` (default `true`) uses the modified FMI and `--rfjackknife_relaxed` (default `true`) applies relaxed pairing criteria. `--rfjackknife_keep_lonelypairs` (default `true`) retains 1-bp helices in the reference, and `--rfjackknife_keep_pseudoknots` (default `false`) excludes pseudoknotted base-pairs. Additional rf-fold parameters used inside the search are passed with `--rfjackknife_rf_fold_params` (default `-md 600`).

By default the pipeline pools XMLs from all sample groups into a single jackknife run to derive one consensus slope/intercept; set `--rfjackknife_pool_all false` to run a separate jackknife per group. `--rfjackknife_only_common` restricts the fit to transcripts present across all replicates, and `--rfjackknife_img` writes an R heatmap of the grid-search results. Set `--stop_after_jackknife` to end the run once calibration completes, skipping rf-fold and all downstream outputs.

For the full list of available options see the [rf-jackknife documentation](https://rnaframework-docs.readthedocs.io/en/latest/rf-jackknife/). Any flag not exposed as a pipeline parameter can be passed directly via `ext.args` in a custom config.

### rf-eval (optional)

`rf-eval` assesses how well the reactivity data agrees with a set of known reference structures, giving a quality-control read on the probing signal itself. It compares each sample's rf-norm reactivity profile against the matching reference structure and reports three agreement metrics (each ranging from 0 to 1): the **unpaired coefficient**, the fraction of highly-reactive bases that are unpaired in the reference; **DSCI**, the probability that a randomly chosen unpaired base is more reactive than a randomly chosen paired base; and **AUROC**, how well reactivity discriminates paired from unpaired bases across all thresholds. Enable it by supplying a reference `.db` structure file with `--rfeval_reference`; when omitted, rf-eval is skipped.

Of the three metrics, only AUROC has a fixed chance level (`0.5`). The baselines of DSCI and the unpaired coefficient are structure-specific, so neither can be called good or bad on its own, nor compared between structures. Every run therefore also scores each structure against circular rotations of its own reactivity profile as a comparable negative control and reports the mean and standard deviation of those scores next to each metric.

`--rfeval_reactivity_cutoff` (default `0.7`) sets the reactivity above which a base is classified as highly-reactive when computing the unpaired coefficient. Terminal base-pairs are excluded from the comparison by default (`--rfeval_ignore_terminal true`); alternatively count them as unpaired with `--rfeval_terminal_as_unpaired true`. `--rfeval_keep_pseudoknots` (default `true`) and `--rfeval_keep_lonelypairs` (default `true`) control whether pseudoknotted and isolated 1-bp base-pairs in the reference are retained in the comparison. Set `--rfeval_img true` to additionally write R metric plots.

`--rfeval_windows` takes an optional manifest that slices reactivities to the sub-region each reference structure covers before rf-eval runs, so a short element (e.g. an Rfam family on a whole genome) is not scored against the full transcript. Coordinates are 1-based inclusive in that reference's own space, `strand` is optional (`+` by default), and each `structure_id` must match an entry id in the `--rfeval_reference` `.db` file:

```
# ref_seq_id      start   end     structure_id    [strand]
NC_012532.1        1       161     zika_5UTR       +
NC_012532.1        10383   10454   zika_xrRNA1     +
NC_012532.1        10467   10535   zika_xrRNA2     +
```

For the full list of available options see the [rf-eval documentation](https://rnaframework-docs.readthedocs.io/en/latest/rf-eval/). Any flag not exposed as a pipeline parameter can be passed directly via `ext.args` in a custom config.

### rf-structextract (optional)

Once structures are folded, `rf-structextract` can pull out the high-confidence structural elements: substructures whose bases show consistently low reactivity and low Shannon entropy, the signature of a well-defined, stably folded region, and that meet thermodynamic and geometric criteria. It runs after rf-fold on each `sample_group`, using the fold output (structures + Shannon entropy) together with the group's rf-norm reactivity profiles. Enable it with `--structextract true`.

Two selection tests are toggled individually, and both run by default, so only regions with probing support _and_ high folding confidence are extracted. `--structextract_ignore_react` (default `false`) keeps the reactivity test on, so an unprobed region with no reactivity is not reported; set it `true` for structure-only mode, extracting well-defined motifs regardless of reactivity coverage. `--structextract_ignore_shannon` (default `false`) keeps the Shannon-entropy test on; set it `true` to skip it and relax the folding-confidence requirement.

A 2D diagram is drawn for each extracted motif with ViennaRNA RNAplot, coloured by SHAPE reactivity in the same style as the rf-fold structure plots (the motif's reactivity is sliced from its parent transcript); disable them with `--structextract_plot false`.

For the full list of available options see the [rf-structextract documentation](https://rnaframework-docs.readthedocs.io/en/latest/rf-structextract/).

## Running the pipeline

Typical usage (genome route, auto-downloads reference from Ensembl):

```bash
nextflow run nf-core/rnastructurome \
  --input samplesheet.csv \
  --outdir results \
  -profile docker
```

With local reference files:

```bash
nextflow run nf-core/rnastructurome \
  --input samplesheet.csv \
  --outdir results \
  --fasta genome.fa.gz \
  --gtf annotation.gtf.gz \
  -profile docker
```

Transcriptome route with jackknife calibration:

```bash
nextflow run nf-core/rnastructurome \
  --input samplesheet.csv \
  --outdir results \
  --transcriptome \
  --jackknife_reference ecoli_rrna_calibration/ecoli_k12_rrna_reference_collab.db \
  -profile docker
```

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

For repeated runs with the same parameters, use a params file:

```bash
nextflow run nf-core/rnastructurome -profile docker -params-file params.yaml
```

Example `params.yaml`:

```yaml
input: /path/to/samplesheet.csv
outdir: results
organism: Homo sapiens
method: SHAPE
principle: RT-stop

# trim 2 nt from 5' end to remove sequence-context bias
rfcount_trim_5prime: 2

# use 2-8% normalisation instead of the default box-plot
rfnorm_norm_method: 1
```

## Reproducibility

Always specify the pipeline version when running on your data using `-r` (one hyphen):

```bash
nextflow run nf-core/rnastructurome -r 1.0.0 -profile docker --input samplesheet.csv --outdir results
```

This pins the exact version of the pipeline code and all software containers, so re-running with the same `-r` tag months later will produce identical results. The version is recorded in the MultiQC report and in the RDAT files so you always know what was used.

For complete reproducibility, save your parameters to a `params.yaml` file (see [Running the pipeline](#running-the-pipeline)) and share it alongside your data. When doing so, avoid including cluster-specific paths or institutional profile names; use relative paths or published dataset identifiers instead so that others can reproduce the run in their own environment.

## Core Nextflow arguments

These options are part of Nextflow itself and use a single hyphen. Pipeline parameters use a double hyphen.

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments or analyses.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods; see below. We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility; when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [nf-core/configs](https://github.com/nf-core/configs) at run time, making institutional cluster profiles available automatically. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Multiple profiles can be combined; they are loaded in sequence so later profiles can overwrite earlier ones:

```bash
-profile test,docker
```

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on `PATH`. This is not recommended as it can lead to different results on different machines.

| Profile              | Description                                                                                                                                                                                                                                                 |
| -------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `test`               | Minimal test using the STAR genome-alignment route. Uses human mitochondrial chromosome (MT-RNR1) test data; no other parameters needed.                                                                                                                    |
| `test_transcriptome` | Minimal test using the Bowtie2 transcriptome route. Uses a single-transcript FASTA (ENST00000389680 / MT-RNR1) to exercise the `--transcriptome` path; no other parameters needed.                                                                          |
| `test_prokaryote`    | Prokaryote transcriptome-route test using E. coli 16S rRNA DMS-MaP data and a 16S reference structure for jackknife calibration; no other parameters needed. Significantly faster with `conda` than with container profiles since it includes rf-jackknife. |
| `test_full`          | Full-size test on rice (_Oryza sativa_) in vivo DMS-MaPseq data (GSE197245, one DMS(-) control and three DMS(+) replicates), using the STAR genome route with the reference resolved from Ensembl Plants; no other parameters needed.                       |
| `docker`             | Use [Docker](https://docs.docker.com/engine/installation/) containers.                                                                                                                                                                                      |
| `singularity`        | Use [Singularity](https://www.sylabs.io/guides/3.0/user-guide/) containers.                                                                                                                                                                                 |
| `podman`             | Use [Podman](https://podman.io/) containers.                                                                                                                                                                                                                |
| `shifter`            | Use [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/) containers.                                                                                                                                                                          |
| `charliecloud`       | Use [Charliecloud](https://hpc.github.io/charliecloud/) containers.                                                                                                                                                                                         |
| `apptainer`          | Use [Apptainer](https://apptainer.org/) containers.                                                                                                                                                                                                         |
| `wave`               | Enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ≥ 24.03.0-edge).                                                                                                                                   |
| `conda`              | Use [Conda](https://conda.io/miniconda.html). Please only use Conda as a last resort when containers are not possible. Note that R2DT structure diagrams are not available under conda/mamba; ViennaRNA RNAplot is used as fallback.                        |
| `arm64`              | Applies overrides supplying ARM-compatible containers and Conda environments. See [Running on Linux ARM architectures](#running-on-linux-arm-architectures).                                                                                                |

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file for process resources, executors, or module-specific `ext.args`. See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18), it will automatically be resubmitted with a higher resource request (2× original, then 3× original). If it still fails after the third attempt then the pipeline execution is stopped.

Computationally intensive steps in this pipeline include STAR alignment, rf-count, rf-norm, and rf-fold. These are labelled `process_high` or `process_medium` and will benefit most from tuning.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) sections of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the [nf-core/configs](https://github.com/nf-core/configs) git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the nf-core/configs repository with the addition of your config file, associated documentation file (see examples in [nf-core/configs/docs](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the [main Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

### Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time. Some HPC setups also allow you to run Nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

### Nextflow memory requirements

In some cases, the Nextflow Java virtual machine can start to request a large amount of memory. We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~/.bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Troubleshooting

### General debugging steps

1. Check `.nextflow.log` for the first error message.
2. Check `pipeline_info/execution_trace_<timestamp>.txt` to identify which process failed and its exit status.
3. Use `-resume` after fixing issues; all previously successful tasks will be cached and skipped.
