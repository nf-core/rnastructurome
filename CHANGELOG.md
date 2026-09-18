# nf-core/rnastructurome: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.1dev - [unreleased]

## v1.0.0 - [8 September 2026]

First release of nf-core/rnastructurome, which analyses chemical high-throughput RNA structure-probing data and predicts RNA secondary structures from it. The pipeline covers **SHAPE** and **DMS** chemistries read out by either the **RT-stop** or **mutational profiling (MaP)** principle, and runs in two modes depending on the reference: a **genome** route aligned with STAR and a **transcriptome** route aligned with Bowtie/Bowtie2.

### `Added`

- Single- and paired-end FASTQ input, with re-sequenced samples concatenated automatically.
- Read quality control with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) before and after adapter and quality trimming with [Cutadapt](https://cutadapt.readthedocs.io/).
- Optional UMI extraction and UMI-aware deduplication with [UMI-tools](https://umi-tools.readthedocs.io/); position-based duplicate removal with [SAMtools markdup](https://www.htslib.org/) is available but off by default.
- Reference genomes and annotations downloaded automatically from [Ensembl](https://www.ensembl.org/), falling back to [NCBI](https://www.ncbi.nlm.nih.gov/datasets/) for bacteria, viruses and other organisms Ensembl does not cover, or supplied directly as local files.
- Genome-route alignment with [STAR](https://github.com/alexdobin/STAR), and transcriptome-route alignment with [Bowtie](https://bowtie-bio.sourceforge.net/) for RT-stop and [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/) for MaP.
- BAM sorting, indexing and alignment statistics with [SAMtools](https://www.htslib.org/).
- Per-base reactivity counting with [rf-count](https://rnaframework-docs.readthedocs.io/en/latest/rf-count/), tallying mutations for MaP and RT-stops for RT-stop on transcript coordinates.
- An alternative genome-coordinate counting path with [rf-count-genome](https://rnaframework-docs.readthedocs.io/en/latest/rf-count-genome/), resolving library strandedness with [BEDOPS](https://bedops.readthedocs.io/) and [RSeQC](https://rseqc.sourceforge.net/) and extracting per-transcript reactivity with [rf-rctools](https://rnaframework-docs.readthedocs.io/en/latest/rf-rctools/).
- Reactivity normalisation with [rf-norm](https://rnaframework-docs.readthedocs.io/en/latest/rf-norm/), pairing treated, untreated and denatured samples automatically and selecting scoring and normalisation methods from the controls present.
- Replicate reproducibility QC with [rf-correlate](https://rnaframework-docs.readthedocs.io/en/latest/rf-correlate/), reporting pairwise Pearson and Spearman correlation of reactivity profiles.
- RNA secondary structure prediction across grouped replicates with [rf-fold](https://rnaframework-docs.readthedocs.io/en/latest/rf-fold/), using chemistry- and reagent-aware folding defaults.
- Reactivity, Shannon entropy and base-pair arc tracks in transcript and genome coordinates, generated with [rf-wiggle](https://rnaframework-docs.readthedocs.io/en/latest/rf-wiggle/) and [UCSC wigToBigWig](https://genome.ucsc.edu/).
- 2D structure diagrams coloured by reactivity, drawn with [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) RNAplot and, where a template model exists for the RNA type, with [R2DT](https://github.com/RNAcentral/R2DT) alongside them. R2DT is container-only, so `-profile conda` draws with ViennaRNA alone.
- RMDB-compatible RDAT export combining per-transcript reactivity and structure.
- Optional folding calibration against known structures with [rf-jackknife](https://rnaframework-docs.readthedocs.io/en/latest/rf-jackknife/), which tunes the slope and intercept passed to rf-fold.
- Optional structural-element extraction of high-confidence, low-reactivity and low-Shannon motifs with [rf-structextract](https://rnaframework-docs.readthedocs.io/en/latest/rf-structextract/).
- Optional structure-accuracy evaluation with [rf-eval](https://rnaframework-docs.readthedocs.io/en/latest/rf-eval/), reporting AUROC, DSCI and the unpaired coefficient against an automatic rotation baseline.
- Aggregated quality-control report with [MultiQC](http://multiqc.info/), including alignment, reactivity and replicate-correlation summaries.
- Test profiles for the genome route, the transcriptome route, a prokaryote run with jackknife calibration, and a full-size rice DMS-MaPseq dataset.

### `Fixed`

### `Dependencies`

### `Deprecated`
