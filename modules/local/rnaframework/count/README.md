# RNAFRAMEWORK_RFCOUNT

Local module wrapping `rf-count` from RNAFramework.

Inputs:

- mapped BAM
- BAM index
- transcript FASTA reference

Behavior:

- counts RT-stop events or MaP mutations depending on upstream metadata
- writes raw count files and optional index sidecars
- emits optional diagnostic files from RNAFramework, including `error.out`, `samtools.log`, and `index.rci`
- emits plot PDFs when RNAFramework produces them
- fails explicitly when RNAFramework reports zero covered transcripts

This module is kept local to the pipeline but is structured to stay close to nf-core module conventions.
