# ENSEMBL_TRANSCRIPTOME

Local module to resolve Ensembl transcript FASTA files for a species and merge:

- `cdna.all.fa.gz`
- optional `ncrna.fa.gz`

Behavior:

- accepts Ensembl species format such as `homo_sapiens`
- also accepts Latin binomials such as `Homo sapiens`
- normalizes whitespace to underscores before lookup
- writes `ensembl_source_url.txt` with the URLs used
- writes `ensembl_warnings.log` when `ncrna` is unavailable

This module is kept local to the pipeline but is structured to stay close to nf-core module conventions.
