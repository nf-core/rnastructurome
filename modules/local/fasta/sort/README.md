# FASTA_SORT

Local module that sorts FASTA records lexicographically by sequence identifier, via `seqkit sort`.

Inputs:

- reference FASTA, plain or gzip-compressed

Behavior:

- sorts records by id (`seqkit sort` default, not the full header)
- always writes a plain-text `.sorted.fa`, wrapped at 80 characters, regardless of whether the
  input was gzip-compressed — STAR/samtools faidx reject gzip FASTA, so this stays decompressed
  rather than mirroring `seqkit sort`'s own default of preserving the input's compression state

This is used to keep transcript references deterministic before index building and downstream counting.
Organism-specific header quirks (e.g. yeast isoform ids, unsafe characters) are handled upstream by
`fasta/sanitize` — see its README.
