# FASTA_SANITIZE

Local module that normalizes bespoke organism-specific FASTA header id quirks, upstream of `fasta/sort`.

Inputs:

- reference FASTA, plain or gzip-compressed

Behavior:

- normalizes yeast systematic isoform transcript ids from `-A/-B/-C` to `_A/_B/_C`
- replaces characters such as parentheses (e.g. yeast tRNA IDs like `tK(UUU)K`) that hang
  RNAframework's XML parser
- leaves the sequence and any header text after the first token untouched

This is only run for organisms known to need it (see `idSanitizeOrganisms()` in
`subworkflows/local/prepare_references/main.nf`, shared with `gtf/sanitize`) — other
organisms skip the process entirely rather than pay for a no-op sanitize pass.
