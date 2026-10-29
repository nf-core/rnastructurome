# GTF_SANITIZE

Local module that strips regex/shell-unsafe characters from GTF `transcript_id`/`gene_id` attributes.

Inputs:

- GTF annotation

Behavior:

- rewrites `transcript_id`/`gene_id` attribute values, replacing characters such as
  parentheses (e.g. yeast tRNA IDs like `tK(UUU)K`) that hang RNAframework's XML parser
- leaves all other attributes and comment lines untouched
- writes a `.sanitized.gtf`

This is only run for organisms known to need it (see `gtfSanitizeOrganisms()` in
`subworkflows/local/prepare_references/main.nf`) — other organisms skip the process
entirely rather than pay for a no-op sanitize pass.
