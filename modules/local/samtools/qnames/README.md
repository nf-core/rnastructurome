# SAMTOOLS_QNAMES

Local module that extracts the sorted, unique set of read names (QNAMEs) present in a BAM file.

Inputs:

- BAM file

Behavior:

- runs `samtools view | cut -f1 | sort -u`
- writes the resulting name list to a plain-text `.qnames.txt`

This is used in the genome route (STAR) to reconcile STAR's `--quantMode TranscriptomeSAM` output
(produced at alignment time, before dedup) against the deduplicated genome-coordinate BAM: the
qnames file becomes the `--qname-file` input to `SAMTOOLS_VIEW`, filtering the transcript BAM down
to only the reads that survived dedup.
