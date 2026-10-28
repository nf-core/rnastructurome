# RNAFRAMEWORK_RFNORM

Local module wrapping `rf-norm` from RNAFramework.

Inputs:

- treated raw count files
- optional untreated raw count files
- optional denatured raw count files
- optional `rci` sidecar files

Behavior:

- normalizes RNAFramework raw count files to XML reactivities
- supports treated-only, treated/untreated, and denatured control workflows
- emits plot PDFs when RNAFramework produces them

This module is kept local to the pipeline but is structured to stay close to nf-core module conventions.
