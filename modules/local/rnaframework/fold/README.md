# RNAFRAMEWORK_RFFOLD

Local module wrapping `rf-fold` from RNAFramework.

Inputs:

- normalized RNAFramework XML reactivity files

Behavior:

- runs RNA secondary-structure prediction from normalized reactivities
- writes the RNAFramework fold output directory
- fails explicitly if RNAFramework reports fold errors or produces no structure files

This module is kept local to the pipeline but is structured to stay close to nf-core module conventions.
