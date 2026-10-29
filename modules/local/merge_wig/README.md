# MERGE_WIG

Local module that concatenates per-transcript rf-wiggle WIG files into a single track.

Behavior:

- merges per-transcript WIG files into a single `*.merged.wig`
- strips per-file `track` headers and emits one combined `track type=wiggle_0` header
- fails explicitly if no `*.wig` files are staged
