process AVERAGE_WIG {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(wigs, stageAs: "inputs/*.wig")

    output:
    tuple val(meta), path("${prefix}.merged.wig"), emit: merged_wig
    tuple val("${task.process}"), val('python'), eval("python3 --version | cut -d' ' -f2"), topic: versions, emit: versions_python

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 - "${prefix}" <<'PY'
import sys
from pathlib import Path


HEADERS = (b"variableStep", b"fixedStep")


def header_fields(line):
    return dict(p.split(b"=", 1) for p in line.split() if b"=" in p)


# Byte offsets of each chromosome's block headers, so one chromosome can be read at a time.
def index_blocks(filepath):
    blocks = {}
    offset = 0
    with filepath.open("rb") as fh:
        for line in fh:
            if line.startswith(HEADERS):
                blocks.setdefault(header_fields(line)[b"chrom"].decode(), []).append(offset)
            offset += len(line)
    return blocks


def read_chrom(filepath, offsets):
    data = {}
    with filepath.open("rb") as fh:
        for offset in offsets:
            fh.seek(offset)
            header = fh.readline()
            parts = header_fields(header)
            is_fixed = header.startswith(b"fixedStep")
            pos = int(parts.get(b"start", 1))
            step = int(parts.get(b"step", 1))
            for line in fh:
                line = line.strip()
                if not line or line.startswith(b"track"):
                    continue
                if line.startswith(HEADERS):
                    break
                if is_fixed:
                    data[pos] = float(line)
                    pos += step
                else:
                    p, v = line.split()
                    data[int(p)] = float(v)
    return data


prefix = sys.argv[1]
wig_files = sorted(Path("inputs").glob("*.wig"))
output_path = Path(f"{prefix}.merged.wig")

if len(wig_files) == 1:
    output_path.write_bytes(wig_files[0].read_bytes())
else:
    # Holding every replicate genome-wide in memory was OOM-killed at 96 GB on 6 human replicates.
    indexes = [(f, index_blocks(f)) for f in wig_files]
    all_chroms = set().union(*(blocks.keys() for _, blocks in indexes))
    with output_path.open("w") as fh:
        fh.write("track type=wiggle_0\\n")
        for chrom in sorted(all_chroms):
            sums = {}
            counts = {}
            for f, blocks in indexes:
                if chrom not in blocks:
                    continue
                for pos, val in read_chrom(f, blocks[chrom]).items():
                    sums[pos] = sums.get(pos, 0.0) + val
                    counts[pos] = counts.get(pos, 0) + 1
            fh.write(f"variableStep chrom={chrom}\\n")
            for pos in sorted(sums):
                fh.write(f"{pos} {sums[pos]/counts[pos]:.6g}\\n")
PY
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged.wig
    """
}
