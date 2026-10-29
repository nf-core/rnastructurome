process WIG_CHROM_SIZES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(wig)

    output:
    tuple val(meta), path("${prefix}.chrom.sizes"), emit: sizes
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/^Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 - "${wig}" "${prefix}.chrom.sizes" <<'PY'
import re
import sys


def _int_header(line, key, default):
    m = re.search(rf"{key}=(\\d+)", line)
    return int(m.group(1)) if m else default


def main():
    wig_path, output_path = sys.argv[1], sys.argv[2]

    sizes = {}

    # State for the current block
    chrom = None
    mode = None   # "fixed" or "variable"
    start = 1
    step = 1
    span = 1
    count = 0
    max_pos = 0

    def _flush():
        if chrom is None:
            return
        if mode == "fixed":
            if count == 0:
                return
            end = start + (count - 1) * step + span - 1
        else:
            if max_pos == 0:
                return
            end = max_pos + span - 1
        sizes[chrom] = max(sizes.get(chrom, 0), end)

    with open(wig_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("track") or line.startswith("#"):
                continue

            if line.startswith("fixedStep"):
                _flush()
                chrom = re.search(r"chrom=(\\S+)", line).group(1)
                mode = "fixed"
                start = _int_header(line, "start", 1)
                step = _int_header(line, "step", 1)
                span = _int_header(line, "span", 1)
                count = 0
                max_pos = 0

            elif line.startswith("variableStep"):
                _flush()
                chrom = re.search(r"chrom=(\\S+)", line).group(1)
                mode = "variable"
                span = _int_header(line, "span", 1)
                count = 0
                max_pos = 0

            else:
                if mode == "fixed":
                    count += 1
                elif mode == "variable":
                    pos = int(line.split()[0])
                    if pos > max_pos:
                        max_pos = pos

    _flush()

    with open(output_path, "w") as out:
        for name in sorted(sizes):
            out.write(f"{name}\\t{sizes[name]}\\n")

    return 0


raise SystemExit(main())
PY
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.chrom.sizes"
    """
}
