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


def parse_wig(filepath):
    data = {}
    chrom = None
    is_fixed = False
    fixed_start = 1
    fixed_step = 1
    pos_counter = 0
    for line in filepath.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("track"):
            continue
        if line.startswith("variableStep"):
            parts = dict(p.split("=") for p in line.split() if "=" in p)
            chrom = parts.get("chrom")
            is_fixed = False
            data.setdefault(chrom, {})
        elif line.startswith("fixedStep"):
            parts = dict(p.split("=") for p in line.split() if "=" in p)
            chrom = parts.get("chrom")
            fixed_start = int(parts.get("start", 1))
            fixed_step = int(parts.get("step", 1))
            is_fixed = True
            pos_counter = 0
            data.setdefault(chrom, {})
        else:
            if is_fixed:
                pos = fixed_start + pos_counter * fixed_step
                pos_counter += 1
                data[chrom][pos] = float(line)
            else:
                p, v = line.split()
                data[chrom][int(p)] = float(v)
    return data


prefix = sys.argv[1]
wig_files = sorted(Path("inputs").glob("*.wig"))
output_path = Path(f"{prefix}.merged.wig")

if len(wig_files) == 1:
    output_path.write_bytes(wig_files[0].read_bytes())
else:
    all_data = [parse_wig(f) for f in wig_files]
    all_chroms = set().union(*(d.keys() for d in all_data))
    with output_path.open("w") as fh:
        fh.write("track type=wiggle_0\\n")
        for chrom in sorted(all_chroms):
            all_pos = set().union(*(d[chrom].keys() for d in all_data if chrom in d))
            fh.write(f"variableStep chrom={chrom}\\n")
            for pos in sorted(all_pos):
                vals = [d[chrom][pos] for d in all_data if chrom in d and pos in d[chrom]]
                fh.write(f"{pos} {sum(vals)/len(vals):.6g}\\n")
PY
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged.wig
    """
}
