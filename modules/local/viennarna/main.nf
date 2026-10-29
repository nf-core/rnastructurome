process VIENNARNA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/viennarna_findutils_llvm-openmp_python:708e32eb1af5e5a1' :
        'community.wave.seqera.io/library/viennarna_findutils_llvm-openmp_python:d5cfe7db8eb4e492' }"

    input:
    tuple val(meta), path(fold_dir), path(xml, stageAs: "xml_input*/*"), path(drawn_ids)

    output:
    tuple val(meta), path("${prefix}_structures/*.svg"), emit: plots, optional: true
    tuple val("${task.process}"), val('viennarna'), eval("RNAplot --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 | grep . || echo unknown"), topic: versions, emit: versions_viennarna

    when:
    task.ext.when == null || task.ext.when

    script:
    def rnaplot = task.ext.rnaplot ?: 'RNAplot'
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_structures

    sort "${drawn_ids}" > .r2dt_drawn.txt

    _process_db() {
        local _db="\$1"
        local _id=\$(basename "\$_db" .db)
        grep -qxF "\$_id" .r2dt_drawn.txt 2>/dev/null && return
        viennarna_extract_xml.py "\${_id}" xml_input*/"\${_id}".xml > "\${_id}.shape" || true
        "${rnaplot}" --output-format=svg < "\$_db" || true
        if [[ -f "\${_id}_ss.svg" ]]; then
            mv "\${_id}_ss.svg" "${prefix}_structures/\${_id}.svg"
            viennarna_colour_svg.py "\${_id}.shape" "${prefix}_structures/\${_id}.svg" || true
        fi
    }
    export -f _process_db

    # Parallelise per-structure rendering with xargs -P (the same find|xargs pattern R2DT uses),
    # rather than manual `&`/`wait -n` bash job control.
    find ${fold_dir}/dotbracket -maxdepth 1 -name '*.db' -print0 \\
        | xargs -0 -P ${task.cpus} -I{} bash -c '_process_db "\$@"' _ {} \\
        || true
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_structures
    touch ${prefix}_structures/stub_ENST00000000001.svg
    """
}
