process RNAFRAMEWORK_RFSTRUCTEXTRACT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/rnaframework_findutils:a2ed531525463ed6' :
        'community.wave.seqera.io/library/rnaframework_findutils:70fa3c51935922d9' }"

    input:
    tuple val(meta), path(fold_dir), path(xmls, stageAs: 'xml_input/*')

    output:
    tuple val(meta), path("${prefix}_structextract/"), emit: motifs, optional: true
    tuple val(meta), path("${prefix}.rfstructextract.log"), emit: log, optional: true
    tuple val("${task.process}"), val('rnaframework'), eval("rf-structextract -h 2>&1 | grep -oE 'v[0-9]+\\.[0-9]+\\.[0-9]+' | sed 's/v//' | head -1 | grep . || echo unknown"), topic: versions, emit: versions_rnaframework

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def rnaplot = task.ext.rnaplot ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    export TERM="\${TERM:-xterm}"

    # rf-structextract reads the structure files from a "structures/" subfolder of the rf-fold output
    # (dot-bracket .db is accepted; CT is not required). The RFFOLD module renames rf-fold's native
    # structures/ folder to dotbracket/ (and repurposes structures/ for RNAplot SVGs), so reconstruct
    # a minimal rf-fold-style directory here that maps the real structure files back under structures/,
    # alongside the shannon/ folder it needs for Shannon-entropy evaluation.
    rffold_in=rffold_input
    mkdir -p "\${rffold_in}"
    ln -s "\$(cd ${fold_dir}/dotbracket && pwd)" "\${rffold_in}/structures"
    [[ -d ${fold_dir}/shannon ]] && ln -s "\$(cd ${fold_dir}/shannon && pwd)" "\${rffold_in}/shannon"
    [[ -d ${fold_dir}/dotplot ]] && ln -s "\$(cd ${fold_dir}/dotplot && pwd)" "\${rffold_in}/dotplot"

    rf-structextract \\
        -p ${task.cpus} \\
        -ro "\${rffold_in}" \\
        -xf xml_input \\
        -o ${prefix}_structextract \\
        -ow \\
        ${args} 2>&1 | tee "${prefix}.rfstructextract.log"

    # Render an SVG diagram per extracted motif, coloured by SHAPE reactivity to match the rf-fold /
    # VIENNARNA structure plots (same RNAplot + viennarna_colour_svg.py styling). rf-structextract
    # writes one multi-record dot-bracket file per transcript (>id_start-end / seq / structure);
    # RNAplot emits <id>_<start>-<end>_ss.svg per record. Each motif is a sub-range of its transcript,
    # so slice the full-transcript reactivity to [start..end] (re-based to 1) before colouring.
    if [[ -n "${rnaplot}" ]]; then
        for db in ${prefix}_structextract/*.db; do
            [[ -e "\${db}" ]] || continue
            tid=\$(basename "\${db}" .db)
            # Full-transcript per-position reactivity (averaged across any matching replicate XMLs).
            viennarna_extract_xml.py "\${tid}" xml_input/"\${tid}".xml >| "\${tid}.full.shape" 2>/dev/null || true
            ( cd ${prefix}_structextract && "${rnaplot}" --output-format=svg < "\$(basename "\${db}")" ) || true
            for svg in ${prefix}_structextract/"\${tid}"_*_ss.svg; do
                [[ -e "\${svg}" ]] || continue
                region=\$(basename "\${svg}" _ss.svg); region="\${region#\${tid}_}"
                s="\${region%-*}"; e="\${region#*-}"
                if [[ -s "\${tid}.full.shape" ]]; then
                    awk -v s="\${s}" -v e="\${e}" 'BEGIN{FS=OFS="\\t"} \$1>=s && \$1<=e {print \$1-s+1, \$2}' \\
                        "\${tid}.full.shape" >| motif.shape
                    viennarna_colour_svg.py motif.shape "\${svg}" || true
                fi
            done
        done
    fi

    # Organize outputs into the published layout: dot-bracket motif files under dotbracket/, their
    # SVG diagrams under images/. The output directory itself is renamed to extracted_structures/ at
    # publish time (see modules.config).
    mkdir -p ${prefix}_structextract/dotbracket
    find ${prefix}_structextract -maxdepth 1 -name '*.db' -exec mv {} ${prefix}_structextract/dotbracket/ \\;
    if find ${prefix}_structextract -maxdepth 1 -name '*.svg' | grep -q .; then
        mkdir -p ${prefix}_structextract/images
        find ${prefix}_structextract -maxdepth 1 -name '*.svg' -exec mv {} ${prefix}_structextract/images/ \\;
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_structextract/dotbracket ${prefix}_structextract/images
    touch ${prefix}_structextract/dotbracket/${prefix}.db
    touch ${prefix}_structextract/images/${prefix}_1-50_ss.svg
    touch ${prefix}.rfstructextract.log
    """
}
