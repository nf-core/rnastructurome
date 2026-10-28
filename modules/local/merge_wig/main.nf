process MERGE_WIG {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/coreutils:04a9693aaeef79aa'
        : 'community.wave.seqera.io/library/coreutils:9a0aa1447088f229'}"

    input:
    tuple val(meta), path(wig)

    output:
    tuple val(meta), path("*.merged.wig"), emit: merged_wig
    tuple val("${task.process}"), val('cat'), eval("cat --version 2>&1 | head -n 1 | sed 's/^.*coreutils) //; s/ .*\$//'"), topic: versions, emit: versions_cat

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    shopt -s nullglob
    wig_files=(*.wig)
    if [[ \${#wig_files[@]} -eq 0 ]]; then
        echo "No WIG files were provided for merging" >&2
        exit 1
    fi

    {
        echo "track type=wiggle_0"
        for f in "\${wig_files[@]}"; do
            while IFS= read -r line || [[ -n "\$line" ]]; do
                [[ "\$line" == "track "* ]] && continue
                printf '%s\\n' "\$line"
            done < "\$f"
        done
    } > "${prefix}.merged.wig"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged.wig
    """
}
