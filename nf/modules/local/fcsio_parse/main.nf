process FCSIO_PARSE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/python:3.12-slim' :
        'docker.io/python:3.12-slim' }"

    input:
    // stageAs sidesteps Nextflow's backslash-escaping of spaces in
    // vendor filenames (e.g. "Demo_CD8 T_cells.fcs")
    tuple val(meta), path(fcs, stageAs: 'input.fcs')

    output:
    tuple val(meta), path("*.events.tsv"), emit: events
    tuple val(meta), path("*.meta.json"),  emit: meta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // fcsio is a sibling package in this monorepo, not on PyPI — expose
    // its src/ on PYTHONPATH inside the task dir
    """
    export PYTHONPATH="${projectDir}/../fcs_io/src:\${PYTHONPATH:-}"
    fcs_to_tsv.py input.fcs -o "${prefix}.events.tsv" --meta "${prefix}.meta.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        fcsio: "monorepo path import"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.events.tsv" "${prefix}.meta.json"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
        fcsio: "stub"
    END_VERSIONS
    """
}
