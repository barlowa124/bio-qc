process FCSIO_DEMO {
    tag "demo"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/python:3.12-slim'

    input:
    val n_events

    output:
    path "*.fcs", emit: fcs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // generate deterministic gated-population FCS fixtures via the
    // sibling packages (fcsio writer + cytof_qc demo module)
    """
    export PYTHONPATH="${projectDir}/../fcs_io/src:\${PYTHONPATH:-}"
    mkdir -p demo_out
    make_demo_fcs.py demo_out ${n_events}
    mv demo_out/*.fcs .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch demo.fcs
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
