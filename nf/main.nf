#!/usr/bin/env nextflow
/* FCS QC demo pipeline.

   FCSIO_DEMO -> FCSIO_PARSE -> FCS_STATS

   With --input unset, FCSIO_DEMO writes deterministic gated-population
   FCS fixtures (fcsio writer + cytof_qc demo generator), then every
   file is parsed and summarized per-channel. With --input pointing at
   a directory of *.fcs the demo step is skipped and real files flow
   through the same parse/stats DAG.
*/

include { FCSIO_DEMO  } from './modules/local/fcsio_demo/main.nf'
include { FCSIO_PARSE } from './modules/local/fcsio_parse/main.nf'
include { FCS_STATS   } from './modules/local/fcs_stats/main.nf'

workflow FCS_QC {
    main:
    if (params.input) {
        ch_fcs = channel.fromPath("${params.input}/*.fcs")
            .map { f -> [[id: f.baseName.replaceAll('\\s+', '_')], f] }
    } else {
        FCSIO_DEMO(params.demo_events)
        ch_fcs = FCSIO_DEMO.out.fcs
            .flatten()
            .map { f -> [[id: f.baseName.replaceAll('\\s+', '_')], f] }
    }

    FCSIO_PARSE(ch_fcs)
    FCS_STATS(FCSIO_PARSE.out.events)

    emit:
    stats   = FCS_STATS.out.stats
    meta    = FCSIO_PARSE.out.meta
}

workflow {
    FCS_QC()

    FCS_QC.out.stats
        .collectFile(
            name: "qc_summary.jsonl",
            storeDir: params.outdir,
            newLine: true) { meta, f ->
                new groovy.json.JsonBuilder(
                    [sample: meta.id,
                     stats: new groovy.json.JsonSlurper().parseText(f.text)]
                ).toString() }
        .map { f ->
            // f is the written summary; emit the provenance manifest
            def outdir = f.getParent()
            def cmd = ["python3",
                       "${projectDir}/../scripts/run_manifest.py",
                       "--pipeline", "fcs_qc",
                       "--engine", "nextflow@${nextflow.version}",
                       "--out", "${outdir}/run_manifest.json",
                       "--artifacts", outdir.toString(),
                       "--param", "input=${params.input ?: 'demo'}",
                       "--param", "demo_events=${params.demo_events}",
                       "--param", "run_name=${workflow.runName}",
                       "--param", "session=${workflow.sessionId}"]
            if (params.input) {
                cmd.addAll(["--inputs", params.input as String])
            }
            def p = cmd.execute()
            p.waitFor()
            println "wrote ${f}" + (p.exitValue() == 0
                ? " + run_manifest.json"
                : " (manifest failed: ${p.err.text})")
            f }

    // Hash-chained execution receipts (scripts/run_receipts.py): one
    // JSONL per run, one record per task binding staged inputs, written
    // outputs, env fingerprint, and timestamps. Best-effort — the trace
    // writer may still be flushing; the deterministic path is the
    // documented post-run command, which is what CI uses.
    // params/workflow/nextflow are not bound inside onComplete's
    // delegate — capture them in locals now.
    def outdirForReceipts = params.outdir
    def runName = workflow.runName
    def sessionId = "${workflow.sessionId}"
    def nfVersion = "${nextflow.version}"
    def receiptsScript = "${projectDir}/../scripts/run_receipts.py"
    workflow.onComplete {
        def cmd = ["python3", receiptsScript,
                   "--trace", "${outdirForReceipts}/trace.txt",
                   "--out", "${outdirForReceipts}/run_receipts.jsonl",
                   "--run-name", runName,
                   "--session-id", sessionId,
                   "--engine", "nextflow@${nfVersion}"]
        def p = cmd.execute()
        p.waitFor()
        println (p.exitValue() == 0 ?
            "wrote run_receipts.jsonl" :
            "receipts deferred (run scripts/run_receipts.py): ${p.err.text}")
    }
}
