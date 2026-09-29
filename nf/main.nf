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
        .map { f -> println "wrote ${f}" }
}
