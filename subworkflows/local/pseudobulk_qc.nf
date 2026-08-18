#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {
    PSEUDOBULK_QC_CALC
    PSEUDOBULK_QC_AGGREGATE as PSEUDOBULK_QC_SUMMARY
    PSEUDOBULK_QC_AGGREGATE as PSEUDOBULK_VFAMILY_SUMMARY
} from '../../modules/local/pseudobulk_qc/main'

/*
 * PSEUDOBULK_QC_SW
 *
 * TCRtoolkit QC gate applied to pseudobulk-derived bulk TCR data
 * (SC_TO_BULK in combined Track B, SC_TO_CDR3 in full-SC mode).
 *
 * For every sample it computes n_clones / n_cells and gates on:
 *     n_clones >= params.pseudobulk_qc_min_clones  AND
 *     n_cells  >= params.pseudobulk_qc_min_cells
 *
 * Behaviour on failing samples (params.pseudobulk_qc_mode):
 *     'drop'      (default) — skip failing samples, continue with survivors
 *     'hard_stop'           — abort the run if any sample fails QC
 *
 * Emits:
 *     sample_map  — [meta, file] for QC-passing samples only
 *     concat_cdr3 — concatenated CDR3 file rebuilt from passing samples
 *                   (consumed by ANNOTATE in full-SC mode)
 *     qc_summary  — per-sample QC table (includes dropped samples)
 *     v_family    — per-sample V gene-family usage table
 */
workflow PSEUDOBULK_QC_SW {
    take:
    sample_map   // channel: [meta, file] pseudobulk per-sample TCR tables

    main:
    def min_clones = params.pseudobulk_qc_min_clones
    def min_cells  = params.pseudobulk_qc_min_cells
    def gate_mode  = (params.pseudobulk_qc_mode ?: 'drop').toLowerCase()

    // Wrap the scalar thresholds as explicit value channels. Nextflow 24.10 does not
    // reliably auto-broadcast a bare value passed alongside a queue-channel input.
    PSEUDOBULK_QC_CALC( sample_map, channel.value(min_clones), channel.value(min_cells) )

    branched = PSEUDOBULK_QC_CALC.out.scored.branch { _meta, _file, qc_pass, _nc, _ce ->
        pass: qc_pass == 'PASS'
        fail: true
    }

    // Log every dropped sample
    branched.fail.view { meta, _file, _qc_pass, nc, ce ->
        "[Pseudobulk QC] DROPPED '${meta.sample}': clones=${nc} (min ${min_clones}), cells=${ce} (min ${min_cells})"
    }

    // Optional hard-stop: abort the run if any sample fails QC
    if (gate_mode == 'hard_stop') {
        branched.fail
            .map { meta, _file, _qc_pass, nc, ce -> "${meta.sample} (clones=${nc}, cells=${ce})" }
            .collect()
            .subscribe { failed ->
                if (failed) {
                    throw new RuntimeException(
                        "[Pseudobulk QC] hard-stop: ${failed.size()} sample(s) failed QC -> ${failed.join('; ')}"
                    )
                }
            }
    }

    // Passing samples only -> [meta, file]
    passed_map = branched.pass.map { meta, file, _qc_pass, _nc, _ce -> [meta, file] }

    // Warn (drop mode) if nothing survived the gate
    passed_map.count().subscribe { n ->
        if (n == 0) {
            log.warn "[Pseudobulk QC] No samples passed QC " +
                     "(min_clones=${min_clones}, min_cells=${min_cells}); downstream steps will be skipped."
        }
    }

    // Rebuild concatenated CDR3 from passing samples (for ANNOTATE)
    concat_cdr3 = passed_map
        .map { _meta, file -> file }
        .collectFile(name: 'concat_cdr3.tsv', keepHeader: true, skip: 1)

    // ── Reporting (non-redundant): QC summary + V gene-family usage ──────────
    PSEUDOBULK_QC_SUMMARY( PSEUDOBULK_QC_CALC.out.qc_csv.collect(), 'pseudobulk_qc_summary.csv' )
    PSEUDOBULK_VFAMILY_SUMMARY( PSEUDOBULK_QC_CALC.out.v_family_csv.collect(), 'pseudobulk_v_family.csv' )

    emit:
    sample_map  = passed_map
    concat_cdr3 = concat_cdr3
    qc_summary  = PSEUDOBULK_QC_SUMMARY.out.aggregated_csv
    v_family    = PSEUDOBULK_VFAMILY_SUMMARY.out.aggregated_csv
}
