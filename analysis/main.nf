#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Processing by sample
include { DROPLETS_TO_CELLS } from './modules/qc'
include { DOUBLET_DETECTION } from './modules/qc'
include { CELL_QC } from './modules/qc'
// Merged analysis
include { MERGE } from './modules/analysis'
include { ANNOTATE } from './modules/analysis'
include { INFERCNV } from './modules/analysis'

workflow {
    seed = Channel.value(params.seed)

    // Create input channels from metadata sheets
    ch_samples = Channel.fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sample_id = row.sample
            def expected_cells = row.expected_cells.toInteger()
            def patient_id = row.patient
            def timepoint = row.timepoint
            def compartment = row.compartment
            def replicate = row.replicate
            def counts = file(row.counts)

            return [sample_id, expected_cells, patient_id, timepoint, compartment, replicate, counts]
        }

    ch_patients = Channel.fromPath(params.patients, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def patient_id = row.patient
            def k_obs_groups = row.k_obs_groups.toInteger()

            return [patient_id, k_obs_groups]
        }

    // Create QC-specific channel from the same CSV
    ch_qc_params = Channel.fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sample_id = row.sample
            // Get individual QC parameters with defaults if not specified
            def nUMI_thresh = row.containsKey('nUMI_thresh') && row.nUMI_thresh?.trim()
                ? row.nUMI_thresh.toInteger()
                : params.qc.nUMI_thresh

            def nGenes_thresh = row.containsKey('nGenes_thresh') && row.nGenes_thresh?.trim()
                ? row.nGenes_thresh.toInteger()
                : params.qc.nGenes_thresh

            def mito_thresh = row.containsKey('mito_thresh') && row.mito_thresh?.trim()
                ? row.mito_thresh.toFloat()
                : params.qc.mito_thresh

            return [sample_id, nUMI_thresh, nGenes_thresh, mito_thresh]
        }

    // Execute QC modules
    DROPLETS_TO_CELLS(
        file("templates/1_droplets_to_cells.Rmd"),
        ch_samples,
        seed,
    )

    DOUBLET_DETECTION(
        file("templates/2_doublets.Rmd"),
        DROPLETS_TO_CELLS.out.sce,
        seed,
    )

    CELL_QC(
        file("templates/3_cell_qc.Rmd"),
        DOUBLET_DETECTION.out.sce.join(ch_qc_params),
        seed,
    )

    MERGE(
        CELL_QC.out.sce.collect()
    )

//    ANNOTATE(
//        file("templates/5_annotation.Rmd"),
//        MERGE.out.merged_sce,
//        seed
//    )
//
//    INFERCNV(
//        ch_patients,
//        ANNOTATE.out.annotated_sce,
//        seed
//    )
//
}
