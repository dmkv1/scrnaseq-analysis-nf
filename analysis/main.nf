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

    ch_patients= Channel.fromPath(params.patients, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def patient_id = row.patient
            def k_obs_groups = row.k_obs_groups.toInteger()

            return [patient_id, k_obs_groups]
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
        DOUBLET_DETECTION.out.sce,
        seed,
    )

    MERGE(
        CELL_QC.out.sce.collect(),
        seed
    )

    ANNOTATE(
        file("templates/5_annotation.Rmd"),
        MERGE.out.merged_sce,
        seed
    )

    INFERCNV(
        ch_patients,
        ANNOTATE.out.annotated_sce,
        seed
    )
}
