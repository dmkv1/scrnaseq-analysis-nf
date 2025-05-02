#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Droplet processing
include { DROPLETS_TO_CELLS } from './modules/qc'
include { DOUBLET_DETECTION } from './modules/qc'
include { CELL_QC } from './modules/qc'
// Cell annotations
include { SEURAT_CLUSTERING } from './modules/analysis'

workflow {
    seed = Channel.value(params.seed)

    // Validate input parameters
    if (params.input == null) {
        exit(1, "Input samplesheet not specified!")
    }

    // Create channel from input samplesheet
    ch_input = Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def sample_id = row.sample
            def expected_cells = row.expected_cells.toInteger()
            def patient_id = row.patient
            def timepoint = row.timepoint
            def compartment = row.compartment
            def counts = file(row.counts)

            return [sample_id, expected_cells, patient_id, timepoint, compartment, counts]
        }

    // Execute QC modules
    DROPLETS_TO_CELLS(
        file("templates/1_droplets_to_cells.Rmd"),
        ch_input,
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

    // Cluster cells
    SEURAT_CLUSTERING(
        file("templates/4_seurat_clustering.Rmd"),
        CELL_QC.out.sce,
        seed,
    )
}
