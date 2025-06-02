#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MATRIX_TO_SCE } from './modules/utils'
include { MATRIX_TO_SCE as CLEAN_MATRIX_TO_SCE } from './modules/utils'
// Processing by sample
include { DROPLETS_TO_CELLS } from './modules/qc'
include { AMBIENT_RNA } from './modules/qc'
include { DOUBLET_DETECTION } from './modules/qc'
include { CELL_QC } from './modules/qc'
// Merged analysis
include { MERGE } from './modules/analysis'
include { LABEL_TUMOR_CELLS } from './modules/analysis'
include { INFERCNV } from './modules/analysis'

workflow {
    seed = Channel.value(params.seed)

    // Create input channels from metadata sheets
    ch_input_samples = Channel.fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sample_id = row.sample
            def patient_id = row.patient
            def timepoint = row.timepoint
            def compartment = row.compartment
            def replicate = row.replicate
            def raw_fbmtx = file(row.raw_feature_bc_matrix)

            return [sample_id, patient_id, timepoint, compartment, replicate, raw_fbmtx]
        }

    ch_metadata = ch_input_samples.map { sample_id, patient_id, timepoint, compartment, replicate, _raw_fbmtx ->
        return [sample_id, patient_id, timepoint, compartment, replicate]
    }

    MATRIX_TO_SCE(ch_input_samples)

    ch_droplets = Channel.fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sample_id = row.sample
            def expected_cells = row.expected_cells.toInteger()

            return [sample_id, expected_cells]
        }

    DROPLETS_TO_CELLS(
        file("templates/1_droplets_to_cells.Rmd"),
        MATRIX_TO_SCE.out.sce.join(ch_droplets),
        seed,
    )

    ch_raw_fbmtx = ch_input_samples.map { sample_id, _patient_id, _timepoint, _compartment, _replicate, raw_fbmtx ->
        return [sample_id, raw_fbmtx]
    }

    AMBIENT_RNA(
        DROPLETS_TO_CELLS.out.filtered_fbmtx.join(ch_raw_fbmtx),
        seed,
    )

    ch_clean_fbmtx = AMBIENT_RNA.out.clean_fbmtx
        .join(ch_metadata)
        .map { sample_id, clean_fbmtx, patient_id, timepoint, compartment, replicate ->
            return [sample_id, patient_id, timepoint, compartment, replicate, clean_fbmtx]
        }

    CLEAN_MATRIX_TO_SCE(ch_clean_fbmtx)

    DOUBLET_DETECTION(
        file("templates/2_doublets.Rmd"),
        CLEAN_MATRIX_TO_SCE.out.sce,
        seed,
    )

    ch_qc_params = Channel.fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sample_id = row.sample

            def nUMI_thresh = row.containsKey('nUMI_thresh') && row.nUMI_thresh?.trim()
                ? row.nUMI_thresh.toInteger()
                : params.qc.nUMI_thresh

            def nGenes_thresh = row.containsKey('nGenes_thresh') && row.nGenes_thresh?.trim()
                ? row.nGenes_thresh.toInteger()
                : params.qc.nGenes_thresh

            def mito_thresh = row.containsKey('mito_thresh') && row.mito_thresh?.trim()
                ? row.mito_thresh.toFloat()
                : params.qc.mito_thresh

            def contam_thresh = row.containsKey('contam_thresh') && row.contam_thresh?.trim()
                ? row.contam_thresh.toFloat()
                : params.qc.contam_thresh

            return [sample_id, nUMI_thresh, nGenes_thresh, mito_thresh, contam_thresh]
        }

    CELL_QC(
        file("templates/3_cell_qc.Rmd"),
        file("bin/process_sce.R"),
        DOUBLET_DETECTION.out.sce.join(ch_qc_params).join(AMBIENT_RNA.out.perCellCont),
        seed,
    )

    MERGE(
        CELL_QC.out.sce.collect()
    )

    if (params.tumor_labelled) {
        LABEL_TUMOR_CELLS(
            file("templates/4_label_tumor_cells.Rmd"),
            file("bin/process_sce.R"),
            MERGE.out.merged_sce,
            seed,
        )

        ch_patients = Channel.fromPath(params.patients, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def patient_id = row.patient_id
                def k_obs_groups = row.k_obs_groups.toInteger()

                return [patient_id, k_obs_groups]
            }
    }
}
