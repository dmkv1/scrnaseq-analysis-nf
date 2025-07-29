#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MATRIX_TO_SCE } from './modules/utils'
include { MATRIX_TO_SCE as CLEAN_MATRIX_TO_SCE } from './modules/utils'

workflow {
    seed = Channel.value(params.seed)

    // Create input channels from metadata sheets
    ch_input_samples = Channel.fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sample_id = row.sample
            def patient_id = row.patient
            def label = row.label
            def replicate = row.replicate
            def raw_fbmtx = file(row.raw_feature_bc_matrix)

            return [sample_id, patient_id, label, replicate, raw_fbmtx]
        }

    ch_metadata = ch_input_samples.map { sample_id, patient_id, label, replicate, _raw_fbmtx ->
        return [sample_id, patient_id, label, replicate]
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

    ch_raw_fbmtx = ch_input_samples.map { sample_id, _patient_id, _label, _replicate, raw_fbmtx ->
        return [sample_id, raw_fbmtx]
    }

    AMBIENT_RNA(
        DROPLETS_TO_CELLS.out.filtered_fbmtx.join(ch_raw_fbmtx),
        seed,
    )

    ch_clean_fbmtx = AMBIENT_RNA.out.clean_fbmtx
        .join(ch_metadata)
        .map { sample_id, clean_fbmtx, patient_id, label, replicate ->
            return [sample_id, patient_id, label, replicate, clean_fbmtx]
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
        file("bin/sc_functions.R"),
        DOUBLET_DETECTION.out.sce.join(ch_qc_params).join(AMBIENT_RNA.out.perCellCont),
        seed,
    )

    MERGE(
        CELL_QC.out.sce.collect()
    )

    JOINT_CLUSTERING(
        file("templates/4_clustering.Rmd"),
        file("bin/sc_functions.R"),
        file(params.samples),
        file("marker_genes/immunological_markers.csv"),
        MERGE.out.merged_sce,
        seed,
    )
}

process DROPLETS_TO_CELLS {
    tag "${sample_id}"
    publishDir "${params.outdir}/QC/${sample_id}/1_droplets", mode: 'copy'

    input:
    path 'droplets_to_cells.Rmd'
    tuple val(sample_id), path(sce), val(expected_cells)
    val seed

    output:
    tuple val(sample_id), path("${sample_id}_droplets_to_cells.nb.html"), emit: report
    tuple val(sample_id), path("SCE_${sample_id}_cells.rds"), emit: sce
    tuple val(sample_id), path('filtered_feature_bc_matrix/'), emit: filtered_fbmtx
    tuple val(sample_id), path("${sample_id}_droplet_metrics.json"), emit: metrics

    script:
    """
    Rscript -e "rmarkdown::render('droplets_to_cells.Rmd',
                output_file = '${sample_id}_droplets_to_cells.nb.html',
                output_format = 'html_notebook',
                output_options = list(
                    self_contained = TRUE,
                    df_print = 'paged',
                    code_folding = 'hide',
                    toc = TRUE,
                    toc_float = TRUE
                    ),
                params = list(
                               sample_id = '${sample_id}',
                               path_sce_input = '${sce}',
                               expected_cells = '${expected_cells}',
                               FDR_thresh = '${params.droplets.FDR_thresh}',
                               path_sce_output = 'SCE_${sample_id}_cells.rds',
                               seed = ${seed}
                             ))"
    """
}

process AMBIENT_RNA {
    tag "${sample_id}"
    publishDir "${params.outdir}/QC/${sample_id}/2_ambientRNA", mode: 'copy'

    input:
    tuple val(sample_id), path(filtered_mtx_path), path(raw_mtx_path)
    val seed

    output:
    tuple val(sample_id), path("decontX_feature_bc_matrix/"), emit: clean_fbmtx
    tuple val(sample_id), path("${sample_id}_perCellCont.rds"), emit: perCellCont

    script:
    """
    run_decontX.R \
        --cells_fbmtx "${filtered_mtx_path}" \
        --droplets_fbmtx "${raw_mtx_path}" \
        --use_empty TRUE \
        --decont_fbmtx "decontX_feature_bc_matrix" \
        --perCell "${sample_id}_perCellCont.rds" \
        --seed ${seed}
    """
}

process DOUBLET_DETECTION {
    tag "${sample_id}"
    publishDir "${params.outdir}/QC/${sample_id}/3_doublets", mode: 'copy'

    input:
    path 'doublets.Rmd'
    tuple val(sample_id), path(sce)
    val seed

    output:
    tuple val(sample_id), path("${sample_id}_doublets.nb.html"), emit: report
    tuple val(sample_id), path("SCE_${sample_id}_singlets.rds"), emit: sce
    tuple val(sample_id), path("${sample_id}_doublet_metrics.json"), emit: metrics

    script:
    """
    Rscript -e "rmarkdown::render('doublets.Rmd',
                output_file = '${sample_id}_doublets.nb.html',
                output_format = 'html_notebook',
                output_options = list(
                    self_contained = TRUE,
                    df_print = 'paged',
                    code_folding = 'hide',
                    toc = TRUE,
                    toc_float = TRUE
                ),
                params = list(
                               sample_id = '${sample_id}',
                               path_sce_input = '${sce}',
                               path_sce_output = 'SCE_${sample_id}_singlets.rds',
                               seed = ${seed}
                             ))"
    """
}

process CELL_QC {
    tag "${sample_id}"
    publishDir "${params.outdir}/QC/${sample_id}/4_cellQC", mode: 'copy'

    input:
    path 'cell_qc.Rmd'
    path 'sc_functions.R'
    tuple val(sample_id), path(sce), val(nUMI_thresh), val(nGenes_thresh), val(mito_thresh), val(contam_thresh), path(perCellCont)
    val seed

    output:
    tuple val(sample_id), path("${sample_id}_cell_QC.nb.html"), emit: report
    path ("SCE_${sample_id}.rds"), emit: sce
    path ("SCE_${sample_id}_cleaned.rds"), emit: sce_cleaned
    tuple val(sample_id), path("${sample_id}_cell_QC_metrics.json"), emit: metrics

    script:
    """
    Rscript -e "rmarkdown::render('cell_qc.Rmd',
                output_file = '${sample_id}_cell_QC.nb.html',
                output_format = 'html_notebook',
                output_options = list(
                    self_contained = TRUE,
                    df_print = 'paged',
                    code_folding = 'hide',
                    toc = TRUE,
                    toc_float = TRUE
                ),
                params = list(
                    sample_id = '${sample_id}',
                    path_sce_input = '${sce}',
                    path_sce_output = 'SCE_${sample_id}.rds',
                    path_sce_output_cleaned = 'SCE_${sample_id}_cleaned.rds',
                    path_perCellCont = '${sample_id}_perCellCont.rds',
                    nUMI_thresh = '${nUMI_thresh}',
                    nGenes_thresh = '${nGenes_thresh}',
                    mito_thresh = '${mito_thresh}',
                    contam_thresh = '${contam_thresh}',
                    cluster_discard_thresh = '${params.qc.cluster_discard_thresh}',
                    process_fun = 'sc_functions.R',
                    seed = ${seed}
                    ))"
    """
}

process MERGE {
    publishDir "${params.outdir}/SCE", mode: 'copy'

    input:
    path qc_passed_sces

    output:
    path "SCE_merged.rds", emit: merged_sce

    script:
    """
    merge.R '${qc_passed_sces.join(',')}' SCE_merged.rds
    """
}

process JOINT_CLUSTERING {
    publishDir "${params.outdir}/SCE", mode: 'copy'

    input:
    path "clustering.Rmd"
    path "sc_functions.R"
    path "samples.csv"
    path "marker_genes.csv"
    path merged_sce
    val seed

    output:
    path "clustering.nb.html", emit: report
    path "SCE_clustered.rds", emit: tumor_annotated_sce

    script:
    """
    Rscript -e "rmarkdown::render('clustering.Rmd',
                output_file = 'clustering.nb.html',
                output_format = 'html_notebook',
                output_options = list(
                    self_contained = TRUE,
                    df_print = 'paged',
                    code_folding = 'hide',
                    toc = TRUE,
                    toc_float = TRUE
                ),
                params = list(
                    seed = ${seed},
                    input_file = '${merged_sce}',
                    output_file = 'SCE_clustered.rds',
                    sc_functions = 'sc_functions.R',
                    marker_genes = 'marker_genes.csv',
                    metadata_file = 'samples.csv'
                    ))"
    """
}