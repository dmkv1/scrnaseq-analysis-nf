process DROPLETS_TO_CELLS {
    tag "${sample_id}"
    publishDir "${params.outdir}/QC/${sample_id}/1_droplets", mode: 'copy'

    input:
    path 'droplets_to_cells.Rmd'
    tuple val(sample_id), path(sce), val(expected_cells)
    val seed

    output:
    tuple val(sample_id), path("${sample_id}_droplets_to_cells.nb.html"), emit: report
    tuple val(sample_id), path("${sample_id}_cells.sce"), emit: sce
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
                               path_sce_output = '${sample_id}_cells.sce',
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
    tuple val(sample_id), path("${sample_id}_singlets.sce"), emit: sce
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
                               path_sce_output = '${sample_id}_singlets.sce',
                               seed = ${seed}
                             ))"
    """
}

process CELL_QC {
    tag "${sample_id}"
    publishDir "${params.outdir}/QC/${sample_id}/4_cellQC", mode: 'copy'

    input:
    path 'cell_qc.Rmd'
    path 'process_sce.R'
    tuple val(sample_id), path(sce), val(nUMI_thresh), val(nGenes_thresh), val(mito_thresh), val(contam_thresh), path(perCellCont)
    val seed

    output:
    tuple val(sample_id), path("${sample_id}_cell_QC.nb.html"), emit: report
    path("${sample_id}.sce"), emit: sce
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
                    path_sce_output = '${sample_id}.sce',
                    path_perCellCont = '${sample_id}_perCellCont.rds',
                    nUMI_thresh = '${nUMI_thresh}',
                    nGenes_thresh = '${nGenes_thresh}',
                    mito_thresh = '${mito_thresh}',
                    contam_thresh = '${contam_thresh}',
                    cluster_discard_thresh = '${params.qc.cluster_discard_thresh}',
                    process_fun = 'process_sce.R',
                    seed = ${seed}
                    ))"
    """
}
