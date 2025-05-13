process DROPLETS_TO_CELLS {
    tag "${sample_id}"
    publishDir "${params.outdir}/QC/${sample_id}/1_droplets", mode: 'copy'

    input:
    path '1_droplets_to_cells.Rmd'
    tuple val(sample_id), val(expected_cells), val(patient_id), val(timepoint), val(compartment), val(replicate), path(counts)
    val seed

    output:
    tuple val(sample_id), path("${sample_id}_droplets_to_cells.nb.html"), emit: report
    tuple val(sample_id), path("${sample_id}_cells.sce"), emit: sce
    tuple val(sample_id), path("${sample_id}_droplet_metrics.json"), emit: metrics

    script:
    """
    Rscript -e "rmarkdown::render('1_droplets_to_cells.Rmd',
                output_file = '${sample_id}_droplets_to_cells.nb.html',
                output_format = 'html_notebook',
                output_options = list(code_folding = 'hide',
                                    toc = TRUE,
                                    toc_float = TRUE),
                params = list(
                               sample_id = '${sample_id}',
                               expected_cells = '${expected_cells}',
                               patient_id = '${patient_id}',
                               timepoint  = '${timepoint}',
                               compartment = '${compartment}',
                               replicate = '${replicate}',
                               counts_dir = '${counts}',
                               FDR_thresh = '${params.droplets.FDR_thresh}',
                               path_sce_output = '${sample_id}_cells.sce',
                               seed = ${seed}
                             ))"
    """
}

process DOUBLET_DETECTION {
    tag "${sample_id}"
    publishDir "${params.outdir}/QC/${sample_id}/2_doublets", mode: 'copy'

    input:
    path '2_doublets.Rmd'
    tuple val(sample_id), path(sce)
    val seed

    output:
    tuple val(sample_id), path("${sample_id}_doublets.nb.html"), emit: report
    tuple val(sample_id), path("${sample_id}_singlets.sce"), emit: sce
    tuple val(sample_id), path("${sample_id}_doublet_metrics.json"), emit: metrics

    script:
    """
    Rscript -e "rmarkdown::render('2_doublets.Rmd',
                output_file = '${sample_id}_doublets.nb.html',
                output_format = 'html_notebook',
                output_options = list(code_folding = 'hide',
                                    toc = TRUE,
                                    toc_float = TRUE),
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
    publishDir "${params.outdir}/QC/${sample_id}/3_cellQC", mode: 'copy'

    input:
    path '3_cell_qc.Rmd'
    tuple val(sample_id), path(sce)
    val seed

    output:
    tuple val(sample_id), path("${sample_id}_cell_QC.nb.html"), emit: report
    tuple val(sample_id), path("${sample_id}_cells.sce"), emit: sce
    tuple val(sample_id), path("${sample_id}_cell_QC_metrics.json"), emit: metrics

    script:
    """
    Rscript -e "rmarkdown::render('3_cell_qc.Rmd',
                output_file = '${sample_id}_cell_QC.nb.html',
                output_format = 'html_notebook',
                output_options = list(code_folding = 'hide',
                                    toc = TRUE,
                                    toc_float = TRUE),
                params = list(
                    sample_id = '${sample_id}',
                    path_sce_input = '${sce}',
                    path_sce_output = '${sample_id}_cells.sce',
                    nUMI_thresh = '${params.qc.nUMI_thresh}',
                    nGenes_thresh = '${params.qc.nGenes_thresh}',
                    mitochondrial_thresh = '${params.qc.mitochondrial_thresh}',
                    cluster_discard_thresh = '${params.qc.cluster_discard_thresh}',
                    seed = ${seed}
                    ))"
    """
}
