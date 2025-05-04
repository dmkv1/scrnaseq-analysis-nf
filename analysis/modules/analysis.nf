process SEURAT_CLUSTERING {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    path('seurat_clustering.Rmd')
    tuple val(sample_id), path(sce)
    val seed

    output:
    tuple val(sample_id), path("${sample_id}_seurat_clustering.nb.html"), emit: report
    tuple val(sample_id), path("${sample_id}.sce"), emit: sce
    tuple val(sample_id), path("${sample_id}_markers.rds"), emit: markers

    script:
    """
    Rscript -e "rmarkdown::render('seurat_clustering.Rmd',
                output_file = '${sample_id}_seurat_clustering.nb.html',
                output_format = 'html_notebook',
                output_options = list(code_folding = 'hide',
                                    toc = TRUE,
                                    toc_float = TRUE),
                params = list(
                    sample_id = '${sample_id}',
                    path_sce_input = '${sce}',
                    path_sce_output = '${sample_id}.sce',
                    path_markers_output = '${sample_id}_markers.rds',
                    FindVariableFeatures_nfeatures = '${params.clustering.FindVariableFeatures_nfeatures}',
                    RunPCA_npcs = '${params.clustering.RunPCA_npcs}',
                    PCA_variance_threshold = '${params.clustering.PCA_variance_threshold}',
                    FindNeighbors_knn = '${params.clustering.FindNeighbors_knn}',
                    FindClusters_res = '${params.clustering.FindClusters_res}',
                    markers_logfc_thresh = '${params.clustering.markers_logfc_thresh}',
                    markers_min_pct = '${params.clustering.markers_min_pct}',
                    markers_test_use = '${params.clustering.markers_test_use}',
                    seed = ${seed}
                    ))"
    """
}

process INTEGRATION {
    publishDir "${params.outdir}/integration", mode: 'copy'

    input:
    path('integration.Rmd')
    val(sample_ids)
    path(sces)
    val seed

    output:
    path("integrated_samples.nb.html"), emit: report
    path("integrated_samples.sce"), emit: sce

    script:
    """
    Rscript -e "rmarkdown::render('5_integration.Rmd',
                output_file = 'integrated_samples.nb.html',
                output_format = 'html_notebook',
                output_options = list(code_folding = 'hide',
                                    toc = TRUE,
                                    toc_float = TRUE),
                params = list(
                    sample_files = 'sces',
                    output_file = 'integrated_samples.sce',
                    seed = ${seed}
                    ))"
    """
}
