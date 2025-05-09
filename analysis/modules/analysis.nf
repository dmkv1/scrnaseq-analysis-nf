process SEURAT_CLUSTERING {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    path('seurat_clustering.Rmd')
    tuple val(sample_id), path(sce)
    val seed

    output:
    path("${sample_id}_seurat_clustering.nb.html"), emit: report
    path("${sample_id}.sce"), emit: sce
    path("${sample_id}_markers.rds"), emit: markers

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

process MERGE {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path qc_passed_sces
    val seed

    output:
    path 'merged.sce', emit: merged_sce

    script:
    """
    merge.R '${qc_passed_sces.join(',')}' merged.sce
    """
}

process ANNOTATE {
    publishDir "${params.outdir}/annotation", mode: 'copy'

    input:
    path('5_annotation.Rmd')
    path merged_sce
    val seed

    output:
    path("annotation.nb.html"), emit: report
    path("annotated.sce"), emit: annotated_sce

    script:
    """
    Rscript -e "rmarkdown::render('5_annotation.Rmd',
                output_file = 'annotation.nb.html',
                output_format = 'html_notebook',
                output_options = list(code_folding = 'hide',
                                    toc = TRUE,
                                    toc_float = TRUE),
                params = list(
                    input_file = '${merged_sce}',
                    output_file = 'annotated.sce',
                    seed = ${seed}
                    ))"
    """
}