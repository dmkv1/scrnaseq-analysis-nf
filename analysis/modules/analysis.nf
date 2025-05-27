process MERGE {
    publishDir "${params.outdir}/SCE", mode: 'copy'

    input:
    path qc_passed_sces

    output:
    path 'merged.sce', emit: merged_sce

    script:
    """
    merge.R '${qc_passed_sces.join(',')}' merged.sce
    """
}

process ANNOTATE {
    publishDir "${params.outdir}/SCE", mode: 'copy'

    input:
    path('annotation.Rmd')
    path 'process_sce.R'
    path merged_sce
    val seed

    output:
    path("annotation.nb.html"), emit: report
    path("annotated.sce"), emit: annotated_sce

    script:
    """
    Rscript -e "rmarkdown::render('annotation.Rmd',
                output_file = 'annotation.nb.html',
                output_format = 'html_notebook',
                output_options = list(code_folding = 'hide',
                                    toc = TRUE,
                                    toc_float = TRUE),
                params = list(
                    seed = ${seed},
                    input_file = '${merged_sce}',
                    output_file = 'annotated.sce',
                    process_fun = 'process_sce.R'
                    ))"
    """
}

process INFERCNV {
    tag "${patient_id}"
    publishDir "${params.outdir}/inferCNV/${patient_id}", mode: 'copy'

    input:
    tuple val(patient_id), val(k_obs_groups)
    path annotated_sce
    path reference_gencode
    val seed

    output:
    path("infercnv.png"), emit: infercnv_plot

    script:
    """
    run_inferCNV.R \
        -i '${annotated_sce}' \
        --patient '${patient_id}' \
        --k_obs_groups '${k_obs_groups}' \
        --hg38_gencode '${reference_gencode}'\
        --seed ${seed}
    """
}