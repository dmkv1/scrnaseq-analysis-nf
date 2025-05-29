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

process LABEL_TUMOR_CELLS {
    publishDir "${params.outdir}/SCE", mode: 'copy'

    input:
    path('label_tumor_cells.Rmd')
    path 'process_sce.R'
    path merged_sce
    val seed

    output:
    path("label_tumor_cells.nb.html"), emit: report
    path("annotated_tumor.sce"), emit: tumor_annotated_sce

    script:
    """
    Rscript -e "rmarkdown::render('label_tumor_cells.Rmd',
                output_file = 'label_tumor_cells.nb.html',
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
                    output_file = 'annotated_tumor.sce',
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
        --tumor_annotation 'MCL cell' \
        --tumor_label 'MCL' \
        --num_threads '${params.max_cpus}' \
        --seed ${seed}
    """
}