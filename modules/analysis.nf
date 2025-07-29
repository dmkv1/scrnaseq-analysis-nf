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
