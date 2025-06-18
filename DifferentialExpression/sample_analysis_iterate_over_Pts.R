#!/usr/bin/env Rscript

setwd("/mnt/data/NGS/Projects/MCL-scrnaseq/MCL-PhanthomMenace/DifferentialExpression")

patients <- c("P009", "P022", "P027", "P087", "P069")

for(patient_i in patients){
  params_list <- list(
    patient = patient_i
  )
  
  rmarkdown::render(
    "sample_analysis.Rmd",
    output_format = "html_notebook",
    params = params_list,
    output_options = list(
      self_contained = TRUE,
      code_folding = "show",
      toc = FALSE,
      df_print = "paged"
    ),
    output_file = file.path(getwd(), "results", paste0(patient_i, "_diffexpranal_report.html"))
  )
}

rmarkdown::render(
  "sample_analysis_P069_RELgut.Rmd",
  output_format = "html_notebook",
  params = list(
    patient = "P069"
  ),
  output_options = list(
    self_contained = TRUE,
    code_folding = "show",
    toc = FALSE,
    df_print = "paged"
  ),
  output_file = file.path(getwd(), "results", paste0("P069_RELgut", "_diffexpranal_report.html"))
)