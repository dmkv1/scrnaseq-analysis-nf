#!/usr/bin/env Rscript
options(scipen = 100)

# path to output folder
workdir <- "/mnt/data/NGS/Projects/MCL-scrnaseq/MCL-PhanthomMenace/CNV_inference/inferCNV"

gencode_ref <- file.path(workdir, "hg38_gencode_v27.txt")

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(scran)
  library(scater)
  library(infercnv)
})

k_obs_groups.list <- list(
  "P009" = 5,
  "P022" = 2,
  "P027" = 3,
  "P087" = 6,
  "P069" = 3
)

sce <- readRDS(file = "/mnt/data/NGS/Projects/MCL-scrnaseq/MCL-PhanthomMenace/scrnaseq-analysis-nxf/results/SCE/annotated_tumor.sce")

sce$cell_type_infercnv <- colData(sce) %>%
  as.data.frame() %>%
  mutate(cell_type_manual = case_when(
    cell_type_manual == "MCL cells" ~ paste("MCL", Timepoint, Compartment),
    TRUE ~ "Reference"
  )) %>% pull(cell_type_manual)

patient_i = "P009"
for (patient_i in names(k_obs_groups.list)) {
  sce_patient <- sce[, sce$Patient == patient_i]
  
  sample_annotations <- data.frame(Barcode = colnames(sce_patient),
                                   Cell_type = sce_patient$cell_type_infercnv)
  
  inferCNV.sample.path <- file.path(workdir, patient_i)
  dir.create(inferCNV.sample.path,
             showWarnings = F,
             recursive = T)
  
  annotations_filename <-
    paste0(patient_i, ".sampleAnnotation.tsv")
  write_tsv(
    sample_annotations,
    file = file.path(inferCNV.sample.path, annotations_filename),
    col_names = FALSE,
    quote = "none"
  )
  
  inferCNVobj <- CreateInfercnvObject(
    raw_counts_matrix = as.matrix(counts(sce_patient)),
    annotations_file = file.path(inferCNV.sample.path, annotations_filename),
    gene_order_file = gencode_ref,
    ref_group_names = "Reference"
  )
  
  inferCNVobj <- infercnv::run(
    inferCNVobj,
    cutoff = 0.1,
    out_dir = inferCNV.sample.path,
    resume_mode = TRUE,
    
    cluster_by_groups = FALSE,
    cluster_references = FALSE,
    denoise = TRUE,
    
    k_obs_groups = k_obs_groups.list[[patient_i]],
    
    HMM = TRUE,
    HMM_type = 'i6',
    analysis_mode = "subclusters",
    hclust_method = 'ward.D2',
    
    tumor_subcluster_partition_method = "random_trees",
    tumor_subcluster_pval = 0.01,
    
    write_phylo = TRUE,
    
    no_prelim_plot = TRUE,
    output_format = "png",
    png_res = 600,
    
    num_threads = 16
  )
  
  print(paste("Sample", patient_i, "done!\n"))
}

for (patient_i in names(k_obs_groups.list)) {
  inferCNV.sample.path <-
    file.path(workdir, patient_i)
  file.copy(
    file.path(inferCNV.sample.path, "infercnv.png"),
    file.path(workdir, paste0(patient_i, ".png"))
  )
}
