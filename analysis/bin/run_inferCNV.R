#!/usr/bin/env Rscript
options(scipen = 100)

args <- commandArgs(trailingOnly = TRUE)

seed <- args[[1]]
input_sce_file <- args[[2]]
infercnv_output_path <- args[[3]]
hg38_gencode <- args[[4]]

seed <- 42
set.seed(seed)

input_sce_file <- "/mnt/data/NGS/Projects/MCL-scrnaseq/MCL-PhanthomMenace/analysis/results/annotation/annotated.sce"
hg38_gencode <- "/mnt/data/NGS/Projects/MCL-scrnaseq/MCL-PhanthomMenace/analysis/tmp/hg38_gencode_v27.txt"

tumor_cells_label <- "MCL"
reference_label <- "Reference"

k_obs_groups.list <- list(
  "P009" = 5,
  "P022" = 2,
  "P027" = 3,
  "P069" = 3,
  "P087" = 6
)

# Setup SCE
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(tidyverse)
  library(infercnv)
})

sce <- readRDS(input_sce_file)

sce$infercnv_groups <- colData(sce) %>%
  as.data.frame() %>%
  mutate(cell_type_manual = case_when(
    cell_type_manual == "MCL cell" ~ paste(tumor_cells_label, Timepoint, Compartment, sep = "_"),
    TRUE ~ reference_label
  )) %>%
  pull(cell_type_manual)


# Construct inferCNV object grouped by patient
patient = "P009"
sample_sce <- sce[, sce$Patient == patient]


sample_annotations <- data.frame(Barcode = colnames(sample_sce),
                                 Cell_type = sample_sce$infercnv_groups)

annotations_file <- paste0(patient, "_sampleAnnotation.tsv")
write_tsv(
  sample_annotations,
  file = annotations_file,
  col_names = FALSE,
  quote = "none"
)

inferCNVobj <- CreateInfercnvObject(
  raw_counts_matrix = as.matrix(counts(sample_sce)),
  annotations_file = annotations_file,
  gene_order_file = hg38_gencode,
  ref_group_names = reference_label
)

inferCNVobj <- infercnv::run(
  inferCNVobj,
  cutoff = 0.1,
  out_dir = getwd(),
  resume_mode = TRUE,
  
  cluster_by_groups = FALSE,
  cluster_references = FALSE,
  denoise = TRUE,
  
  k_obs_groups = k_obs_groups.list[[patient]],
  
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






