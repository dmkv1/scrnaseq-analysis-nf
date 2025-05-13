#!/usr/bin/env Rscript
options(scipen = 100)

suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
  library(tidyverse)
  library(infercnv)
})

option_list <- list(
  make_option(
    c("-i", "--input_sce_file"),
    type = "character",
    default = NULL,
    help = "Input SingleCellExperiment object stored as RDS"
  ),
  make_option(
    c("--patient"),
    type = "character",
    default = NULL,
    help = "Patient ID"
  ),
  make_option(
    c("--tumor_annotation"),
    type = "character",
    default = "Tumor cell",
    help = "How the tumor cells are annotated in the SCE object"
  ),
  make_option(
    c("--tumor_label"),
    type = "character",
    default = "Tumor",
    help = "How to label the tumor cells in the inferCNV object and plots"
  ),
  make_option(
    c("--k_obs_groups"),
    type = "integer",
    default = 4,
    help = "k_obs_groups parameter"
  ),
  make_option(
    c("--hg38_gencode"),
    type = "character",
    default = NULL,
    help = "hg38 gencode reference"
  ),
  make_option(
    c("--seed"),
    type = "integer",
    default = 42,
    help = "Random seed"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

seed <- opts$seed
set.seed(seed)

input_sce_file <- opts$input_sce_file
patient <- opts$patient

tumor_annotation <- opts$tumor_annotation
tumor_label <- opts$tumor_label
reference_label <- "Reference"

k_obs_groups <- opts$k_obs_groups
hg38_gencode <- opts$hg38_gencode

# Setup SCE
sce <- readRDS(input_sce_file)
sample_sce <- sce[, sce$Patient == patient]

sample_sce$infercnv_groups <- colData(sample_sce) %>%
  as.data.frame() %>%
  mutate(
    cell_type_manual = case_when(
      cell_type_manual == tumor_annotation ~ paste(tumor_label, Timepoint, Compartment, sep = "_"),
      TRUE ~ reference_label
    )
  ) %>%
  pull(cell_type_manual)

# Construct inferCNV object
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
  
  k_obs_groups = k_obs_groups,
  
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
