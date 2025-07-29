#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
  library(DropletUtils)
})

option_list <- list(
  make_option(
    c("--fbmtx_path"),
    type = "character",
    default = NULL,
    help = "Path to a 'feature_bc_matrix' directory"
  ),
  make_option(
    c("--sample_id"),
    type = "character",
    default = NULL,
    help = "Sample ID"
  ),
  make_option(
    c("--patient_id"),
    type = "character",
    default = NULL,
    help = "Patient ID"
  ),
  make_option(
    c("--label"),
    type = "character",
    default = NULL,
    help = "Additional label(s)"
  ),
  make_option(
    c("--replicate"),
    type = "character",
    default = NULL,
    help = "Replicate"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

fbmtx_path <- opts$fbmtx_path
sample_id <- opts$sample_id
patient_id <- opts$patient_id
label <- opts$label
replicate <- opts$replicate

sce <- DropletUtils::read10xCounts(
  fbmtx_path,
  sample.names = sample_id,
  col.names = TRUE
)

# Assign sample column values
sce[["patient"]] <- patient_id
sce[["label"]] <- label
sce[["replicate"]] <- replicate

# Add metadata
metadata(sce)$sample <- sample_id
metadata(sce)$patient <- patient_id
metadata(sce)$label <- label
metadata(sce)$replicate <- replicate

cat("Writing SCE:")
print(sce)

saveRDS(sce, paste0(sample_id, ".sce"))
