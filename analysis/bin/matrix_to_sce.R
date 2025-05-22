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
    c("--timepoint"),
    type = "character",
    default = NULL,
    help = "Timepoint"
  ),
  make_option(
    c("--compartment"),
    type = "character",
    default = NULL,
    help = "Compartment"
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
timepoint <- opts$timepoint
compartment <- opts$compartment
replicate <- opts$replicate

sce <- DropletUtils::read10xCounts(
  fbmtx_path,
  sample.names = sample_id,
  col.names = TRUE
)

# Assign sample column values
sce[["Patient"]] <- patient_id
sce[["Timepoint"]] <- timepoint
sce[["Compartment"]] <- compartment
sce[["Replicate"]] <- replicate

# Add metadata
metadata(sce)$Sample <- sample_id
metadata(sce)$Patient <- patient_id
metadata(sce)$Timepoint <- timepoint
metadata(sce)$Compartment <- compartment
metadata(sce)$Replicate <- replicate

cat("Writing SCE:")
print(sce)

saveRDS(sce, paste0(sample_id, ".sce"))
