#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SingleCellExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
sce_paths <- args[1]
sce_paths <- strsplit(sce_paths, ",")[[1]]

print(sprintf(
  "Reading %s SingleCellExperiment's",
  length(sce_paths))
)

sce_list <- lapply(sce_paths, readRDS)


for(i in seq_along(sce_list)){
  sce <- sce_list[[i]]

  # Remove row data so the SCEs could be column-bound
  rowData(sce)[["scDblFinder.selected"]] <- NULL
  rowData(sce)[["subsets_Mito"]] <- NULL
  rowData(sce)[["subsets_Ribo"]] <- NULL
  rowData(sce)[["vst.variance.standardized"]] <- NULL
  rowData(sce)[["vst.variable"]] <- NULL
  
  sce_list[[i]] <- sce
  rm(sce)
}

sce <- do.call(cbind, sce_list)

path_sce_output <- args[2]
print(sprintf(
  "Writing %s SingleCellExperiment",
  path_sce_output)
)
saveRDS(sce, path_sce_output)

