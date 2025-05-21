#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(optparse)
  library(celda)
})

option_list <- list(
  make_option(c("--cells_fbmtx"), type = "character", default = NULL),
  make_option(c("--droplets_fbmtx"), type = "character", default = NULL),
  make_option(c("--use_empty"), type = "logical", default = FALSE),
  make_option(c("--decont_fbmtx"), type = "character", default = "decontX_feature_bc_matrix"),
  make_option(c("--perCell"), type = "character", default = "perCell.rds")
)
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

run_DecontX <- function(cells_fbmtx,
                        droplets_fbmtx,
                        use_empty,
                        decont_fbmtx,
                        perCell) {
  counts_cells <- Read10X(cells_fbmtx)
  
  z = NULL
  
  # Background profile
  if (use_empty) {
    print("Using empty droplet profile")
    raw <- Read10X(droplets_fbmtx)
  } else {
    print(paste0("use_empty = ", use_empty, ": empty droplets not used"))
    raw <- NULL
  }
  
  res <- decontX(x = counts_cells,
                 z = z,
                 background = raw)
  
  perCell_cont <- data.frame(cell = colnames(counts_cells),
                             cont = res$contamination)
  saveRDS(perCell_cont, perCell)
  
  corcounts <- ceiling(res$decontXcounts)
  dir.create(decont_fbmtx, recursive = TRUE)
  
  # 1. Copy barcodes
  file.copy(
    file.path(cells_fbmtx, "barcodes.tsv.gz"),
    file.path(decont_fbmtx, "barcodes.tsv.gz")
  )
  # 2. Copy genes
  file.copy(
    file.path(cells_fbmtx, "features.tsv.gz"),
    file.path(decont_fbmtx, "features.tsv.gz")
  )
  # 3. Extract and save counts
  Matrix::writeMM(corcounts, file.path(decont_fbmtx, "matrix.mtx"))
  # gzip counts
  system(paste0("gzip ", file.path(decont_fbmtx, "matrix.mtx")))
}

run_DecontX(
  cells_fbmtx = opts$cells_fbmtx,
  droplets_fbmtx = opts$droplets_fbmtx,
  use_empty = opts$use_empty,
  decont_fbmtx = opts$decont_fbmtx,
  perCell = opts$perCell
)