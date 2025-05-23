#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(tidyverse)
})

process_sce <- function(sce,
                        HVF_method = "seurat_vst",
                        n_HVGs = 2000,
                        HVG_fdr_threshold = 0.1,
                        exclude_HVGs = NULL,
                        vars_to_regress = NULL,
                        RunPCA_npcs = 40,
                        pca_variance_threshold = 2,
                        chosen_pcs = NULL,
                        UMAP_min_dist = 0.3,
                        UMAP_n_neighbors = 30,
                        FindNeighbors_knn = 20,
                        FindClusters_res = 0.5,
                        seed = 42,
                        return_obj = "sce"
                        ) {
  if (!(return_obj %in% c("sce", "seurat"))) {
    stop("process_sce supports only these types of return_obj: sce, seurat")
  }
  
  if (!(HVF_method %in% c("seurat_vst", "scran_mgv"))) {
    stop("process_sce supports only these types of HVF_method: seurat_vst, scran_mgv")
  }
  
  set.seed(seed)
  
  rownames(sce) <- rowData(sce)[["Symbol"]]
  
  reducedDims(sce)[["PCA"]] <- NULL
  reducedDims(sce)[["UMAP"]] <- NULL
  
  seurat <- as.Seurat(sce, counts = "counts", data = NULL)
  seurat <- NormalizeData(seurat,
                          normalization.method = "LogNormalize",
                          scale.factor = 10000)
  assay(sce, "logcounts", withDimnames = FALSE) <- GetAssayData(seurat, layer = "data", assay = "originalexp")
  
  if (HVF_method == "seurat_vst") {
    seurat <- FindVariableFeatures(
      seurat,
      selection.method = "vst",
      loess.span = 0.3,
      nfeatures = n_HVGs,
      verbose = F
    )
    top_HVG <- VariableFeatures(seurat)
    HVG_plot <- LabelPoints(plot = VariableFeaturePlot(seurat), points = top_HVG,
                            repel = TRUE, xnudge = 0, ynudge = 0)
  } else if (HVF_method == "scran_mgv") {
    variance_decomposition <- modelGeneVar(sce)
    top_HVG <- getTopHVGs(variance_decomposition, fdr.threshold = HVG_fdr_threshold)
    VariableFeatures(seurat) <- top_HVG
    
    variance_df <- as.data.frame(variance_decomposition) %>%
      rownames_to_column("Symbol") %>%
      left_join(., as.data.frame(rowData(sce)), by = "Symbol") %>%
      mutate(HVG = case_when(Symbol %in% top_HVG ~ TRUE, TRUE ~ FALSE))
    
    HVG_plot <- ggplot(variance_df, aes(x = mean, y = bio, color = HVG)) +
      scale_color_manual(values = c("black", "red")) +
      ggrepel::geom_text_repel(
        data = dplyr::filter(variance_df, bio > 0.75),
        aes(label = Symbol),
        min.segment.length = 0,
        color = "gray30"
      ) +
      geom_point(shape = 1) +
      cowplot::theme_cowplot() +
      labs(
        title = "Variance decomposition plot",
        subtitle = paste0(
          "Detected ",
          length(top_HVG),
          " variable genes at FDR = ",
          HVG_fdr_threshold
        ),
        x = "Mean of log-expression",
        y = "Variance of log-expression",
        color = "Biological\nvariance"
      )
  }
  
  if (!is.null(exclude_HVGs)) {
    top_HVG <- top_HVG[!top_HVG %in% exclude_HVGs]
    VariableFeatures(seurat) <- top_HVG
  }
  
  if (is.null(vars_to_regress)) {
    seurat <- ScaleData(seurat, model.use = "linear")
  } else {
    seurat <- ScaleData(seurat, model.use = "linear", vars.to.regress = vars_to_regress)
  }
  
  seurat <- RunPCA(seurat, seed.use = seed, npcs = RunPCA_npcs)
  
  # Extract variance explained by each PC
  pca_variance <- (seurat@reductions$pca@stdev)^2
  percent_var <- 100 * pca_variance / sum(pca_variance)
  
  # Choose PCs that explain more than the threshold percentage of variance
  if(is.null(chosen_pcs)){
    chosen_pcs <- length(which(percent_var > pca_variance_threshold))
  }
  
  # Create a data frame for visualization
  PCvariance_df <- data.frame(
    PC = 1:length(percent_var),
    PercentVariance = percent_var,
    Selected = ifelse(percent_var > pca_variance_threshold, "Yes", "No")
  )
  
  # Plot the variance explained by each PC
  PCA_variance_plot <- ggplot(PCvariance_df, aes(x = PC, y = PercentVariance)) +
    geom_bar(stat = "identity", aes(fill = Selected)) +
    geom_hline(yintercept = pca_variance_threshold,
               linetype = "dashed",
               color = "red") +
    scale_fill_manual(values = c("No" = "gray", "Yes" = "blue")) +
    labs(
      title = paste0("Chosen PCs: ", chosen_pcs, " (var threshold: ", pca_variance_threshold, "%)"),
      x = "Principal Component",
      y = "Percent Variance Explained"
    ) +
    theme_minimal()
  
  seurat <- RunUMAP(
    seurat,
    reduction = "pca",
    dims = 1:chosen_pcs,
    min.dist = UMAP_min_dist,
    n.neighbors = UMAP_n_neighbors
  )
  
  seurat <- FindNeighbors(
    seurat,
    k.param = FindNeighbors_knn,
    reduction = "pca",
    dims = 1:chosen_pcs
  )
  
  seurat <- FindClusters(
    seurat,
    algorithm = 4,
    resolution = FindClusters_res,
    random.seed = seed
  )
  
  # Add data back to SCE
  # Normalization
  assay(sce, "logcounts", withDimnames = FALSE) <- GetAssayData(seurat, layer = "data", assay = "originalexp")
  
  # VST - write HVF data to sce if seurat_vst was used
  rowData_vst <- seurat@assays$originalexp@meta.features
  if (all(rowData_vst[["ID"]] == rowData(sce)[["ID"]])) {
    rowData(sce)[["vst.variance.standardized"]] <- rowData_vst[["vst.variance.standardized"]]
    rowData(sce)[["vst.variable"]] <- rowData_vst[["vst.variable"]]
  } else {
    stop("Variable features do not match!")
  }
  
  # PCA
  reducedDims(sce)[["PCA"]] <- Embeddings(seurat, reduction = "pca")
  attr(reducedDims(sce)[["PCA"]], "rotation") <- Loadings(seurat, reduction = "pca")
  pca_variance <- seurat@reductions$pca@stdev^2
  attr(reducedDims(sce)[["PCA"]], "percentVar") <- pca_variance / sum(pca_variance) * 100
  
  # UMAP
  reducedDims(sce)[["UMAP"]] <- Embeddings(seurat, reduction = "umap")
  
  # clusters
  cluster_colnames <- grep("^originalexp_snn_res", colnames(seurat@meta.data), value = TRUE)
  for (cluster_colname in cluster_colnames) {
    snn_res <- word(cluster_colname, 2, sep = "res\\.")
    snn_res <- paste0("clusters_res_", snn_res)
    
    sce[[snn_res]] <- seurat[[cluster_colname]][[cluster_colname]]
  }
  
  print(HVG_plot)
  print(PCA_variance_plot)
  
  if (return_obj == "sce") {
    return(sce)
  } else if (return_obj == "seurat") {
    return(seurat)
  }
}
