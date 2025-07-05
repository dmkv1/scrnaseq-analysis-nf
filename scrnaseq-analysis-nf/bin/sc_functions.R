# Reduce SCE dimentions and find clusters
process_sce <- function(sce,
                        rowdata = "Symbol",
                        HVF_method = "scran_mgv",
                        n_HVFs = 2000,
                        HVF_fdr_threshold = 0.1,
                        exclude_HVFs = NULL,
                        HVFs = NULL,
                        vars_to_regress = NULL,
                        RunPCA_npcs = 40,
                        pca_variance_threshold = 2,
                        chosen_pcs = NULL,
                        UMAP_min_dist = 0.3,
                        UMAP_n_neighbors = 30,
                        FindNeighbors_knn = 20,
                        FindClusters_res = 0.5,
                        seed = 42,
                        return_obj = "sce") {
  suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
    library(tidyverse)
    library(scran)
  })
  
  if (!(return_obj %in% c("sce", "seurat"))) {
    stop("process_sce supports only these types of return_obj: sce, seurat")
  }
  
  if (!(HVF_method %in% c("seurat_vst", "scran_mgv"))) {
    stop("process_sce supports only these types of HVF_method: seurat_vst, scran_mgv")
  }
  
  set.seed(seed)
  
  rownames(sce) <- rowData(sce)[[rowdata]]
  rowData(sce)[["Symbol"]] <- rownames(sce)
  
  reducedDims(sce)[["PCA"]] <- NULL
  reducedDims(sce)[["UMAP"]] <- NULL
  
  seurat <- as.Seurat(sce, counts = "counts", data = NULL)
  seurat <- NormalizeData(seurat,
                          normalization.method = "LogNormalize",
                          scale.factor = 10000)
  assay(sce, "logcounts", withDimnames = FALSE) <- GetAssayData(seurat, layer = "data", assay = "originalexp")
  
  if (HVF_method == "seurat_vst" & is.null(HVFs)) {
    seurat <- FindVariableFeatures(
      seurat,
      selection.method = "vst",
      loess.span = 0.3,
      nfeatures = n_HVFs,
      verbose = F
    )
    top_HVF <- VariableFeatures(seurat)
    HVF_plot <- LabelPoints(
      plot = VariableFeaturePlot(seurat),
      points = top_HVF,
      repel = TRUE,
      xnudge = 0,
      ynudge = 0
    )
  } else if (HVF_method == "scran_mgv" & is.null(HVFs)) {
    variance_decomposition <- modelGeneVar(sce)
    top_HVF <- getTopHVGs(variance_decomposition, fdr.threshold = HVF_fdr_threshold)
    VariableFeatures(seurat) <- top_HVF
    
    variance_df <- as.data.frame(variance_decomposition) %>%
      rownames_to_column("Symbol") %>%
      left_join(., as.data.frame(rowData(sce)), by = "Symbol") %>%
      mutate(HVF = case_when(Symbol %in% top_HVF ~ TRUE, TRUE ~ FALSE))
    
    HVF_plot <- ggplot(variance_df, aes(x = mean, y = bio, color = HVF)) +
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
          length(top_HVF),
          " variable genes at FDR < ",
          HVF_fdr_threshold
        ),
        x = "Mean of log-expression",
        y = "Variance of log-expression",
        color = "Biological\nvariance"
      )
  } else if (!is.null(HVFs)) {
    HVFs <- HVFs[HVFs %in% rowData(sce)[["Symbol"]]]
    if (length(HVFs) == 0) {
      stop("No HVFs found in the features of the object!")
    }
    VariableFeatures(seurat) <- HVFs
    HVF_plot <- NULL
  }
  
  if (!is.null(exclude_HVFs)) {
    top_HVF <- top_HVF[!top_HVF %in% exclude_HVFs]
    VariableFeatures(seurat) <- top_HVF
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
  if (is.null(chosen_pcs)) {
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
      title = paste0(
        "Chosen PCs: ",
        chosen_pcs,
        " (var threshold: ",
        pca_variance_threshold,
        "%)"
      ),
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
  
  # Remove previous clusters
  cluster_cols <- grep("^clusters_res", colnames(sce@colData), value = TRUE)
  colData(sce) <- colData(sce)[, !(colnames(colData(sce)) %in% cluster_cols)]
  
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
  
  print(HVF_plot)
  print(PCA_variance_plot)
  
  if (return_obj == "sce") {
    return(sce)
  } else if (return_obj == "seurat") {
    return(seurat)
  }
}

# Standardized marker detection based on the effect sizes
effectSizeMarkers <- function(sce,
                              group,
                              AUC_up_thresh = 0.7,
                              AUC_down_thresh = 0.3,
                              logFC.cohen_thresh = 0.5,
                              group_id = "cluster") {
  suppressPackageStartupMessages({
    library(tidyverse)
    library(SingleCellExperiment)
    library(scater)
    library(scran)
  })
  
  groups <- sce[[group]]
  markers <- scoreMarkers(sce, groups, full.stats = F)
  
  markers_sig <- lapply(markers, function(x)
    x %>%
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      dplyr::select(
        gene,
        "self.average",
        "other.average",
        "mean.logFC.cohen",
        "mean.AUC",
        "self.detected"
      ) %>%
      mutate(
        change = case_when(
          mean.AUC > AUC_up_thresh &
            mean.logFC.cohen > logFC.cohen_thresh ~ "up",
          mean.AUC < AUC_down_thresh &
            mean.logFC.cohen < logFC.cohen_thresh ~ "down",
          TRUE ~ "ns",
        ),
        change = factor(change, levels = c("up", "down", "ns"))
      ) %>%
      filter(change != "ns") %>%
      arrange(change, desc(abs(mean.logFC.cohen)))) %>%
    bind_rows(., .id = group_id)
  return(markers_sig)
}

# Enhanced dot plots
plotDots_twotailed <- function(sce,
                               group,
                               high = c(),
                               low = c(),
                               zlim = NULL,
                               title = NULL,
                               column_names_angle = 45,
                               scale = TRUE,
                               center = TRUE
                               ) {
  suppressPackageStartupMessages({
    library(tidyverse)
    library(SingleCellExperiment)
    library(scater)
    library(scran)
  })
  
  marker_genes <- list("high" = high, "low" = low)
  
  all_markers_raw <- unlist(marker_genes)
  all_markers <- unique(all_markers_raw) %>% rev()
  
  # Check if unique() removed any duplicates
  if (length(all_markers) < length(all_markers_raw)) {
    stop("Duplicated markers found in marker_genes lists")
  }
  
  if (length(marker_genes$high) > 0 &
      length(marker_genes$low) > 0) {
    line_position <- length(marker_genes$low) + 0.5
  } else {
    line_position <- NULL
  }
  
  suppressMessages(
    plot <- plotDots(
      sce,
      features = all_markers,
      scale = scale,
      center = center,
      group = group,
      zlim = zlim
    ) +
      scale_color_gradient2(
        low = "blue",
        mid = "gray90",
        high = "red",
        midpoint = 0,
        name = "Average"
      ) +
      labs(x = NULL, y = NULL, title = title) +
      theme(
        axis.text.x = element_text(
          size = 12,
          angle = column_names_angle,
          hjust = 1
        ),
        axis.text.y = element_text(size = 12)
      ) +
      geom_hline(yintercept = line_position, linetype = "solid")
  )
  
  return(plot)
}

# Turn cluster markers into text stream for LLM
markers_to_text <- function(markers,
                            direction = "up",
                            max_length = Inf) {
  stopifnot(direction %in% c("up", "down"))
  
  clusters <- unique(markers$cluster)
  list_markers <- list()
  for (i in clusters) {
    list_markers[[i]] <- markers %>%
      filter(cluster == i) %>%
      filter(change == direction) %>%
      slice_head(n = max_length) %>%
      pull(gene) %>%
      paste(collapse = ", ") %>%
      paste0("Cluster ", i, ": ", .)
  }
  cat(unlist(list_markers), sep = "\n\n")
}

# Сontingency heatmap between two cell labels
plotContingency <- function(sce, columns = "cell_type_manual", rows = "cell_type_conv", column_names_angle = 45) {
  suppressPackageStartupMessages(library(ComplexHeatmap))
  
  contingency_matrix  <- colData(sce) %>%
    as.data.frame() %>%
    dplyr::select(cols = any_of(columns), rows = any_of(rows)) %>%
    with(table(rows, cols))
  
  color_fun <- circlize::colorRamp2(
    breaks = c(0, 10, 100, 1000, 10000, 30000),
    colors = c("white", "#FFE5E5", "#FFB3B3", "#FF6666", "red", "#FF3366")
  )
  
  ht <- Heatmap(
    contingency_matrix,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(contingency_matrix[i, j], x, y, gp = gpar(fontsize = 10, col = "black"))
    },
    col = color_fun,
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    name = "Count",
    column_names_rot = column_names_angle,
    show_heatmap_legend = FALSE
  )
  
  draw(ht)
}

############ Comparison between conditions

CreatePseudoBulkData <- function(raw.data,
                                 normalized.data,
                                 sample.id,
                                 fun = "sum")
{
  # fun can be "sum" or "mean"
  if (fun == "sum")
  {
    sce <- SingleCellExperiment(assays = list(counts = raw.data),
                                colData = DataFrame(sample_id = sample.id))
  } else if (fun == "mean") {
    sce <- SingleCellExperiment(
      assays = list(counts = normalized.data),
      colData = DataFrame(sample_id = sample.id)
    )
  }
  pb <- muscat::aggregateData(sce,
                              assay = "counts",
                              fun = fun,
                              by = c("sample_id"))
  return(pb)
}

RunPseudobulkMethod <- function(raw.data,
                                normalized.data,
                                individual,
                                group,
                                test = "ROTS",
                                sum.or.mean = "sum")
{
  suppressPackageStartupMessages({
    library(muscat)
    library(SingleCellExperiment)
    library(DESeq2)
    library(edgeR)
    library(limma)
    library(ROTS)
  })
  
  pb <- CreatePseudoBulkData(
    raw.data = raw.data,
    normalized.data = normalized.data,
    sample.id = individual,
    fun = sum.or.mean
  )
  
  if (test == "ROTS")
  {
    group <- apply(as.matrix(table(individual, group)), 1, function(x)
      x[1] == 0)
    group_names <- names(group)
    group <- as.numeric(group)
    names(group) <- group_names
    group <- group[colnames(pb)]
    
    if (sum.or.mean == "sum")
    {
      y <- DGEList(assay(pb), remove.zeros = T)
      y <- calcNormFactors(y)
      logcpm <- edgeR::cpm(
        y,
        normalized.lib.sizes = T,
        prior.count = 1,
        log = T
      )
      resrots <- ROTS(data = logcpm,
                      groups = group,
                      seed = 1234)
    }
    else if (sum.or.mean == "mean")
    {
      y <- assay(pb)
      print(dim(y))
      y <- y[apply(y, 1, sum) != 0, ]
      print(dim(y))
      resrots <- ROTS(data = y,
                      groups = group,
                      seed = 1234)
    }
    rotsdf <- as.data.frame(resrots$logfc)
    rotsdf$gene <- rownames(rotsdf)
    rotsdf <- cbind(rotsdf, resrots$pvalue)
    rotsdf <- cbind(rotsdf, resrots$FDR)
    rownames(rotsdf) <- 1:nrow(rotsdf)
    colnames(rotsdf) <- c("logFC", "gene", "pvalue", "FDR")
    rotsdf <- rotsdf[, c("logFC", "pvalue", "FDR", "gene")]
    colnames(rotsdf) <- c("logFC", "pvalue", "padj", "gene")
    return(rotsdf)
  }
  
  else if (test == "Limma")
  {
    group <- apply(as.matrix(table(individual, group)), 1, function(x)
      x[1] == 0)
    group_names <- names(group)
    group <- as.numeric(group)
    names(group) <- group_names
    group <- group[colnames(pb)]
    group <- plyr::mapvalues(group, from = c(0, 1), to = c("A", "B"))
    
    
    #Make the model and contrast for statistical testing
    mm_highcells <- model.matrix(~ 0 + factor(as.character(group)))
    dimnames(mm_highcells) <- list(names(group), levels(factor(as.character(group))))
    mm_highcells <- mm_highcells[colnames(pb), ]
    contrast_highcells <- makeContrasts("B-A", levels = mm_highcells)
    
    if (sum.or.mean == "sum")
    {
      dge <- DGEList(counts = assay(pb), remove.zeros = T)
      dge <- calcNormFactors(dge)
      v <- voom(dge, design = mm_highcells, plot = F)
      fit <- lmFit(v, design = mm_highcells)
      fit2 <- contrasts.fit(fit, contrasts = contrast_highcells)
      fit2 <- eBayes(fit2)
      limma_out <- topTable(fit2, number = nrow(dge))
      
    } else {
      v <- assay(pb)
      v <- v[apply(v, 1, sum) != 0, ]
      fit <- lmFit(v, design = mm_highcells)
      fit2 <- contrasts.fit(fit, contrasts = contrast_highcells)
      fit2 <- eBayes(fit2)
      limma_out <- topTable(fit2, number = nrow(v))
      
    }
    df <- limma_out[, c("logFC", "P.Value", "adj.P.Val")]
    df$gene <- rownames(limma_out)
    colnames(df) <- c("logFC", "pvalue", "padj", "gene")
    
    return(df)
  }
  
  else if (test == "edgeR")
  {
    group <- apply(as.matrix(table(individual, group)), 1, function(x)
      x[1] == 0)
    group_names <- names(group)
    group <- as.numeric(group)
    names(group) <- group_names
    group <- group[colnames(pb)]
    group <- plyr::mapvalues(group, from = c(0, 1), to = c("A", "B"))
    
    
    #Make the model and contrast for statistical testing
    mm_highcells <- model.matrix(~ 0 + factor(as.character(group)))
    dimnames(mm_highcells) <- list(names(group), levels(factor(as.character(group))))
    mm_highcells <- mm_highcells[colnames(pb), ]
    contrast_highcells <- makeContrasts("B-A", levels = mm_highcells)
    
    
    dge <- DGEList(counts = assay(pb),
                   group = group,
                   remove.zeros = T)
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, design = mm_highcells)
    fit <- glmQLFit(dge, design = mm_highcells)
    fit2 <- glmQLFTest(fit, contrast = contrast_highcells)
    tt <- topTags(fit2, n = nrow(dge))
    edger_out <- tt$table
    
    df <- edger_out[, c("logFC", "PValue", "FDR")]
    df$gene <- rownames(edger_out)
    colnames(df) <- c("logFC", "pvalue", "padj", "gene")
    
    
    return(df)
  }
  
  else if (test == "DESeq2")
  {
    group <- apply(as.matrix(table(individual, group)), 1, function(x)
      x[1] == 0)
    group_names <- names(group)
    group <- as.numeric(group)
    names(group) <- group_names
    group <- group[colnames(pb)]
    group <- plyr::mapvalues(group, from = c(0, 1), to = c("A", "B"))
    
    #Make the model and contrast for statistical testing
    mm_highcells <- model.matrix(~ 0 + factor(as.character(group)))
    dimnames(mm_highcells) <- list(names(group), levels(factor(as.character(group))))
    mm_highcells <- mm_highcells[colnames(pb), ]
    contrast_highcells <- makeContrasts("B-A", levels = mm_highcells)
    
    pb_counts <- assay(pb)
    mode(pb_counts) <- "integer"
    dge <- DESeqDataSetFromMatrix(pb_counts, colData = colData(pb), design = mm_highcells)
    dds <- DESeq(dge)
    res <- results(dds, contrast = contrast_highcells)
    deseq_out <- as.data.frame(res@listData, row.names = res@rownames)
    
    df <- deseq_out[, c("log2FoldChange", "pvalue", "padj")]
    df$gene <- rownames(deseq_out)
    colnames(df) <- c("logFC", "pvalue", "padj", "gene")
    return(df)
  }
}

# tableau colors for plots (used by scater)
get_palettes <- function (palette_name)
{
  switch(
    palette_name,
    tableau20 = c(
      "#1F77B4",
      "#AEC7E8",
      "#FF7F0E",
      "#FFBB78",
      "#2CA02C",
      "#98DF8A",
      "#D62728",
      "#FF9896",
      "#9467BD",
      "#C5B0D5",
      "#8C564B",
      "#C49C94",
      "#E377C2",
      "#F7B6D2",
      "#7F7F7F",
      "#C7C7C7",
      "#BCBD22",
      "#DBDB8D",
      "#17BECF",
      "#9EDAE5"
    ),
    tableau10medium = c(
      "#729ECE",
      "#FF9E4A",
      "#67BF5C",
      "#ED665D",
      "#AD8BC9",
      "#A8786E",
      "#ED97CA",
      "#A2A2A2",
      "#CDCC5D",
      "#6DCCDA"
    ),
    colorblind10 = c(
      "#006BA4",
      "#FF800E",
      "#ABABAB",
      "#595959",
      "#5F9ED1",
      "#C85200",
      "#898989",
      "#A2C8EC",
      "#FFBC79",
      "#CFCFCF"
    ),
    colourblind10 = c(
      "#006BA4",
      "#FF800E",
      "#ABABAB",
      "#595959",
      "#5F9ED1",
      "#C85200",
      "#898989",
      "#A2C8EC",
      "#FFBC79",
      "#CFCFCF"
    ),
    trafficlight = c(
      "#B10318",
      "#DBA13A",
      "#309343",
      "#D82526",
      "#FFC156",
      "#69B764",
      "#F26C64",
      "#FFDD71",
      "#9FCD99"
    ),
    purplegray12 = c(
      "#7B66D2",
      "#A699E8",
      "#DC5FBD",
      "#FFC0DA",
      "#5F5A41",
      "#B4B19B",
      "#995688",
      "#D898BA",
      "#AB6AD5",
      "#D098EE",
      "#8B7C6E",
      "#DBD4C5"
    ),
    bluered12 = c(
      "#2C69B0",
      "#B5C8E2",
      "#F02720",
      "#FFB6B0",
      "#AC613C",
      "#E9C39B",
      "#6BA3D6",
      "#B5DFFD",
      "#AC8763",
      "#DDC9B4",
      "#BD0A36",
      "#F4737A"
    ),
    greenorange12 = c(
      "#32A251",
      "#ACD98D",
      "#FF7F0F",
      "#FFB977",
      "#3CB7CC",
      "#98D9E4",
      "#B85A0D",
      "#FFD94A",
      "#39737C",
      "#86B4A9",
      "#82853B",
      "#CCC94D"
    ),
    cyclic = c(
      "#1F83B4",
      "#1696AC",
      "#18A188",
      "#29A03C",
      "#54A338",
      "#82A93F",
      "#ADB828",
      "#D8BD35",
      "#FFBD4C",
      "#FFB022",
      "#FF9C0E",
      "#FF810E",
      "#E75727",
      "#D23E4E",
      "#C94D8C",
      "#C04AA7",
      "#B446B3",
      "#9658B1",
      "#8061B4",
      "#6F63BB"
    )
  )
}
