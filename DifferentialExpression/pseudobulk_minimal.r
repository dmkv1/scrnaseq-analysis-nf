###
# Essential version of `pseudobulk.R` script
# Edited by D. Manakov

library(muscat)
library(SingleCellExperiment)
library(DESeq2)
library(edgeR)
library(limma)
library(ROTS)

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
    mm_highcells <- model.matrix( ~ 0 + factor(as.character(group)))
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
    mm_highcells <- model.matrix( ~ 0 + factor(as.character(group)))
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
    mm_highcells <- model.matrix( ~ 0 + factor(as.character(group)))
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
