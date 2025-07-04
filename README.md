# Analysis overview

- Setup
    1. Docker env - `container`.
    2. Reference - `reference`, not used.
- scRNAseq part
    1. Alignment - `scrnaseq-alignment`, performed in STARsolo.
    2. QC, filtering, tumor cell annotation - `scrnaseq-analysis-nxf`, canned in nextflow. Result is SingleCellExperiment.
    3. CNV detection - `CNV_inference/inferCNV`, random trees clustering.
    4. Sub-clone comparison - `DifferentialExpression`, effect size comparison between DGs and DGex -> ranked gene list, GSEA.
    5. Normal cell typa annotation - `CellTypeAnnotation_final`, using effect size markers to annotate cells, separately in blood and in tissues.
- WES part
    1. Alignment and SNV detection - `WES-snakemake`, latest main branch
    2. CNV detection - `WES-CNV-nf`, WIP to can CNVkit bash scripts
- Other
    1. `tests` - tried benchmarking different decontamination methods. Turns out decontX is the most practical.

# Alignment

STAR version 2.7.10b
Reference genome by 10X genomics, GRCh38-2020-A
key arguments from `scrnaseq-alignment/STARsolo.sh` mapping script:

```bash
STAR \
        --outFilterType BySJout \
        --limitBAMsortRAM 128000000000 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes NH HI AS nM CR CY UR UY GX GN CB UB \
        --outFilterScoreMin 30 \
        --outFilterMultimapNmax 10 \
        --outSAMmultNmax 1 \
        --soloMultiMappers EM \
        --outMultimapperOrder Random \
        --clipAdapterType CellRanger4 \
        --soloType CB_UMI_Simple \
        --soloStrand Forward \
        --soloCBstart 1 \
        --soloCBlen 16 \
        --soloUMIstart 17 \
        --soloUMIlen 10 \
        --soloBarcodeReadLength 0 \
        --soloUMIdedup 1MM_CR \
        --soloUMIfiltering MultiGeneUMI_CR \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --soloCellFilter None \
        --soloFeatures Gene GeneFull
```

Only uniquiely mapped transcripts were used downstream.

## HVG
Seurat with ~2000 HVG - for strong separation by the broad traits.
ModelGeneVar with 0.1 FDR - for revealing intricate structure based on specific co-expression patterns.