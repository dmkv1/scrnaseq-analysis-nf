# Analysis overview

- Setup
    1. Docker env - `container`.
    2. Reference - `reference`, not used.
- scRNAseq part
    1. Alignment - `scrnaseq-alignment`, performed in STARsolo.
    2. QC, filtering, tumor cell annotation - `scrnaseq-analysis-nf`, canned in nextflow. Result is SingleCellExperiment.
    3. CNV detection - `CNV_inference/inferCNV`, random trees clustering.
    4. Sub-clone comparison - `DifferentialExpression`, effect size comparison between DGs and DGex -> ranked gene list, GSEA.
- WES part
    1. Alignment and SNV detection - `WES-snakemake`, latest main branch.
    2. CNV detection - `WES-CNV-nf`, WIP to can CNVkit bash scripts.
- Other
    1. `tests` - tried benchmarking different decontamination methods. Turns out decontX is the most practical.

# Container

Built using the dockerfile, `rebuild.sh` is used for convenience. Should be available at docker: `dmkv/rstudio-server-scrnaseq:1.0`

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

Run `scrnaseq-alignment/zip_outputs.Rmd` - decontX expects zipped feature/barcode/matrix files.

# Process and cluster cells

From the `scrnaseq-analysis-nf` directory, run nextflow:

```bash
nextflow run main.nf
```

This would produce `scrnaseq-analysis-nf/results` where the cleaned and clustered SCE object would be stored.

# Cell type annotation

Cell types were annotated using pre-compiled set of immunological markers `scrnaseq-analysis-nf/marker_genes/Marker_genes.csv`

# CNV inferennce on single cell level

# Differential gene expression between clones

