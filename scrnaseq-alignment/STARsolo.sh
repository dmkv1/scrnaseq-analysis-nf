#!/bin/bash
# Increase soft limit for open files which STAR can temporarily require
ulimit -n 65536

WD="/mnt/data/NGS/Projects/Clonal2/scRNAseq"
FASTQ_DIR="$WD/fastq/pooled"
OUTDIR="$WD/alingment/output"

REFDIR="/mnt/data/NGS/refs/STAR_indexed_10X_GRCh38/starsolo-reference"
TRANSCRIPTSGTF="/mnt/data/NGS/refs/10X_GRCh38/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
CBWHITELIST="/mnt/data/NGS/refs/10X_barcodes/3M-february-2018.txt"

while read -r SAMPLE R1_CB_FQ R2_cDNA_FQ; do
    echo -e "Processing sample $SAMPLE\nBarcodes: $R1_CB_FQ\ncDNA: $R2_cDNA_FQ\n"
    
    mkdir -p $OUTDIR/$SAMPLE

    /usr/local/bin/STAR-2.7.10b/bin/Linux_x86_64/STAR \
        --genomeLoad NoSharedMemory \
        --runThreadN 32 \
        --readFilesCommand zcat \
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
        --soloFeatures Gene GeneFull \
        --genomeDir $REFDIR \
        --sjdbGTFfile $TRANSCRIPTSGTF \
        --soloCBwhitelist $CBWHITELIST \
        --readFilesIn $FASTQ_DIR/$R2_cDNA_FQ $FASTQ_DIR/$R1_CB_FQ \
        --outFileNamePrefix "${OUTDIR}/${SAMPLE}/${SAMPLE}_"
        echo ""
done < samples
