library(tidyverse)

wd <- "/mnt/data/NGS/Projects/MCL-scrnaseq/MCL-PhanthomMenace/scrnaseq-alignment"

fastq.dir <- "/mnt/data/NGS/Projects/MCL-scrnaseq/data/scRNAseq/fastq"

samples <- data.frame(
  files = list.files(fastq.dir, pattern = ".fastq.gz$", full.names = T)
) %>% 
  mutate(
    file.name = basename(files),
    sample_read = str_remove(file.name, ".fastq.gz"),
    read = word(sample_read, -1, sep = "_"),
    sample = str_remove(sample_read, "_R[1,2]"),
    sample_read = NULL,
    files = NULL
  ) %>% 
  pivot_wider(values_from = file.name, names_from = read)

write_tsv(samples, file.path(wd, "samples"), col_names = FALSE)
