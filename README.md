# Instalation

- nextflow
- docker

# Inputs

- samples.csv

 # Process and cluster cells

From the `scrnaseq-analysis-nf` directory, run nextflow:

```bash
nextflow run main.nf
```

This would produce `scrnaseq-analysis-nf/results` where the cleaned and clustered SCE object would be stored.
