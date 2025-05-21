#!/bin/bash
set -e

nextflow run main.nf \
  --samples metadata_samples.csv \
  --patients metadata_patients.csv \
  --outdir results \
  -resume \
  -bg > nextflow.log 2>&1
