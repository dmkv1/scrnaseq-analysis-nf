#!/bin/bash
set -e

nextflow run main.nf \
  -resume \
  --samples metadata_samples.csv \
  --patients metadata_patients.csv \
  --outdir results \
  -bg > nextflow.log 2>&1
