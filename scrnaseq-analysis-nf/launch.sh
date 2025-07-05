#!/bin/bash
set -e

nextflow run main.nf \
  --samples metadata_samples.csv \
  --patients metadata_patients.csv \
  --outdir results \
  -resume \
  -bg > run.log 2>&1
