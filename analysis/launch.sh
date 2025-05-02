#!/bin/bash
set -e

nextflow run main.nf \
  -resume \
  --input samplesheet.csv \
  --outdir results \
  -bg > nextflow.log 2>&1
