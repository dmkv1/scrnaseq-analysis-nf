#!/bin/bash
set -e

nextflow run main.nf \
  -resume \
  -bg > run.log 2>&1
