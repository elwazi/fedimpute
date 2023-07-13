#!/usr/bin/env bash

nextflow \
  run elwazi/fedimpute/main.nf -r main \
  -c conf/test_b38.config \
  -resume -profile slurm,singularity
