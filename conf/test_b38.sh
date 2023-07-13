#!/usr/bin/env bash

cd /cbio/users/mamana

nextflow \
  run /users/mamana/federated-imputation/main.nf \
  -c /users/mamana/federated-imputation/conf/test_b38.config \
  -resume -profile slurm,singularity
