#!/usr/bin/env bash

set -xe

SCH=$1
MY_ENV=wf

eval "$(conda shell.bash hook)"
conda activate $MY_ENV
conda list
cd snakemake && ls -alh
snakemake --use-conda --configfile analysis.yaml \
  --latency-wait 60 --jobs \
  --cluster "xenon -vvv scheduler $SCH --location local:// submit \
    --name smk.{rule} --inherit-env --max-run-time 5 --working-directory . \
    --stderr stderr-%j.log --stdout stdout-%j.log"

echo -e "\nLog files:"
ls *.log
for f in *.log; do
  echo -e "\n### $f ###\n"
  cat $f
done
