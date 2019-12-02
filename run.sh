#!/bin/bash -xe

source ~/.profile
export BRANCH=dev
git clone -b $BRANCH https://github.com/GooglingTheCancerGenome/sv-gen.git
cd sv-gen/snakemake
snakemake --version
xenon --version

SCH=$1
snakemake --use-conda --configfile analysis.yaml \
  --latency-wait 60 --jobs \
  --cluster "xenon -vvv scheduler $SCH --location local:// submit \
  --name smk.{rule} --cores-per-task 1 --inherit-env \
  --max-run-time 5 --working-directory . \
  --stderr stderr-%j.log --stdout stdout-%j.log"

echo -e "\nLog files:"
ls *.log
for f in $(ls stderr-*.log); do
  echo -e "\n### $f ###\n"
  cat $f
done
