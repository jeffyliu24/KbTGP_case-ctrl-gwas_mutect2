#!/usr/bin/env bash
set -euo pipefail

### Tools
locuszoom=./software/locuszoom-1.4/bin/locuszoom

### Input
region_list=${1:-./data/locus/regions.tsv}
assoc_bed=${2:-./results/gwas/allchr.bed.gz}
ld_vcf_json=${3:-./data/locus/ld_vcf.json}

### Output
outdir=${4:-./results/locuszoom}
mkdir -p ${outdir}

while read -r chr start end; do
  [ -z "${chr}" ] && continue
  ${locuszoom} \
    --epacts ${assoc_bed} \
    --chr ${chr} \
    --start ${start} \
    --end ${end} \
    --ld-vcf ${ld_vcf_json} \
    --build hg38 \
    --ignore-vcf-filter \
    --plotonly \
    --prefix ${outdir}/${chr}_${start}_${end}
done < ${region_list}
