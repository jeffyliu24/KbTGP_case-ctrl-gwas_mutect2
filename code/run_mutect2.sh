#!/usr/bin/env bash
set -euo pipefail

### Tools
gatk=./software/gatk-4.1.9.0/gatk
bcftools=./software/bcftools-1.10.2/bcftools
bgzip=./software/htslib-1.12/bgzip
tabix=./software/htslib-1.12/tabix
annovar=./software/annovar-20200607/annotate_variation.pl

### Reference
reference=./reference/hg38/Homo_sapiens_assembly38.fasta
germline_resource=./reference/hg38/germline_resource.vcf.gz
annovar_db=./database/annovar/humandb

### Input
tumor_id=$1
tumor_bam=$2
out_prefix=$3
whitelist_table=${4:-./data/mutect2/chip_whitelist.tsv}

mkdir -p $(dirname ${out_prefix})
tmpdir=$(dirname ${out_prefix})/tmp
mkdir -p ${tmpdir}

### Mutect2
${gatk} Mutect2 \
  -R ${reference} \
  -I ${tumor_bam} \
  -tumor ${tumor_id} \
  --germline-resource ${germline_resource} \
  --native-pair-hmm-threads 2 \
  -O ${out_prefix}.vcf

${bgzip} -f ${out_prefix}.vcf
${tabix} -f ${out_prefix}.vcf.gz

if [ -f ${out_prefix}.vcf.stats ]; then
  mv ${out_prefix}.vcf.stats ${out_prefix}.vcf.gz.stats
fi

${gatk} FilterMutectCalls \
  -R ${reference} \
  -V ${out_prefix}.vcf.gz \
  -O ${out_prefix}.filtered.vcf.gz \
  --min-allele-fraction 0.02 \
  --unique-alt-read-count 3 \
  --max-events-in-region 2

${tabix} -f ${out_prefix}.filtered.vcf.gz

${bcftools} view -f PASS ${out_prefix}.filtered.vcf.gz \
  -Oz -o ${out_prefix}.pass.vcf.gz
${tabix} -f ${out_prefix}.pass.vcf.gz

${bcftools} view -i 'INFO/DP>30' ${out_prefix}.pass.vcf.gz \
  | ${bcftools} filter -i 'FMT/AD[1]>5 & FMT/AF[0]>0.02' \
  | ${bcftools} view -m2 -M2 -Oz -o ${out_prefix}.candidate.vcf.gz
${tabix} -f ${out_prefix}.candidate.vcf.gz

${bcftools} view -H ${out_prefix}.candidate.vcf.gz \
  | awk 'BEGIN{OFS="\t"} {print $1, $2, $2, $4, $5, "VCF"}' \
  > ${out_prefix}.avinput

${annovar} ${out_prefix}.avinput ${annovar_db} \
  -build hg38 \
  -out ${out_prefix}.annovar \
  -geneanno \
  -dbtype refGene

awk 'BEGIN{FS=OFS="\t"} NR>1 && $1 != "" && $2 != "" {print $1 ":" $2}' \
  ${whitelist_table} > ${tmpdir}/exonic_patterns.txt
awk 'BEGIN{FS=OFS="\t"} NR>1 && $2 != "" {print $2}' \
  ${whitelist_table} > ${tmpdir}/transcript_patterns.txt

grep -F -f ${tmpdir}/exonic_patterns.txt \
  ${out_prefix}.annovar.exonic_variant_function \
  > ${out_prefix}.exonic.whitelist.txt || true

grep -F -f ${tmpdir}/transcript_patterns.txt \
  ${out_prefix}.annovar.variant_function \
  > ${out_prefix}.splicing.whitelist.txt || true

awk -F'\t' '$2 ~ /stopgain|nonsynonymous SNV|frameshift insertion|frameshift deletion|frameshift block substitution/' \
  ${out_prefix}.exonic.whitelist.txt \
  > ${out_prefix}.exonic.damage.txt

awk -F'\t' '$1 ~ /splicing/' \
  ${out_prefix}.splicing.whitelist.txt \
  > ${out_prefix}.splicing.damage.txt

awk 'BEGIN{OFS="\t"} {print "exonic", $0}' ${out_prefix}.exonic.damage.txt \
  > ${out_prefix}.damage.txt
awk 'BEGIN{OFS="\t"} {print "splicing", $0}' ${out_prefix}.splicing.damage.txt \
  >> ${out_prefix}.damage.txt
