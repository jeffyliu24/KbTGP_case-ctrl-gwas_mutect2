#!/usr/bin/env bash
set -euo pipefail

### Tools
bcftools=./software/bcftools-1.10.2/bcftools
bgzip=./software/htslib-1.12/bgzip
tabix=./software/htslib-1.12/tabix
vcftools=./software/vcftools-0.1.15/vcftools
plink=./software/plink-1.90b6.9/plink
plink2=./software/plink2-v2.00a5.10LM/plink2
python=./software/python-3.8.12/bin/python

### Reference
reference=./reference/hg38/Homo_sapiens_assembly38.fasta

### Input
chr=${1:-22}
threads=${2:-15}
memory=${3:-100000}

raw_case_vcf=./data/raw/case/chr${chr}.vcf.gz
case_gwas_vcf=./data/merge/gwas/case/chr${chr}.vcf.gz
kgp_gwas_vcf=./data/merge/gwas/1kgp/chr${chr}.vcf.gz
ctrl_gwas_vcf=./data/merge/gwas/nyuwa/chr${chr}.vcf.gz
case_pca_vcf=./data/merge/pca/case/chr${chr}.vcf.gz
kgp_pca_vcf=./data/merge/pca/1kgp/chr${chr}.vcf.gz
ctrl_pca_vcf=./data/merge/pca/nyuwa/chr${chr}.vcf.gz
pheno=./data/phenotype/pheno.txt
covar=./data/phenotype/covar.txt

### Output
chr_work=./work/gwas/chr${chr}
all_work=./work/gwas/allchr
pca_work=./work/gwas/pca
gwas_out=./results/gwas

mkdir -p ${chr_work} ${all_work} ${pca_work} ${gwas_out}

### Step 1: case VCF QC
${bcftools} view -f PASS ${raw_case_vcf} -Oz -o ${chr_work}/case.pass.vcf.gz --threads ${threads}
${tabix} -f ${chr_work}/case.pass.vcf.gz

mkdir -p ${chr_work}/tmp
${bcftools} norm -m -both --threads ${threads} -f ${reference} ${chr_work}/case.pass.vcf.gz \
  | ${bcftools} sort -T ${chr_work}/tmp \
  | ${bcftools} annotate --threads ${threads} --set-id "%CHROM\:%POS\[case_hg38\]%REF\,%FIRST_ALT" -e 'ALT="*"' \
    -Oz -o ${chr_work}/case.norm.vcf.gz
${tabix} -f ${chr_work}/case.norm.vcf.gz

${bcftools} filter -i 'QUAL > 30 && INFO/MQ > 20 && INFO/InbreedingCoeff > -0.3' -g 10 \
  ${chr_work}/case.norm.vcf.gz -Oz -o ${chr_work}/case.siteqc.vcf.gz
${tabix} -f ${chr_work}/case.siteqc.vcf.gz

${vcftools} --gzvcf ${chr_work}/case.siteqc.vcf.gz \
  --min-meanDP 15 --minDP 10 --minGQ 20 --recode --stdout \
  | ${bgzip} -c --threads ${threads} > ${chr_work}/case.gtqc.vcf.gz
${tabix} -f ${chr_work}/case.gtqc.vcf.gz

${plink2} --vcf ${chr_work}/case.gtqc.vcf.gz \
  --snps-only --max-alleles 2 --min-alleles 2 \
  --mind 0.1 --geno 0.1 --hwe 0.0001 \
  --recode vcf --out ${chr_work}/case.filtered \
  --threads ${threads} --memory ${memory}
${bgzip} -f ${chr_work}/case.filtered.vcf
${tabix} -f ${chr_work}/case.filtered.vcf.gz

${plink2} --vcf ${chr_work}/case.filtered.vcf.gz \
  --maf 0.002 \
  --recode vcf --out ${chr_work}/case.gwas \
  --threads ${threads} --memory ${memory}
${bgzip} -f ${chr_work}/case.gwas.vcf
${tabix} -f ${chr_work}/case.gwas.vcf.gz
mkdir -p $(dirname ${case_gwas_vcf})
cp ${chr_work}/case.gwas.vcf.gz ${case_gwas_vcf}
cp ${chr_work}/case.gwas.vcf.gz.tbi ${case_gwas_vcf}.tbi

${plink2} --vcf ${chr_work}/case.filtered.vcf.gz \
  --maf 0.005 \
  --recode vcf --out ${chr_work}/case.pca \
  --threads ${threads} --memory ${memory}
${bgzip} -f ${chr_work}/case.pca.vcf
${tabix} -f ${chr_work}/case.pca.vcf.gz
mkdir -p $(dirname ${case_pca_vcf})
cp ${chr_work}/case.pca.vcf.gz ${case_pca_vcf}
cp ${chr_work}/case.pca.vcf.gz.tbi ${case_pca_vcf}.tbi

### Step 2: merge GWAS cohorts
${bcftools} merge ${case_gwas_vcf} ${kgp_gwas_vcf} ${ctrl_gwas_vcf} \
  --missing-to-ref -Oz -o ${chr_work}/merged.gwas.vcf.gz --threads ${threads}
${tabix} -f ${chr_work}/merged.gwas.vcf.gz

${bcftools} norm --threads ${threads} -m -both ${chr_work}/merged.gwas.vcf.gz \
  -Oz -o ${chr_work}/merged.gwas.norm.vcf.gz
${tabix} -f ${chr_work}/merged.gwas.norm.vcf.gz

${python} code/select_overlap.py \
  --input ${chr_work}/merged.gwas.norm.vcf.gz \
  --output ${chr_work}/merged.gwas.overlap.vcf
${bgzip} -f ${chr_work}/merged.gwas.overlap.vcf
${tabix} -f ${chr_work}/merged.gwas.overlap.vcf.gz

### Step 3: merge PCA cohorts
${bcftools} merge ${case_pca_vcf} ${kgp_pca_vcf} ${ctrl_pca_vcf} \
  --missing-to-ref -Oz -o ${chr_work}/merged.pca.vcf.gz --threads ${threads}
${tabix} -f ${chr_work}/merged.pca.vcf.gz

${bcftools} norm --threads ${threads} -m -both ${chr_work}/merged.pca.vcf.gz \
  -Oz -o ${chr_work}/merged.pca.norm.vcf.gz
${tabix} -f ${chr_work}/merged.pca.norm.vcf.gz

${python} code/select_overlap.py \
  --input ${chr_work}/merged.pca.norm.vcf.gz \
  --output ${chr_work}/merged.pca.overlap.vcf
${bgzip} -f ${chr_work}/merged.pca.overlap.vcf
${tabix} -f ${chr_work}/merged.pca.overlap.vcf.gz

### Step 4: PCA and GWAS
${bcftools} concat ./work/gwas/chr*/merged.gwas.overlap.vcf.gz \
  -Oz -o ${all_work}/gwas.vcf.gz --threads ${threads}
${tabix} -f ${all_work}/gwas.vcf.gz

${bcftools} concat ./work/gwas/chr*/merged.pca.overlap.vcf.gz \
  -Oz -o ${all_work}/pca.vcf.gz --threads ${threads}
${tabix} -f ${all_work}/pca.vcf.gz

${plink2} --vcf ${all_work}/gwas.vcf.gz \
  --make-bed --out ${all_work}/gwas_input \
  --const-fid 0 --threads ${threads} --memory ${memory}

${plink2} --vcf ${all_work}/pca.vcf.gz \
  --make-bed --out ${all_work}/pca_input \
  --const-fid 0 --threads ${threads} --memory ${memory}

${plink} --bfile ${all_work}/pca_input \
  --bp-space 1000 --make-bed --keep-allele-order \
  --out ${pca_work}/ld_space --threads ${threads} --memory ${memory}

${plink} --bfile ${pca_work}/ld_space \
  --indep-pairwise 50 5 0.1 \
  --out ${pca_work}/ld_prune --threads ${threads} --memory ${memory}

${plink} --bfile ${pca_work}/ld_space \
  --extract ${pca_work}/ld_prune.prune.in \
  --keep-allele-order --make-bed \
  --out ${pca_work}/pca_input --threads ${threads} --memory ${memory}

${plink} --bfile ${pca_work}/pca_input \
  --pca 20 tabs var-wts --make-rel --keep-allele-order \
  --out ${pca_work}/pca --threads ${threads} --memory ${memory}

${plink2} --bfile ${all_work}/gwas_input \
  --logistic hide-covar --adjust \
  --pheno ${pheno} \
  --covar ${covar} \
  --covar-col-nums 3-15 \
  --no-sex --allow-no-sex \
  --out ${gwas_out}/gwas \
  --threads ${threads} --memory ${memory}
