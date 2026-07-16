# Dipin NC Genetic Analysis Workflow

## Overview
This repository contains standardized code templates for the beta-thalassemia genetic analysis workflow. The study builds on the 1000 beta-Thalassemia Genomes Project and matched control resources to investigate phenotype-associated loci and clonal hematopoiesis mutations in Chinese beta-thalassemia genomes.

The code is organized into two main parts. The first part supports the Figure 4 case-control GWAS workflow, including per-chromosome VCF filtering, cohort merge, PCA, and genome-wide association analysis. The second part supports the Figure 6 clonal hematopoiesis variant workflow, including Mutect2 calling, PASS extraction, depth/allele-fraction filtering, ANNOVAR annotation, and candidate damaging variant selection.


## 1. Variant QC, cohort merge, PCA, and GWAS

The VCF workflow can be run by chromosome. The examples below use chromosome 22.

### 1.1. Case VCF QC
code: `code/run_gwas.sh`

Main inputs:
```bash
./data/raw/case/chr22.vcf.gz
./reference/hg38/Homo_sapiens_assembly38.fasta
```

Main outputs:
```bash
./work/gwas/chr22/case.filtered.vcf.gz
./work/gwas/chr22/case.gwas.vcf.gz
./work/gwas/chr22/case.pca.vcf.gz
./data/merge/gwas/case/chr22.vcf.gz
./data/merge/pca/case/chr22.vcf.gz
```

Before running the workflow, chromosome names should use the same style across all cohorts. Removing the `chr` prefix, for example from `chr22` to `22`, is recommended to make the merge step stable.

The workflow applies MAF filtering before cohort merge. GWAS uses MAF 0.002 and PCA uses MAF 0.005.

### 1.2. Merge case, NYUWA control, and 1KGP control
code: `code/run_gwas.sh`

Main inputs:
```bash
./data/merge/gwas/case/chr22.vcf.gz
./data/merge/gwas/1kgp/chr22.vcf.gz
./data/merge/gwas/nyuwa/chr22.vcf.gz
./data/merge/pca/case/chr22.vcf.gz
./data/merge/pca/1kgp/chr22.vcf.gz
./data/merge/pca/nyuwa/chr22.vcf.gz
```

Main outputs:
```bash
./work/gwas/chr22/merged.gwas.norm.vcf.gz
./work/gwas/chr22/merged.gwas.overlap.vcf.gz
./work/gwas/chr22/merged.pca.norm.vcf.gz
./work/gwas/chr22/merged.pca.overlap.vcf.gz
```

### 1.3. Select shared cohort-overlap sites
code: `code/select_overlap.py`

usage:
```bash
./software/python-version/bin/python code/select_overlap.py \
  --input ./work/gwas/chr22/merged.gwas.norm.vcf.gz \
  --output ./work/gwas/chr22/merged.gwas.overlap.vcf
```

### 1.4. PCA and GWAS
code: `code/run_gwas.sh`

Main inputs:
```bash
./work/gwas/chr*/merged.gwas.overlap.vcf.gz
./work/gwas/chr*/merged.pca.overlap.vcf.gz
./data/phenotype/pheno.txt
./data/phenotype/covar.txt
```

Main outputs:
```bash
./work/gwas/allchr/gwas.vcf.gz
./work/gwas/allchr/pca.vcf.gz
./work/gwas/allchr/gwas_input
./work/gwas/allchr/pca_input
./work/gwas/pca/pca
./results/gwas/gwas
```

usage:
```bash
bash code/run_gwas.sh 22
```

The full argument order is:
```bash
bash code/run_gwas.sh chr threads memory
```

For example:
```bash
bash code/run_gwas.sh 22 15 100000
```

## 2. LocusZoom
code: `code/run_locuszoom.sh`

usage:
```bash
bash code/run_locuszoom.sh ./data/locus/regions.tsv
```

The region list should contain three columns:
```text
chr start end
```

## 3. KEGG and GO enrichment
code:
```bash
code/run_enrichment.sh
code/enrichment.R
```

usage:
```bash
bash code/run_enrichment.sh ./data/enrichment/top_genes.tsv ./results/enrichment
```

## 4. Mutect2 somatic variant calling
code: `code/run_mutect2.sh`

usage:
```bash
bash code/run_mutect2.sh tumor_id ./data/bam/tumor.bam ./results/mutect2/tumor \
  ./data/mutect2/chip_whitelist.tsv
```

The clonal hematopoiesis workflow uses Mutect2 tumor-only calling with a cohort-matched germline resource, followed by PASS extraction and candidate-level depth and allele-fraction filtering. Candidate variants are retained when they pass DP, AD, AF, and biallelic-site filters.

After ANNOVAR gene annotation, coding and splicing candidates are filtered separately because ANNOVAR writes them to different files. A single whitelist table is used to derive both matching formats: `gene:transcript` patterns for `exonic_variant_function`, and transcript ID patterns for `variant_function`.

Whitelist format:
```text
gene_symbol	transcript_id	category
DNMT3A	NM_022552	CH_driver
TET2	NM_001127208	CH_driver
ASXL1	NM_015338	CH_driver
PPM1D	NM_003620	CH_driver
```

For exonic records, the script matches entries such as `ASXL1:NM_015338`. For splicing or other non-exonic records, the script matches transcript IDs such as `NM_015338`.

The final candidate damaging set combines:
```text
stopgain
nonsynonymous SNV
frameshift insertion
frameshift deletion
frameshift block substitution
splicing
```

Main outputs:
```bash
./results/mutect2/tumor.candidate.vcf.gz
./results/mutect2/tumor.annovar.exonic_variant_function
./results/mutect2/tumor.annovar.variant_function
./results/mutect2/tumor.exonic.whitelist.txt
./results/mutect2/tumor.splicing.whitelist.txt
./results/mutect2/tumor.exonic.damage.txt
./results/mutect2/tumor.splicing.damage.txt
./results/mutect2/tumor.damage.txt
```
