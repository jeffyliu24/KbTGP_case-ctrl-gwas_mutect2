# Software and Issue Checklist

## Software paths to check

| Tool | Template path | Original line(s) | Note |
|---|---|---:|---|
| bcftools | `bcftools=./software/bcftools-1.10.2/bcftools` | 4, 6, 9, 11, 23, 25, 27, 59, 209, 212-214, 217 | Checked from `/home2/Tools/bin/bcftools --version`. Source mixes plain `bcftools` and `/home2/Tools/bin/bcftools`. |
| bgzip | `bgzip=./software/htslib-1.12/bgzip` | 6, 13, 16, 195 | Checked from `bgzip --version`. |
| tabix | `tabix=./software/htslib-1.12/tabix` | 7, 24, 196, 208 | Checked from `tabix --version`. |
| vcftools | `vcftools=./software/vcftools-0.1.15/vcftools` | 13 | Checked from `vcftools --version`. |
| plink2 | `plink2=./software/plink2-v2.00a5.10LM/plink2` | 15, 62, 76 | Checked from `plink2 --version`: 5 Jan 2024 build. |
| plink | `plink=./software/plink-1.90b6.9/plink` | 68-70, 73 | Checked from `plink --version`: 4 Mar 2019 build. |
| Python | `python=./software/python-3.8.12/bin/python` | 31 | Source uses `/usr/bin/python` 2.6.6, but the cleaned helper script is Python 3 compatible and should use Python 3.8.12. |
| Rscript | `rscript=./software/R-4.3.2/bin/Rscript` | 100 | Checked from `/Parastor300s_G30S/xiaoyu/software/miniconda3/envs/KEGG/bin/Rscript --version`. |
| R packages | R library paths or environment | 102-113 | ggplot2, ggsci, cowplot, ggpubr, CMplot, clusterProfiler, org.Hs.eg.db, enrichplot, stringr, tidyverse, createKEGGdb, KEGG.db. |
| locuszoom | `locuszoom=./software/locuszoom-1.4/bin/locuszoom` | 85 | Checked remotely at `/home2/jeffy/WGS_tools/locuszoom`; `CHANGELOG.md` top version is 1.4, `dist/GITHASH` is `5f46d27`. This installation needs Python 2.7-compatible runtime. |
| GATK | `gatk=./software/gatk-4.1.9.0/gatk` | 181, 187, 200 | Source uses GATK 4.1.9.0 in one variable but plain `gatk` later. |
| ANNOVAR | `annovar=./software/annovar-20200607/annotate_variation.pl` | 218 | Source path suggests version 20200607. |
| awk | system awk | 217, 227, 228 | System tool. |
| grep | system grep | 224, 225 | System tool. |

## Items to confirm before final execution

1. Source line 7 indexes `1kgp.CHS_KHV_CDX.nofam305.chr22.PASS.delchr.vcf.gz`, but line 6 writes `1kgp.chr22.PASS.delchr.vcf.gz`. The template now uses case as the QC example and indexes the file produced by the previous command.

2. Source lines 30-49 write only non-header VCF records. The template writes a complete VCF with headers for downstream `bcftools concat` and PLINK loading.

3. Source lines 42-46 compare REF/ALT by taking fixed characters from the ID string. The template checks cohort tags in the ID and uses the VCF REF/ALT fields directly. Please confirm this logic.

4. Source line 54 starts `for i in {1..22}; do`, but there is no matching `done` before section five. The template removes that unused outer loop and performs one all-chromosome concat step.

5. Source lines 59-76 mix `maf002` and `maf005` input prefixes. The cleaned workflow treats this as a fixed rule before merge: GWAS uses MAF 0.002 and PCA uses MAF 0.005.

6. Source line 76 uses `${phenodir}` and `${newpcadir}`, but these variables are not defined in `daima.txt`. The template defines phenotype and covariate files under `./data/phenotype`.

7. Source lines 187-209 use one tumor ID/input path but later filter a different output prefix. The template uses one `tumor_id`, one BAM, and one output prefix throughout.

8. Source line 212 comment says `DO30`; this is likely `DP30`. The template uses `DP30` in output names.

9. Source line 213 uses `FMT/AD[1]` and `FMT/AF[0]`. Please confirm these fields and indices match the Mutect2 VCF format used in your data.

10. Source line 228 uses `splcing`, likely a typo for `splicing`. The template uses `splicing`.
