#!/usr/bin/env bash
set -euo pipefail

### Tools
rscript=./software/R-4.3.2/bin/Rscript

### Input
gene_table=${1:-./data/enrichment/top_genes.tsv}
outdir=${2:-./results/enrichment}

mkdir -p ${outdir}

${rscript} code/enrichment.R ${gene_table} ${outdir}
