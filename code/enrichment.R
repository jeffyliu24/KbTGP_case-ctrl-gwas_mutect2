args <- commandArgs(trailingOnly = TRUE)
gene_table <- args[1]
outdir <- args[2]

library(ggplot2)
library(ggsci)
library(cowplot)
library(ggpubr)
library(CMplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(stringr)
library(tidyverse)
library(createKEGGdb)
library(KEGG.db)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

top_genes <- read.table(gene_table, header = TRUE, sep = "\t")

codings <- top_genes %>%
  filter(gene_type %in% c("protein_coding", "IG_V_gene", "TR_C_gene",
                          "TR_D_gene", "TR_J_gene", "TR_V_gene")) %>%
  select(gene_name)

case_gene <- codings$gene_name
case_gene_set <- bitr(case_gene, "SYMBOL",
                      toType = c("ENTREZID"), "org.Hs.eg.db",
                      drop = FALSE)

case_kegg <- enrichKEGG(case_gene_set$ENTREZID,
                        keyType = "kegg",
                        pvalueCutoff = 0.2,
                        qvalueCutoff = 0.2,
                        pAdjustMethod = "BH",
                        minGSSize = 5,
                        maxGSSize = 500,
                        organism = "hsa",
                        use_internal_data = TRUE)
case_kegg_readable <- setReadable(case_kegg, "org.Hs.eg.db", "ENTREZID")

case_go <- enrichGO(case_gene_set$ENTREZID,
                    OrgDb = "org.Hs.eg.db",
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    keyType = "ENTREZID",
                    pvalueCutoff = 0.2,
                    qvalueCutoff = 0.2,
                    readable = TRUE)

go_bp <- enrichGO(case_gene_set$ENTREZID,
                  OrgDb = "org.Hs.eg.db",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  keyType = "ENTREZID",
                  pvalueCutoff = 1,
                  qvalueCutoff = 1,
                  readable = TRUE)

write.csv(as.data.frame(case_kegg_readable),
          file.path(outdir, "kegg.csv"),
          row.names = FALSE)
write.csv(as.data.frame(case_go),
          file.path(outdir, "go.csv"),
          row.names = FALSE)
write.csv(as.data.frame(go_bp),
          file.path(outdir, "go_bp.csv"),
          row.names = FALSE)

pdf(file.path(outdir, "kegg_barplot.pdf"), width = 8, height = 6)
barplot(case_kegg, showCategory = 15, title = "KEGG enrichment")
dev.off()

pdf(file.path(outdir, "go_dotplot.pdf"), width = 8, height = 6)
dotplot(case_go, showCategory = 20, title = "GO enrichment")
dev.off()
