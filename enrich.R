#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("do GO enrichment")

# Add command line arguments
p <- add_argument(p, "--de_result", help="input de_result file, from run_DE_analysis.pl", type="character")
p <- add_argument(p, "--de_log2FoldChange", help="log2FoldChange cutoff", type="numeric", default = 1)
p <- add_argument(p, "--de_padj", help="adjust pvalue cutoff", type="numeric", default = 0.05)

# Parse the command line arguments
argv <- parse_args(p)
out_prefix <- argv$de_result

script_dir <- dirname(strsplit(commandArgs(trailingOnly = FALSE)[4],"=")[[1]][2])
# load library ------------------------------------------------------------
library(tidyverse)
library(clusterProfiler)
library(DOSE)
library(pathview)

#dir.create("R_library", recursive = TRUE)
#install.packages("org.My.eg.db", repos = NULL, lib = "R_library")

library(org.My.eg.db, lib.loc = "R_library")
select(org.My.eg.db, keys = head(keys(org.My.eg.db), n = 2), columns = c('GID', 'GO', 'COG', 'Ko', 'Pathway', 'GENENAME'))

# load gene list or de_result ---------------------------------------------
de_result <- read.delim(argv$de_result)
geneList <- de_result$log2FoldChange
names(geneList) <- rownames(de_result)


de_result_filter <- mutate(de_result, GID = rownames(de_result)) %>%
  filter(abs(log2FoldChange) > argv$de_log2FoldChange, pvalue < argv$de_padj)

deg <- as.character(de_result_filter$GID)


# do GO enrich -------------------------------------------------------------
ego <- enrichGO(gene          = deg,
                keyType       = "GID",
                OrgDb         = org.My.eg.db,
                ont           = "CC",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                pAdjustMethod = "BH",
                readable      = FALSE)

ego_results<-as.data.frame(ego)
write.table(ego_results, file = paste(out_prefix, "ego_results.txt", sep = "."),
            quote = F, sep = "\t")

pdf(file = paste(out_prefix, "ego_barplot.pdf", sep = "."))
barplot(ego, showCategory=20, x = "GeneRatio")
dev.off()

pdf(file = paste(out_prefix, "ego_dotplot.pdf", sep = "."))
dotplot(ego)
dev.off()

pdf(file = paste(out_prefix, "ego_emapplot.pdf", sep = "."))
emapplot(ego)
dev.off()

pdf(file = paste(out_prefix, "ego_goplot.pdf", sep = "."))
goplot(ego)
dev.off()


# Do Pathway enrich ------------------------------------------------------
pathway2gene <- select(org.My.eg.db, keys = keys(org.My.eg.db), columns = c("Pathway")) %>%
  na.omit() %>%
  dplyr::select(Pathway, GID)

load(paste(script_dir, "kegg_info.RData", sep = "/"))

ekp <- enricher(deg, 
                TERM2GENE = pathway2gene, 
                TERM2NAME = pathway2name, 
                pvalueCutoff = 1, 
                qvalueCutoff = 1,
                pAdjustMethod = "BH",
                minGSSize = 1)

ekp_results <- as.data.frame(ekp)
write.table(ekp_results, file = paste(out_prefix, "ekp_results.txt", sep = "."),
            quote = F, sep = "\t")

pdf(file = paste(out_prefix, "ekp_barplot.pdf", sep = "."))
barplot(ekp, showCategory=20, x = "GeneRatio")
dev.off()

pdf(file = paste(out_prefix, "ekp_dotplot.pdf", sep = "."))
dotplot(ekp)
dev.off()

pdf(file = paste(out_prefix, "ekp_emapplot.pdf", sep = "."))
emapplot(ekp)
dev.off()


# pathway view ------------------------------------------------------------

# id.map <- select(org.My.eg.db, keys = names(geneList), columns = "Ko")
# gene.ko <- mol.sum(mol.data = geneList, id.map = id.map)
# 
# sig.pathway <- as.character(filter(ekp_results, p.adjust < 0.05)$ID)
# 
# work_dir <- getwd()
# pathview_dir <- paste(out_prefix, 'pathwiew', sep = "_")
# dir.create(pathview_dir, recursive=T)
# setwd(pathview_dir)
# pathview(gene.data  = gene.ko,
#                      pathway.id = sig.pathway,
#                      species    = "ko")
# setwd(work_dir)
