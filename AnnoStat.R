#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("Annotation result statistics")

# Add command line arguments
p <- add_argument(p, "fasta", help="fasta file, used to get sequence name and number", type="character")

# Parse the command line arguments
argv <- parse_args(p)
script_dir <- dirname(strsplit(commandArgs(trailingOnly = FALSE)[4],"=")[[1]][2])

# test --------------------------------------------------------------------

# setwd("/home/zhxd/software/emcp/example_data")
# argv <- list()
# argv$fasta <- "GDDH13_1-1_prot.fasta"
# script_dir <- "/home/zhxd/software/emcp"

# function:make OrgDB -----------------------------------------------------
library(tidyverse)
library(clusterProfiler)
library(formattable)
library(seqinr)

install.packages("org.My.eg.db", repos = NULL, lib = "R_library")
library("org.My.eg.db", character.only = TRUE, lib.loc = "R_library")

annoStat <- function() {
  all_gene <- getName.list(read.fasta(file = argv$fasta, 
                         seqtype = 'AA'))
  
  load("gene_annotation.RData")
  
  # number and percentage ---------------------------------------------------
  total_gene = length(all_gene)
  eggnog_anno = length(gene_info$GID)
  go_anno = length(unique(gene2go$GID))
  cog_anno = length(unique(gene2cog$GID))
  pathway_anno = length(unique(gene2pathway$GID))
  
  anno_stat <- tibble(
    database = c("EggNOG", "GO", "COG/KOG", "KEGG Pathway"),
    number = comma(c(eggnog_anno, go_anno, cog_anno, pathway_anno), digits = 0),
    percentage = percent(c(eggnog_anno, go_anno, cog_anno, pathway_anno)/total_gene)
  )
  
  write.table(anno_stat, "anno_stat.txt", quote = F, row.names = F, sep = "\t")
  
  # GO statistics and plot --------------------------------------------------
  
  go_bp <- groupGO(gene     = all_gene,
                   OrgDb    = org.My.eg.db,
                   keyType  = "GID",
                   ont      = "BP",
                   level    = 2,
                   readable = FALSE)
  
  go_bp <- as.data.frame(go_bp)
  go_bp$GO_Class <- "Biological Process"
  
  go_cc <- groupGO(gene     = all_gene,
                   OrgDb    = org.My.eg.db,
                   keyType  = "GID",
                   ont      = "CC",
                   level    = 2,
                   readable = FALSE)
  
  go_cc <- as.data.frame(go_cc)
  go_cc$GO_Class <- "Cellular Component"
  
  go_mf <- groupGO(gene     = all_gene,
                   OrgDb    = org.My.eg.db,
                   keyType  = "GID",
                   ont      = "MF",
                   level    = 2,
                   readable = FALSE)
  go_mf <- as.data.frame(go_mf)
  go_mf$GO_Class <- "Molecular Function"
  
  go_all <- rbind(go_bp, go_cc, go_mf)
  write.table(go_all, "go.txt", sep = "\t", quote = F)
  p <- ggplot(go_all) + 
    geom_bar(aes(x = Description, 
                       y = Count,
                       fill = GO_Class),
                   stat = "identity") + facet_wrap(~GO_Class, scales = "free_x") + 
    labs(title = "GO function classification", y = "Number of genes") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
          axis.title.x = element_blank(),
          legend.position = "none")
  
  ggsave("go.pdf", p, width = 20, height = 7)
  
  
  # Pathway statistics and plot ---------------------------------------------
  pathway_stat <- dplyr::select(gene2pathway, GID, Pathway_Class, Pathway_Subclass) %>% 
    distinct() %>% 
    group_by(Pathway_Class, Pathway_Subclass) %>%
    summarise(Count = n(), Percentage = percent(n()/pathway_anno))
  
  pathway_stat$Pathway_Subclass <- ordered(pathway_stat$Pathway_Subclass, levels = pathway_stat$Pathway_Subclass) 
  
  ggplot(pathway_stat, aes(x = Pathway_Subclass, y = Percentage)) +
    geom_bar(aes(fill = Pathway_Class), stat = 'identity') +
    geom_text(aes(label = Count), nudge_y = 0.005) +
    scale_y_continuous(labels=percent) + 
    labs(y = "Percent of genes(%)", x ="", fill = "Class") +
    coord_flip() +
    theme_classic()

  ggsave("pathway.pdf", p, width = 20, height = 7)
  write.table(gene2pathway, file = "pathway.txt", sep = "\t", quote = F)
  write.table(pathway_stat, file = "pathway_stat.txt", sep = "\t", quote = F, row.names = F)
  
  
  # COG statistics and plot -------------------------------------------------
  gene2cog$COG_Name = paste("(", gene2cog$COG, ")", gene2cog$COG_Name, sep = " ")
  
  write.table(gene2cog, file = "cog.txt", sep = "\t", quote = F, row.names = F)
  
  p <- ggplot(data = gene2cog) + 
    geom_bar(aes(x = COG, 
                 fill = COG_Name)) +
    labs(title = "COG/KOG Function Classification ", 
         x = "",
         y = "Number of genes") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(),
          legend.key.size=unit(1,"line"),
          legend.text = element_text(size = 7.5)) +
    guides(fill=guide_legend(ncol=1))
  ggsave("cog.pdf", p, width = 16, height = 7)

}

annoStat()