#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("make OrgDb from emapper")

# Add command line arguments
p <- add_argument(p, "annotation", help="emapper annotation result", type="character")

# Parse the command line arguments
argv <- parse_args(p)
script_dir <- dirname(strsplit(commandArgs(trailingOnly = FALSE)[4],"=")[[1]][2])

# function:make OrgDB -----------------------------------------------------
library(tidyverse)
library(stringr)
library(AnnotationForge)
library(formattable)
library(clusterProfiler)
#' make OrgDB Package From eggnog-mapper result
#'
#' @param f_emapper_anno eggnog-mapper annotation result
#' @param author Who is the creator of this package? like "xxx <xxx@xxx.xx>"
#' @param tax_id The Taxonomy ID that represents your organism. (NCBI has a nice online browser for finding the one you need)
#' @param genus Single string indicating the genus
#' @param species Single string indicating the species
#'
#' @return OrgDb name
#' @export
#'
#' @examples
makeOrgPackageFromEmapper <- function(f_emapper_anno, 
                                      author, 
                                      tax_id = "0", 
                                      genus = "default", 
                                      species = "default") {
  # test
  # setwd("/home/zhxd/software/emcp/example_data")
  # f_emapper_anno <- "my.emapper.annotations"
  # author <- "zxd"
  # tax_id = "0"
  # genus = "default"
  # species = "default"
  # script_dir <- "/home/zhxd/software/emcp"
  
  # read emapper result
  emapper <- read_delim(f_emapper_anno,
                        "\t", 
                        escape_double = FALSE, 
                        trim_ws = TRUE,
                        skip = 3)
  names(emapper)[1] <- "query_name"
  # extract gene name from emapper ------------------------------------------
  gene_info <- emapper %>%
    dplyr::select(GID = query_name, GENENAME = `eggNOG annot`) %>%
    na.omit()
  
  # extract go annotation from emapper --------------------------------------
  gos <- emapper %>%
    dplyr::select(query_name, GO_terms) %>%
    na.omit()
  
  gene2go = data.frame(GID = character(),
                       GO = character(),
                       EVIDENCE = character())
  
  df_temp <- list()
  for (row in 1:nrow(gos)) {
    the_gid <- gos[row, "query_name"][[1]]
    the_gos <- str_split(gos[row,"GO_terms"], ",", simplify = FALSE)[[1]]
    
    df_temp[[row]] <- tibble(GID = rep(the_gid, length(the_gos)),
                                 GO = the_gos,
                                 EVIDENCE = rep("IEA", length(the_gos)))
  }
  
  gene2go <- bind_rows(df_temp)

  # extract kegg pathway annotation from emapper ----------------------------
  kos <- emapper %>%
    dplyr::select(query_name, KEGG_KOs) %>%
    na.omit()
  
  gene2ko = data.frame(GID = character(),
                       Ko = character())
  
  df_temp <- list()
  for (row in 1:nrow(kos)) {
    the_gid <- kos[row, "query_name"][[1]]
    the_kos <- str_split(kos[row,"KEGG_KOs"], ",", simplify = FALSE)[[1]]
    
    df_temp[[row]] <- data_frame(GID = rep(the_gid, length(the_kos)),
                                 Ko = the_kos)
  }
  
  gene2ko <- bind_rows(df_temp)
  
  load(file = paste(script_dir, "kegg_info.RData", sep = "/"))
  gene2pathway <- gene2ko %>% left_join(ko2pathway, by = "Ko") %>%
    left_join(pathway2name, by = "Pathway") %>%
    dplyr::select(GID, Ko, Pathway, Pathway_Name, Pathway_Class, Pathway_Subclass) %>%
    distinct() %>%
    na.omit()
  
  # extract COG annotation from emapper -------------------------------------
  cog_info <- read_delim(paste(script_dir, "cog_funclass.tab", sep = "/"), 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
  
  cogs <- emapper %>%
    dplyr::select(query_name, COG = `COG cat`) %>%
    na.omit()
  
  gene2cog = data.frame(GID = character(),
                        COG = character())
  
  df_temp <- list()
  for (row in 1:nrow(cogs)) {
    the_gid <- cogs[row, "query_name"][[1]]
    the_cogs <- str_trim(str_split(cogs[row,"COG"], ",", simplify = FALSE)[[1]])
    
    df_temp[[row]] <- data_frame(GID = rep(the_gid, length(the_cogs)),
                                 COG = the_cogs)
  }
  gene2cog <- bind_rows(df_temp)
  
  gene2cog <- gene2cog %>% left_join(cog_info, by = "COG")
  
  save(gene_info, gene2go, gene2pathway, gene2cog, file = "gene_annotation.RData")
  # make OrgDb --------------------------------------------------------------
  makeOrgPackage(gene_info=gene_info,
                 go=gene2go,
                 #ko=gene2ko,
                 pathway=gene2pathway,
                 cog=gene2cog,
                 maintainer=author,
                 author=author,
                 outputDir="./",
                 tax_id=tax_id,
                 genus=genus,
                 species=species,
                 goTable="go",
                 version="1.0")
  
  my_orgdb <- str_c("org.", str_to_upper(str_sub(genus, 1, 1)) , species, ".eg.db", sep = "")
  return(my_orgdb)
}

# makeOrgPackage ----------------------------------------------------------
my_orgdb <- makeOrgPackageFromEmapper(argv$annotation, 
                                      "test <test@genek.tv>", 
                                      tax_id = "0000", 
                                      genus = "M", 
                                      species = "y")
dir.create("R_library")
install.packages(my_orgdb, repos = NULL, lib = "R_library")

library(my_orgdb, character.only = TRUE, lib.loc = "R_library")

select(org.My.eg.db, keys = head(keys(org.My.eg.db), n = 2), columns = c('GID', 'GO', 'COG', 'Ko', 'Pathway', 'GENENAME'))
