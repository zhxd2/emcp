library(jsonlite)
library(tidyverse)

update_kegg <- function(json = "ko00001.json") {
  pathway2name <- tibble(Pathway = character(), Pathway_Name = character(), Pathway_Class = character())
  ko2pathway <- tibble(Ko = character(), Pathway = character())
  
  kegg <- fromJSON(json)
  
  for (a in seq_along(kegg[["children"]][["children"]])) {
    A <- kegg[["children"]][["name"]][[a]]
    
    for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
      B <- kegg[["children"]][["children"]][[a]][["name"]][[b]]
      
      for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
        pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
        
        if(str_detect(pathway_info, "PATH:ko[0-9]{5}")){
          pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
          pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "")
          pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Pathway_Name = pathway_name, Pathway_Class = str_sub(A,7), Pathway_Subclass = str_sub(B, 7)))
          
          kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
          
          kos <- str_match(kos_info, "K[0-9]*")[,1]
          
          ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))
        }
        next
      }
    }
  }
  pathway2name <- pathway2name %>% na.omit()
  ko2pathway <- ko2pathway %>% na.omit()
  save(pathway2name, ko2pathway, file = "kegg_info.RData")
}

update_kegg(json = "ko00001.json")