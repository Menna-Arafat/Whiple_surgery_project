#tutorial: https://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html , 
# https://learn.gencore.bio.nyu.edu/rna-seq-analysis/over-representation-analysis/

#set wd
setwd("C:/Users/USER/Documents/pancreas/ppi")
library(plyr)
library(dplyr)
library(igraph)
library(tidyr)
library(readxl)
library(openxlsx)
library(org.Hs.eg.db)
library(scales)
library(ggplot2)
library(tibble)
library("GSEABase")
library("EnrichmentBrowser")
library(clusterProfiler)
library(enrichplot)
library(biomaRt)
system.file(package = "GSEABase")
list.files()


#load data
communities= read.xlsx("communities_louvain_grey.xlsx")[-1,]
#get and prepare gmt 
gmt = geneIds(getGmt("Final Custom copy.gmt" , geneIdType=SymbolIdentifier()))
#convert id
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

id= gmt %>% unlist() %>% unname() %>% .[. != ""]

# Query the database
gene_symbols <- getBM(
  attributes = c("uniprot_gn_id", "hgnc_symbol"),
  filters = "uniprot_gn_id",
  values = id,
  mart = ensembl
)

mapping=  gene_symbols %>%  filter(!is.na(gene_symbols$hgnc_symbol) & gene_symbols$hgnc_symbol != "")
#write.csv(mapping, "idmapping.csv", row.names = F)
mapping= read.csv("idmapping.csv")


#mapping function
map_list= function (mapping_file, list){
  mapping= as.data.frame(mapping_file)
  names(mapping)= c("query", "name")
  list_mapped= list()
  for(i in 1:length(list)){
    list_mapped[[i]]= mapping$name[match(list[[i]] , mapping$query)] %>% .[!is.na(.)] %>% unique()
  }
  names(list_mapped)= names(list)
  return(list_mapped)
}

mapped_gmt= map_list(mapping, gmt)
#convert gmt to long format so that every term has its corresponding gene
gmt_df= ldply(mapped_gmt, rbind)

gmt_long=  gmt_df %>%
  pivot_longer(
    cols= -1,
    names_to = "term",
    values_to = "gene"
  ) %>%
  filter(gene != "" & !is.na(gene)) %>% 
  dplyr::select(-"term")


gmt_long$.id= gsub("\\|.*", "", gmt_long$.id)
gmt_long$.id= gsub("^hsa\\d+\\_", "",gmt_long$.id )
gmt_long$.id= gsub("_/_", "_",gmt_long$.id )
gmt_long$.id= gsub("^GOBP_", "",gmt_long$.id )
gmt_long$.id= gsub(".v\\d+.*", "",gmt_long$.id )

head(gmt_long, n=100)
#----------------------------------------
communities= lapply(communities, function(x) as.list(x))
geneClusters= lapply(communities, function(x) x[x != ""])


#enrichment multiple lists
clProf= llply(geneClusters, function(i) {
  x= enricher(i, TERM2GENE =gmt_long)
  return(as.data.frame(x))})

clusters.levels = names(geneClusters)
clProf.df <- ldply(clProf, rbind)
clProf.df <- plyr::rename(clProf.df, c(.id = "Cluster"))
clProf.df$Cluster = factor(clProf.df$Cluster, levels = clusters.levels)
write.csv(clProf.df, "enrichment_louvain_communities_grey.csv", row.names = F)


