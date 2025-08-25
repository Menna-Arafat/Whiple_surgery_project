#title:: Automated WGCNA pipeline

#Author: Menna Arafat

#Date: 3/8/2024

#Description: this script automates the process of conducting Weighted Gene Co-expression Network Analysis (WGCNA) to
# to identify modules of highly correlated genes and relate them to external traits

#project name: tTime series proteomic profiling for patients underwent whipple surgery 

#Input:
#data : expression data matrix where samples in columns and features (genes/ preoteins) in rows
#metadata: matrix that store phenotype in a binary (numeric) format; 0 for control and 1 for disease.
#mapping_file: is a dataframe that contain id in a column and corresponding gene symbol (the one used here is from uniprot)

#output:
#"modules_gene_id.csv" : file that contain modules of genes
#"modules_gene_name.csv": file that contain modules of genes mapped to preferred symbol
#"gene.trait.corr.csv": file that store the correlattion between every gene and the phenotype of interest (the file is used for prioritizing hub genes)
#"module_membership.csv" : file that store the intra-modular connectivity (module membership) (sth similar to node degree,and more precisely it is defined as correlation between gene expression profile and the module eigengene; the first principal component or the mean expression profile of a module )
#"module_hubs.csv": hub genes based on high intra-modular connectivity. (you can specify the threshold using the argument: """"hub_threshold"""")
#"phenotype_cor_hubs.csv": hub genes based on the correlation with the phenotype of interest. (you can specify the threshold using the argument: """"""GS_threshold"""""" )

#-------------------------------------------------------------
#pre-requiste libraries
# install.packages("WGCNA")
# install.packages("flashClust")
# install.packages("ggplot2")
# install.packages("data.table")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("biomaRt")


setwd("C:/Users/USER/Documents/pancreas")
list.files()
#create directtories
dir.create("input")
dir.create("output")
dir.create("plots")

#load libraries
library(WGCNA)
library(stats)
library(flashClust)
library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(tidyverse)
library(gridExtra)
library(biomaRt)
library(expm)
library(gplots)
library(ggplot2)
library(readxl)
allowWGCNAThreads()          # allow multi-threading (optional)

#load data
#expression data
list.files(paste0(getwd(), "/input"))
data= read.csv("input/FOR_METABO_FINAL_7days.csv" )[-1,] 
row.names(data)= NULL
data= data %>% column_to_rownames("Sample")
data %>% names()

#for data it is better to be normalized, corrected for batch effect and filtered from low count genes
data= data[rowSums(data == 0) < 4, ] #keep only rows that have number of zeros < 4
#write.csv(data, "input/pancreas_input333.csv")
#metadata
metadata= data.frame(sample= colnames(data) ,
                     condition= c( rep("Before", 3), rep("After.2d",3), rep("After.7d", 3)))
#metadata= read.csv("input/metadata_test.csv" )
metadata = metadata %>%  column_to_rownames("sample")
metadata$condition= factor(metadata$condition, levels = c("Before", "After.2d", "After.7d"))
design= model.matrix(~ 0+condition , metadata)
design
#-------------------------------------------------------------------------------
##Applying WGCNA pipeline
#choose the power that the adjacency matrix raised for that allow for scale free network
# Choose a set of soft thresholding powers
powers = c(1:15)  # in practice this should include powers up to 20.
# choose power based on SFT criterion
sft = pickSoftThreshold(t(data), powerVector = powers, networkType = "signed")
# Plot the results:
png("plots/elbow_plot.png", width = 1700, height = 1000, res=300)
par(mfrow = c(1, 2))
# SFT index as a function of different powers
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     xlab = "Soft Threshold (power)", ylab = "SFT, signed R^2", type = "n", main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     labels = powers, col = "red", cex=.7)
# this line corresponds to using an R^2 cut-off of h
abline(h = 0.8, col = "red")
# Mean connectivity as a function of different powers
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], type = "n", xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red", cex=.7)

dev.off()

#-------------------------------------------------------------------------------
#set your parameters
power= 10
minModuleSize = 30
metadata_binary= design
hub_threshold= 0.75
GS_threshold= 0.75 #(.2) default
networkType = "signed"
mapping_file= read.table("idmapping_uniprot.tsv", header = T)


#------------------------------------------------------------------------------  
# Main Function:
Run_WGCNA= function(data, metadata_binary, power,minModuleSize, hub_threshold ,GS_threshold, networkType,
                    mapping_file=NULL){
  
  
  data= as.data.frame(data)
  data[]= lapply(data, as.numeric)
  datExpr= t(data) #so that samples become in row
  #This function checks data for missing entries, entries with weights below a threshold, and zero-variance genes, 
  goods <- goodSamplesGenes(datExpr, verbose = 3)
  datExpr= datExpr[goods$goodSamples== TRUE, goods$goodGenes == TRUE ]
  
  net = blockwiseModules(datExpr, corType = "pearson", maxBlockSize = 5000, 
                         networkType = networkType , power = power, minModuleSize =minModuleSize,
                         mergeCutHeight = 0.25, 
                         numericLabels = F, saveTOMs = TRUE, 
                         pamRespectsDendro = FALSE, saveTOMFileBase = "TOM")
  #table(net$colors)
  print("module detection is done")
  
  # A data frame with module eigengenes 
  module_eigengenes <- net$MEs
  write.csv(module_eigengenes, "output/module_eigengenes.csv")
  names(module_eigengenes)= gsub("ME", "", names(module_eigengenes) )
  print("module_eigengenes file is saved")
  
  #gene module assignment
  module.gene.assign= net$colors #vector of colors that assign each gene to its corresponding module
  write.csv(module.gene.assign, "output/gene_module_assignment.csv")
  print("gene_module_assignment file is saved")
  
  #Extract proteins for each modules
  module_names= names(module_eigengenes)
  module_lists= list()
  
  for (i in module_names){
    module_lists[[i]]=  module.gene.assign[module.gene.assign== i] %>% names()
  }
  
  df <- ldply(module_lists, rbind) %>% t()
  df[is.na(df)] <- ""
  df= as.data.frame(df)
  names(df) = df[1,] %>% as.vector()
  modules= df[-1,]
  write.csv(modules, paste0("output/", "modules_gene_id.csv"), row.names = FALSE)
  print("modules file is saved!")
  
  
  #calculate gene significance (a measure that show to what degree the gene is correlated with the phenotype of interest based on pearson correlation, and the authors consider GS > .2 is good enough) 
  nSamples <- nrow(datExpr)
  nGenes <- ncol(datExpr)
  gene.signf.corr <- cor(datExpr, metadata_binary, use = 'p')
  gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)
  write.csv(gene.signf.corr, "output/gene.trait.corr.csv")
  gene.signf.corr= gene.signf.corr  %>% as.data.frame()
  #Module membership (intra-modular connectivity)
  #For a node i. its module membership kME(q) specifies how close node i is to module q. 
  #KME (module membership) defined as correlation between gene expression profile and the module eigengene (the first principal component or the mean expression profile of a module)
  #For example, if the module membership of a gene i (MMblue.i) is close to 1 or -1, this indicate that the gene has a positive or a negative relationship with the module eigengene, i.e., this gene is highly connected to the blue module genes.
  
  datKME = signedKME(datExpr, module_eigengenes)
  write.csv(datKME, "output/module_membership.csv")
  
  
  #extract hub genes
  #Extract hubs for each modules
  module_names= names(module_eigengenes)
  hub_lists <- vector("list", length(module_names))
  names(hub_lists)= names(module_eigengenes)
  threshold= hub_threshold
  
  for(j in 1:length(module_names)){
    hub_mod= row.names(datKME)[datKME[,j] >= threshold]
    hub_lists[[j]]= intersect(unlist(modules[j]),  hub_mod)
  }
  
  df= ldply(hub_lists, rbind) %>% t() 
  df[is.na(df)] <- ""
  df= as.data.frame(df)
  names(df) = paste0(unlist(df[1,]), "_hub") 
  hubs= df[-1,]
  write.csv(hubs, paste0("output/", "module_hubs.csv"), row.names = FALSE)
  print("hubs file is saved!")
  
  #Extract hubs correlated with phenotype
  hub_pheno <- vector("list", length= ncol(gene.signf.corr))
  names(hub_pheno)= names(gene.signf.corr)
  GS_threshold= GS_threshold
  
  for(j in 1:ncol(gene.signf.corr)){
    
    all_hubs= hub_lists %>% unlist() %>% .[. !="" & !is.na(.)]
    hub_pheno[[j]]= row.names(gene.signf.corr)[gene.signf.corr[,j] > GS_threshold] %>%  
      .[. %in% all_hubs]
  }
  df= ldply(hub_pheno, rbind) %>% t() 
  df[is.na(df)] <- ""
  df= as.data.frame(df)
  names(df) = paste0(unlist(df[1,]), "_hub") 
  hubs= df[-1,]
  write.csv(hubs, paste0("output/", "phenotype_cor_hubs.csv"), row.names = FALSE)
  print("hubs_pheno file is saved!")
  
  #MAPPING 
  if(!is.null(mapping_file)){
    #mapping of hubs for modules and for phenotype
    mapping= mapping_file
    names(mapping)= c("query", "name")
    
    modules= module_lists
    modules_mapped= vector("list", length= length(module_lists))
    names(modules_mapped) = names(module_lists)
    
    for(i in seq_along(modules)){
      modules_mapped[[i]]= mapping$name[match(modules[[i]], mapping$query) ]
    }
    df= ldply( modules_mapped, rbind) %>% t() %>% as.data.frame()
    df[is.na(df)]= ""
    names(df)= as.character(df[1,])
    write.csv(df[-1,], "output/modules_gene_name.csv", row.names = F)
    
    
    combined_hubs= list(mod= hub_lists, pheno= hub_pheno)
    
    mapped_lists = list(mod= setNames(vector("list", length = length(hub_lists)), names(hub_lists)),
                        pheno= setNames(vector("list", length = length(hub_pheno)), names(hub_pheno)))
    
    unmapped_lists = list(mod= setNames(vector("list", length = length(hub_lists)), names(hub_lists)),
                          pheno= setNames(vector("list", length = length(hub_pheno)), names(hub_pheno)))
    
    for(i in seq_along(combined_hubs)){
      for (j in seq_along(combined_hubs[[i]])) {
        mapped_lists[[i]][[j]] <- mapping$name[match(combined_hubs[[i]][[j]], mapping$query)] 
        unmapped_lists[[i]][[j]] = combined_hubs[[i]][[j]][is.na(match(combined_hubs[[i]][[j]], mapping$query)) ]
      }
      df_mapped= ldply(mapped_lists[[i]], rbind) %>% t() %>%  as.data.frame()
      df_mapped[is.na( df_mapped)]= ""
      names(df_mapped)=  as.vector(df_mapped[1,])
      mapped_hubs=  df_mapped[-1,]
      write.csv(mapped_hubs, paste0("output/", names(combined_hubs)[i],"_hubs_names.csv"), row.names = FALSE)
      print("hubs_name file is saved!")
      
      df_unmapped= ldply(unmapped_lists[[i]], rbind) %>% t() %>%  as.data.frame()
      df_unmapped[is.na( df_unmapped)] <- ""
      names(df_unmapped)=  as.vector(df_unmapped[1,])
      unmapped_hubs=  df_unmapped[-1,]
      write.csv(unmapped_hubs, paste0("output/", names(combined_hubs)[i],"_hubs_unmapped.csv"), row.names = FALSE)
      print("hubs_unmapped file is saved!")
      
    }
  }
  else 
  { print("Supply a mapping id file from uniprot database! and Set mapping_source= \"uniprot\"")}
  
  #export networks for cytoscape
  module_names= names(module_eigengenes)
  for( i in   module_names) {
    datexpr_mod = datExpr[,module.gene.assign == i]
    TOM_mod = TOMsimilarityFromExpr(datexpr_mod, power = power, networkType = "signed", TOMType="signed Nowick");
    genes = colnames(datexpr_mod)
    dimnames(TOM_mod) = list(genes, genes)
    # nTop = 30;
    # IMConn = softConnectivity(datExpr[, probes]);
    # top = (rank(-IMConn) <= nTop)
    cyt = exportNetworkToCytoscape(TOM_mod,
                                   edgeFile = paste("output/CytoscapeInput-edges-.02tom", paste(i, collapse="-"), ".txt", sep=""),
                                   #nodeFile = paste("output/CytoscapeInput-nodes-", paste(i, collapse="-"), ".txt", sep=""),
                                   weighted = TRUE,
                                   threshold = 0.02, #threshold for including edges in the output, Default is 0.02
                                   nodeNames = genes,
                                   altNodeNames = genes);
  }
  print("networks exported!")
  print("ALL IS DONE!")
  return(net)
  
}

#-------------------------------------------------------------------------------
#Run the function
res= Run_WGCNA(data, metadata_binary=metadata_binary, power=power, minModuleSize= minModuleSize,networkType = networkType,
               hub_threshold= hub_threshold, GS_threshold=GS_threshold, mapping_file=mapping_file )

#######################################################################################################################################