#tutorial link: https://www.biostars.org/p/461026/   ,  https://support.bioconductor.org/p/9145644/
#set directory
getwd()
setwd("C:/Users/USER/Documents/pancreas")

#load libraries
BiocManager::install("sva")
library(sva)
library(data.table)
library(dplyr)
library(limma)
#load data
list.files(paste0(getwd(), "/input"))
data= read.csv("input/pancreas_alldata.csv" ) %>% column_to_rownames("ID")
data %>% names()
#metadata
metadata= data.frame(sample= colnames(data) ,
                     condition= c( rep("Before", 10), rep("After.2d",10), rep("Before",3),
                                   rep("After.2d", 3), 
                                   rep("After.7d", 3)),
                     batch= c(rep("batch.1", 20), rep("batch.2", 9)))

#write.csv(metadata,"input/pancreas_metadata.csv")
metadata %>% names()
metadata = metadata %>%  column_to_rownames("sample")


#PCA
#For PCA one commonly uses the most variable genes in the dataset, here we use the top-1000 most variable genes
# Calculate rowwise variance
# var <- apply(data, 1, var)
# 
# # Sort decreasingly and take top 1000
# o <- order(var, decreasing=TRUE)
# top1000 <- head(o, 1000)
# 
# # subset data to include only the most variable genes
# data_top1000 <- data[top1000,]

# Run PCA
pca <- prcomp(t(data))

# Combine PCA coordinates with the metadata from the DGEList
to_plot <- data.frame(pca$x, metadata)

# Calculate how many % of total variance is explained by each principal component
#sdev is the standard deviation sowe square it to get back the variance
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
percentVar
# We focus here on PC1 and PC2
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
labs

p2= ggplot(to_plot, aes(x=PC1, y=PC2, color=batch, shape=condition)) + 
  geom_point(size=3) +
  theme_bw()+
  xlab(labs[1]) + ylab(labs[2])

ggsave( "plots/PCA_after_batch.correction.png",p2, width= 9, height= 8, dpi=600)


#apply batch correction with combat-seq
#data_corrected <- limma::removeBatchEffect(data, batch = factor(metadata$batch))

#apply batch effect correction with batch-seq from sva package
design = model.matrix(~ condition, data= metadata)
design
data_corrected <- ComBat_seq(counts = data, batch=factor(metadata$batch),
                             group=NULL, covar_mod = design)

#--------------------------------------------------------
# Run PCA
pca <- prcomp(t(data_corrected))

# Combine PCA coordinates with the metadata from the DGEList
to_plot <- data.frame(pca$x, metadata)

# Calculate how many % of total variance is explained by each principal component
#sdev is the standard deviation sowe square it to get back the variance
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
percentVar
# We focus here on PC1 and PC2
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
labs

p2= ggplot(to_plot, aes(x=PC1, y=PC2, color=batch, shape=condition)) + 
  geom_point(size=3) +
  theme_bw()+
  xlab(labs[1]) + ylab(labs[2])

ggsave( "plots/PCA_before_batch.correction.png",p1, width= 9, height= 8, dpi=600)

