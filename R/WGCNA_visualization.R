#Auther: Menna Arafat
# This code is based on the tutorial "Corrected R code from chapter 12 of the book"
# by Steve Horvath, available at https://pages.stat.wisc.edu/~yandell/statgen/ucla/WGCNA/wgcna.html


#-------------------------------------------------
#set workind directory
getwd()
setwd("C:/Users/USER/Documents/pancreas")

#load packages
library(gtools)
library(pROC)
library(ape)
library(ggdendro)
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
library(gplots)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
allowWGCNAThreads()          # allow multi-threading (optional)

#load data
list.files(getwd())
data <- read.csv("input/pancreas_input333.csv") %>% column_to_rownames("X")
#set row names
#row.names(data)= NULL
#data = data %>% column_to_rownames(var = "X")

#change type to numeric
data[]= lapply(data, as.numeric)
datExpr= t(data) #samples become in row
#This function checks data for missing entries, entries with weights below a threshold, and zero-variance genes, 
goods <- goodSamplesGenes(datExpr, verbose = 3)
datExpr= datExpr[goods$goodSamples== TRUE, goods$goodGenes == TRUE ]
#------------------------------------------------------------------------------
#metadata
metadata= data.frame(sample= colnames(data) ,
                     condition= c( rep("Before", 3), rep("After.2d",3), rep("After.7d", 3)))
metadata = metadata %>%  column_to_rownames("sample")
metadata$condition= factor(metadata$condition, levels = c("Before", "After.2d", "After.7d"))
#binarize the categorical columns
design= model.matrix(~ 0+condition , metadata)
design
#------------------------------------------------------------------------------
#parameters for WGCNA
power= 10
minModuleSize = 30
metadata_binary= design

#Run WGCNA
net = blockwiseModules(datExpr, corType = "pearson", maxBlockSize = 5000, 
                       networkType = "signed", power = power, minModuleSize =minModuleSize,
                       mergeCutHeight = 0.25, 
                       numericLabels = F, saveTOMs = TRUE, 
                       pamRespectsDendro = FALSE, saveTOMFileBase = "TOM")
#------------------------------------------------------------------------------
#visualizations
#Quality check of samples/ check outlier samples
#hierarchical clustering for samples/ samples should be in rows
# Finding distance matrix
distance_mat <- dist(t(data), method = 'euclidean')
distance_mat
# Fitting Hierarchical clustering Model 
# to training dataset
set.seed(240) # Setting seed
Hierar_cl <- hclust(distance_mat, method = "average")
# Plotting dendrogram
dir.create("plots")
png("plots/dendrogram_hclust_samples.png", width = 8000, height = 6000, res= 600) 
plot(Hierar_cl)
# # Cutting tree by height
#abline(h = 1.5e11, col = "green")
# # Cutting tree by no. of clusters
#rect.hclust(Hierar_cl, k = 3, border = "blue")
dev.off()

#-------------------------------------------------------------------------------
#tree and modules and heatmap of associated traits/ phenotypes
#Hierarchical clustering of samples, detect outlier samples,and association of sample with certain trait
#with heatmap of such trait where red indicate high value

#Build adjacency matrix for samples
A = adjacency(data, type = "distance")
# this calculates the whole network connectivity
k = as.numeric(apply(A, 2, sum)) - 1
# standardized connectivity
Z.k = scale(k)
# Designate samples as outlying if their Z.k value is below the threshold
thresholdZ.k = -5  # often -2.5

# the color vector indicates outlyingness (red)
outlierColor = ifelse(Z.k < thresholdZ.k, "red", "black")

# calculate the cluster tree using flahsClust or hclust
sampleTree = flashClust(as.dist(1 - A), method = "average")
# Convert traits to a color representation: where red indicates high
# values
traitColors = data.frame(numbers2colors(as.numeric(metadata$cond_binary), signed = TRUE))
#dimnames(traitColors)[[2]] = "Inflammation_lvl"
datColors = data.frame(outlier_Samples = outlierColor, Condition= traitColors)
colnames(datColors)[2]= "Condition"
# Plot the sample dendrogram and the colors underneath.
png("plots/WGCNA_dendrogram.png", width = 8000, height = 6000, res= 600)
plotDendroAndColors(sampleTree, groupLabels = names(datColors), colors = datColors, cex.rowText = 5,
                    main = "Sample dendrogram and Homogeneity of samples heatmap")

#grid.arrange(a1, a2, nrow = 2)
dev.off()

#-----------------------------------------------------------------
# Plot the dendrogram and the module colors before and after merging underneath
png("plots/dendrogram_merged_modules.png", width = 2200, height = 2500, res= 600)
plotDendroAndColors(net$dendrograms[[1]],  net$colors,
                    paste0("Modules"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    cex.colorLabels = 0.6,
                    guideHang = 0.05)
dev.off()
#for error, Error in .plotOrderedColorSubplot(order = order, colors = colors, rowLabels = rowLabels,  : 
#Length of colors vector not compatible with number of objects in 'order'.
#set good samples and good genes
#-----------------------------------------------------------
#  TOM plot/ heatmap of modules for all proteins
dissTOM= 1 - TOMsimilarityFromExpr(datExpr, power= power) #datExpr samples in rows
dendro= net$dendrograms[[1]]
moduleColorsAutomatic= net$colors

#visualizations
png("plots/TOM_PLOT_module_heatmap_proteins_all.png", width = 800, height = 600)
#myheatcol = colorpanel(250,'gold',"orange",'darkred')
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
# Transform dissTOM with a power to enhance visibility
TOMplot(dissTOM, dendro, moduleColorsAutomatic,col= myheatcol, 
        main = "Module Heatmap Plot, All Proteins")
dev.off()

#--------------------------------------------------------------------
#Module trait correlation
# Next use a single trait/variable or the whole metadata binarized to define the module significance 
#what module associated to what phenotype
#trait= metadata$diabetes
traits= design %>% as.data.frame()
head(traits)
# Define numbers of genes and samples
nSamples <- nrow(datExpr)
nGenes <- ncol(datExpr)
module_eigengenes= read.csv("output/module_eigengenes.csv") %>% column_to_rownames("X")
module.trait.corr <- WGCNA::cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

#module_trait heatmap of WGCNA package
# correlations and their p-values
png("plots/heat2_module_trait_WGCNA.png", width = 4000, height = 4500, res= 600) 
textMatrix = paste(signif(module.trait.corr, 2), "\n(", signif(module.trait.corr.pvals, 1), ")", 
                   sep = "")
dim(textMatrix) = dim(module.trait.corr)
par(mar = c(6, 6, 4, 6))
color= colorpanel(250, "#F7F5F4","#80CDC1")
#get shades of a color
colfunc <-colorRampPalette(c("#BA6756","#F7F5F4" ))
colfunc(10)

#color= greenWhiteRed(50)
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = module.trait.corr, xLabels = gsub("condition", "", names(traits)), 
               yLabels = names(module_eigengenes), 
               ySymbols = names(module_eigengenes), colorLabels = FALSE, colors = color, 
               textMatrix = textMatrix, setStdMargins = T, cex.text = 0.8,
               zlim = c(-1, 1),xColorWidth = 1 * strheight("M"),
                yColorWidth = 1.5 * strwidth("M"),xColorOffset = strheight("M")/6, 
               yColorOffset = strwidth("M")/6, font.lab.x = 2, cex.legendLabel = 2,
               font.lab.y = 2, xLabelsAngle = 75,
               main = paste("Module-Condition Relationship"), plotLegend= TRUE)

dev.off()
#------------------------------------------------------------------
#MDS plot
png("plots/MDS_plot.png",  width = 2800, height = 3300, res= 600) 
dissTOM= 1 - TOMsimilarityFromExpr(datExpr, power= 5)
cmd1=cmdscale(as.dist(dissTOM),2)
par(mfrow=c(1,1))
plot(cmd1,col=moduleColorsAutomatic,main="MDS plot",
     xlab="Scaling Dimension 1",ylab="Scaling Dimension 2")
dev.off()

#-------------------------------------------------------------
#-------------------------------------------------------------
#-------------------------------------------------------------
#  TOM plot/ heatmap of modules for selected proteins based on fold change
#calculate fold change for proteins 
names(data)
FC_A2.B= foldchange(apply(data[,grepl("A2",names(data))],1, mean),apply(data[,grepl("B",names(data))],1, mean)) %>%
          as.data.frame() %>% filter(. >= 1.5 | .<= -1.5) %>% row.names()
FC_A7.A2= foldchange(apply(data[,grepl("A7",names(data))],1, mean),apply(data[,grepl("A2",names(data))],1, mean)) %>% 
          as.data.frame() %>% filter(. >= 1.5 | .<= -1.5) %>% row.names()
FC_A7.B= foldchange(apply(data[,grepl("A7",names(data))],1, mean),apply(data[,grepl("B",names(data))],1, mean)) %>% 
         as.data.frame() %>% filter(. >= 1.5 | .<= -1.5) %>% row.names()


proteins= c(FC_A2.B, FC_A7.A2, FC_A7.B) %>% unique()
# subset data to have only selected proteins
datExpr_subset= datExpr[,colnames(datExpr) %in% proteins]
dissTOM_subset= 1 - TOMsimilarityFromExpr(datExpr_subset, power= power)
dendro_subset = hclust(as.dist(dissTOM_subset), method = "average")
module.gene.assign= net$colors
moduleColors_subset= module.gene.assign[names(module.gene.assign)%in% proteins]

png("plots/module_heatmap_TOM_PLOT_selected.png", width = 2800, height = 3300, res= 600)
#myheatcol = colorpanel(250,'gold',"orange",'darkred')
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
# Transform dissTOM with a power to enhance visibility
TOMplot(dissTOM_subset, 
        dendro_subset, 
        moduleColors_subset,col= myheatcol,
        main = "Module Heatmap Plot")
dev.off()


#-----------------------------------------------------------
#Bar plot for phenotype discrimination by module eigengene
list.files(paste0(getwd(), "/output"))
module_eigengenes= read.csv("output/module_eigengenes.csv") %>% 
                        column_to_rownames("X") %>% lapply(., as.vector)
metadata= read.csv("input/metadata_binarized.csv" )
#Module eigengenes of WGCNA to statistically distinguish between the phenotypes of my data via Kruskal-Wallis test
shapiro.test(module_eigengenes$MEturquoise)

mods= names(module_eigengenes)
pvalues1= numeric(length = length(mods))

#check the differentiation between 2 groups (cancer/control) (non parametric) >> mann whitny
for (i in seq_along(mods)){
  wilcox_res1 <- wilcox.test(module_eigengenes[[i]] ~ metadata$cancer_status , paired= F)
  pvalues1[i] <-  wilcox_res1$p.value
}
fdr1 <- p.adjust(pvalues1, method = "BH") # Benjamini
names(fdr1)= mods
fdr1
res1= data.frame(pvalue= pvalues1,
                    fdr= fdr1)

res1$sig= ifelse(res1$fdr <= .05, "***", "")
res1 = res1%>% rownames_to_column("modules")
res1$modules= gsub("ME", "",res1$modules )
write.csv(res1, "output/module_sig_distinguish_cancer_control.csv")
#check the differentiation between 2 groups alveolar/ embryonal  (non parametric) >> wilcoxin
module_eigengenes= read.csv("output/module_eigengenes.csv") %>% 
                               filter(grepl("Alv|Emb" ,X))  %>%
                               column_to_rownames("X")  %>%
                               lapply(., as.vector)
mods= names(module_eigengenes)
pvalues2= numeric(length = length(mods))
for (i in seq_along(mods)){
  wilcox_res <- wilcox.test(module_eigengenes[[i]] ~ metadata[1:35, "con"], paired= F) #2 numeric vectors of moduleigengene for two groups
  pvalues2[i] <- wilcox_res$p.value
}
fdr2 <- p.adjust(pvalues2, method = "BH") # Benjamini
names(fdr2)= mods
fdr2
res2= data.frame(pvalue= pvalues2,
                    fdr= fdr2)
res2
res2$sig= ifelse(res2$fdr <= .05, "***", "")
res2 = res2%>% rownames_to_column("modules")
res2$modules= gsub("ME", "",res2$modules )
 write.csv(res2, "output/module_sig_distinguish_Alv_Emb.csv")


#visualization
p1=ggplot(res1) +
  geom_bar(aes(x = modules, y = -log10(fdr), fill = modules), stat = "identity") +
  coord_flip() + # Flip coordinates to make the bar plot horizontal
  scale_fill_identity() + # Use actual colors specified in the data frame
  geom_text(aes(x =  modules  , y =-log10(fdr) , label = sig), 
            position = position_dodge(width = 0.9), vjust = -0.5, size=.5 ,color = "black") +
  theme_minimal() + 
  labs(y = "-log10(FDR)", x = "Modules") +
  ggtitle("RMS/Ctrl Discrimination by Module Eigengene") + # Label axes and provide a title
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size=12),
        axis.title=element_text(size=13) ) + 
  ylim(0,12) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", lwd= 1)

p1
-log10(res1$fdr)
#+
# annotate("text", x = .98, y = Inf, label = paste0("-------", "\n", "Sig. Cut-off"),
#   hjust = 1.05, vjust = 0, color = "black", size = 4.5, angle = 0)
p2=ggplot(res2) +
  geom_bar(aes(x = modules, y = -log10(fdr), fill = modules), stat = "identity") +
  coord_flip() + # Flip coordinates to make the bar plot horizontal
  scale_fill_identity() + # Use actual colors specified in the data frame
  geom_text(aes(x =  modules  , y =-log10(fdr) , label = sig), 
            position = position_dodge(width = 0.9), vjust = -0.5, size=.5, color = "black") +
  theme_minimal() + # Use a minimal theme for the plot
  labs(y = "-log10(FDR)", x = "Modules") +
  ggtitle("Alv/Emb RMS Discrimination by Module Eigengene") + # Label axes and provide a title
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size=12),
        axis.title=element_text(size=13) ) + 
  ylim(0, 10) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", lwd= 1)

p2
p= grid.arrange(p1, p2, ncol = 2)

ggsave("plots/modules_discrimination_new.png", p,width= 10, height = 7, dpi=900 )
#dev.off()
#-----------------------------------------------------------------------------
#chord plot of key drivers
list.files(paste0(getwd(), "/output"))
hubs= read.csv("output/mod_hubs_names.csv")
hubs

#build similarity/ design matrix
keydrivers= unlist(hubs) %>% unique(.) %>% .[. != ""]
mtx= matrix(nrow= ncol(hubs), ncol = length( keydrivers)) #number of all key drivers 
row.names(mtx)= paste0(names(hubs), " mod" , "\n" ,c("/after_7d", "/before_opr", "/after_2d"))
colnames(mtx)= paste0(keydrivers)

#build similarity matrix 
#colnames of matrix included in keydrivers specified for certain module/ phenotype(rows) then put in 1
mod= apply(hubs, 2, function(x) as.list(x))

for (i in seq_along(mod)){
  for (j in 1:ncol(mtx)) {
    if ( colnames(mtx)[j] %in% mod[[i]] ) {
      mtx[i, j] <- 1
    } else {
      mtx[i, j] <- 0
    }
  }
}

library(circlize)

png("plots/chordplot_hubproteins.png", width = 11000, height = 11000, res= 600)

par(cex = .8, mar = c(0, 0, 0, 0))
circos.par(gap.degree = 1, 
           track.margin = c(0.05, 0.05), 
           points.overflow.warning = FALSE
) 

chordDiagram(mtx, 
             annotationTrack = "grid",
             transparency = 0.5)

# List of labels to add an asterisk to "CETP"  "ZMYM6" "SBSN" 
labels_to_asterisk <- c("")
# Labels to color red
labels_red <- c("")

# Customize the labels to be perpendicular, add asterisks, and color specific labels red
circos.track(track.index = 1, panel.fun = function(x, y) {
  label <- CELL_META$sector.index
  # Append asterisk to all specified labels
  modified_label <-ifelse(label %in% labels_to_asterisk , paste0(label, " ***"), label )
  # Check if the label should also be colored red
  label_color <- ifelse(label %in% labels_red, "red", "black")
  
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2]*3.5, # Adjust position as needed
              modified_label, col = label_color, facing = "bending" , #"outside", 
              niceFacing = TRUE, adj = c(0.5, 0))
}, bg.border = NA)

circos.clear()
dev.off()

#-------------------------------------------------------------------------------
#Intramodular analysis: identifying genes with high GS and MM
png("mm_vs_sig.png", width = 800, height = 600) 
colorOfColumn = substring(names(datKME), 4)
par(mar = c(5, 4, 4, 2) + 0.1) 
par(mfrow = c(2, 2))
selectModules = c( "brown","blue", "turquoise", "grey")
par(mfrow = c(2, length(selectModules)/2))
for (module in selectModules) {
  column = match(module, colorOfColumn)
  restModule = moduleColorsAutomatic == module
  verboseScatterplot(datKME[restModule, column], GS.lvl[restModule],
                     xlab = paste("Module Membership ",
                                  module, "module"), ylab = "pertubation_lvl", main = paste("kME.", module,
                                                                                            "vs. Protein Sig."), col = module)
  
}             
dev.off() 


