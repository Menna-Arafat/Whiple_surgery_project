setwd("C:/Users/USER/Documents/pancreas/ppi")



#BiocManager::install("estrogen")
# library(affy)
# library(estrogen)
# library(vsn)
# library(genefilter)
library(igraph)
library(tidyr)
library(readxl)
library(openxlsx)
library(scales)
library(ggplot2)
library(tibble)
library(plyr)
library(dplyr)

#load data
list.files()
for (file in list.files(pattern= "^string")){
  net= read.delim(file)[1:2] %>% as.data.frame()
  names(net)= c("node1", "node2")
  net= net[!apply(net, 1, function(row) any(grepl("^(ENS|LOC)", row))), ]
  name= strsplit(file, "_")[[1]][5] %>% gsub(".tsv", "_ppi", .)
  assign(name, net, envir=.GlobalEnv)
}


for (file in list.files(pattern= "^jaspar")){
  net= read.delim(file, sep= ",")[,c(1,3)] %>% as.data.frame()
  names(net)= c("node1", "node2")
  net= net[!apply(net, 1, function(row) any(grepl("^(ENS|LOC)", row))), ]
  name= strsplit(file, "_")[[1]][2] %>% gsub(".csv", "_tf", .)
  assign(name, net, envir=.GlobalEnv)
}

com_net= rbind(grey_tf, grey_ppi)

hubs= read.csv("module_hubs_names.csv")

modules= read.csv( "modules_gene_name.csv" )
#--------------------------------
#build igraph object from edge list
graph =graph_from_edgelist(com_net %>%as.matrix(), directed = FALSE)

#centrality measures
# Degree Centrality, edges cconnected to a node
degree_centrality <- igraph::degree(graph)  %>% sort(.,decreasing = T) 
print(names(degree_centrality )[1:5])
DC= names(degree_centrality)[1:5] %>% as.vector()

# Betweenness Centrality, The number of shortest paths that pass through a vertex
betweenness_centrality <- betweenness(graph, normalized = TRUE)  %>% sort(.,decreasing = T) 
print(names(betweenness_centrality)[1:5])
BC= names(betweenness_centrality)[1:5] %>% as.vector()

# Eigenvector Centrality, The importance of a node based on the importance of its neighbors.
eigenvector_centrality <- eigen_centrality(graph)$vector %>% sort(.,decreasing = T) 
print(eigenvector_centrality)
eigc= names(eigenvector_centrality)[1:5] %>% as.vector()

#closeness centrality
closeness_centrality= closeness(graph) %>% sort(.,decreasing = T)
closeness_centrality
cc= names(closeness_centrality)[1:5] %>% as.vector()

# PageRank,Nodes are considered important if they are linked to by other important nodes, with a damping factor accounting for random jumps.
pagerank <- page_rank(graph)$vector %>% sort(.,decreasing = T) 
print(pagerank)
pr= names(pagerank)[1:5] %>% as.vector()

all_drivers= c(DC, BC, eigc, cc, pr, hubs$turquoise ) %>% unique() %>% .[.!=""]
all_drivers

#visualization
col <- colorRampPalette(c("#F6E8C3", "#1FA187FF")) 
col(10)
# Create a color vector, defaulting to black
node_colors <- rep("gray80", vcount(graph))

# Set the color of specific nodes to red
node_colors[V(graph)$name %in% all_drivers] <- "orange"
node_colors[V(graph)$name %in% tf] <- "#1FA187FF"
node_colors
#change label position
radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}
n= vcount(graph)
lab.locs <- radian.rescale(x=1:n, direction=-1, start=0)


png("igraph_ppi_tf_blue.png", width=20000, height=20000, res= 1000)
plot(graph, edge.arrow.size=.5, vertex.color= node_colors, vertex.size=5, 
     vertex.frame.color="gray", vertex.label.color="black", 
     vertex.label.cex=1, vertex.label.dist= 1.2, vertex.label.degree=lab.locs,
     layout=layout_in_circle, #layout_with_fr
     main="Integrated PPI and TF Networks for Blue Module") 

dev.off()
#----------------------------------------------------------------------------
graph_blue= graph
graph_turquoise= graph
graph_grey= graph

all_drivers_grey= all_drivers
all_drivers_turquoise= all_drivers 
all_drivers_blue= all_drivers

graph= graph_grey
all_drivers= all_drivers_grey
#---------------------------------------------------------------------------
#Cluster modularity (louvain)
clusterlouvain <- cluster_louvain(graph)
clusterleiden = cluster_leiden(graph, objective_function = "modularity", n_iteration=10)

#see the distance between clustering of two algorithms whaere small distance indicate hat those algorithms gave more similar results
compare(clusterlouvain, clusterleiden)

group1 <- V(graph)$name[which(clusterlouvain$membership == 1)]
#export communities
communities= communities(clusterlouvain)
communities= ldply(communities, rbind) %>% t()
communities[is.na(communities)]= ""
communities= communities %>% as.data.frame()
names(communities)= c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5")
communities
#write.xlsx(communities, "communities_louvain_grey.xlsx", rowNames= F)

#----------------------
#to draw communities splitted
coGrph <- delete_edges(graph, E(graph)[igraph::crossing(clusterlouvain, graph)])
plot(coGrph)
?crossing
#Stretch = 20
col= grDevices::adjustcolor(c("pink", "#B0E0E6", "#BF812D", "#46337EFF", "tan"), alpha=0.3) 
node_col= rep("lightgrey", length(V(graph)$name)) 
#node_col[ V(graph)$name  %in% blue_tf$node1 ]= "#1FA187FF"
#node_col[ V(graph)$name %in% modules$blue] = "#7B68EE"
node_col[  V(graph)$name %in% all_drivers] = "#FF4"
node_col
names(igraph:::.igraph.shapes)
node_shape= rep("circle" , length(V(graph)$name))
node_shape[ V(graph)$name %in% grey_tf$node1]= "square"
node_shape= as.character(node_shape)
node_shape
#annotation that result from enrichment analysis 
terms= read.xlsx("communities_louvain_grey.xlsx", sheet= "term")
#add text
layout <- layout_with_kk(coGrph)
# Calculate centroids for each community
centroids <- sapply(1:4, function(i) {
  colMeans(layout[clusterlouvain$membership == i, ])
}) %>% t()

#blue
png("communities_ppi_tf_blue.png", width=19000, height=18000, res=1000)
plot(graph, layout=layout_with_kk(coGrph) ,
     vertex.color= node_col, vertex.size=5,  #vertex.color= clusterlouvain$membership
     vertex.label.color="black", vertex.shape=node_shape,
     vertex.label.cex=1.2, vertex.label.dist= .9,
     mark.groups = communities(clusterlouvain), #list(group1,group2, group3, group4),
     mark.shape=1, 
     mark.col=col,
     mark.border= col,
     mark.expand = 30)

text(x= c(centroids[,1][1]*.9, centroids[,1][2]*.12, centroids[,1][3]*.12, centroids[,1][3]*.13),
     y= c(centroids[,2][1]*.28, centroids[,2][2]*.22, centroids[,2][3]*.32, centroids[,2][1]*.05),
     labels = terms$term, col = "black", cex = 1.5)

legend(c(centroids[,1][3]*.1, centroids[,2][1]*.22), legend = c("Proteins"), 
       pch = c(21),                 # 21 for circle, 24 for triangle, 22 for square
       pt.bg = "lightgrey",
       col = c("black", "black", "black"),  
       pt.cex = 5,                        
       bty = "n",                          
       title = "")
legend(c(centroids[,1][3]*.1, centroids[,2][1]*.24), legend = c( "TF"), 
       pch = c(22),                 # 21 for circle, 24 for triangle, 22 for square
       pt.bg = "lightgrey",
       col = c("black", "black", "black"),  
       pt.cex = 5,                        
       bty = "n",                          
       title = "")
legend(c(centroids[,1][3]*.1, centroids[,2][1]*.25), legend = "Hubproteins", 
       pch = 21,                    
       pt.bg = "#FF4",              # Background color of the circle
       col = "black",             
       pt.cex = 5,                 
       bty = "n")
dev.off()

#grey
#annotation that result from enrichment analysis 
terms= read.xlsx("communities_louvain_grey.xlsx", sheet= "term")

png("communities_ppi_tf_grey.png", width=19000, height=19000, res=1000)
plot(graph, layout=layout_with_kk(coGrph) ,
     vertex.color= node_col, vertex.size=5,  #vertex.color= clusterlouvain$membership
     vertex.label.color="black", vertex.shape=node_shape,
     vertex.label.cex=1.2, vertex.label.dist= .9,
     mark.groups = communities(clusterlouvain), #list(group1,group2, group3, group4),
     mark.shape=1, 
     mark.col=col,
     mark.border= col,
     mark.expand = 30)

text(x= c( centroids[,1][2]*1.7, centroids[,1][2]*1.6,  centroids[,1][3]*.14, centroids[,1][3]*.13, centroids[,1][3]*.1),
     y= c(centroids[,2][2]*.04, centroids[,2][2]*.27, centroids[,2][2]*.17, centroids[,2][2]*.01, centroids[,2][3]*.24),
     labels = terms$term, col = "black", cex = 1.5)

legend(c(centroids[,1][3]*.1, centroids[,2][2]*.26), legend = c("Proteins"), 
       pch = c(21),                 # 21 for circle, 24 for triangle, 22 for square
       pt.bg = "lightgrey",
       col = c("black", "black", "black"),  
       pt.cex = 5,                        
       bty = "n",                          
       title = "")
legend(c(centroids[,1][3]*.1, centroids[,2][2]*.23), legend = c( "TF"), 
       pch = c(22),                 # 21 for circle, 24 for triangle, 22 for square
       pt.bg = "lightgrey",
       col = c("black", "black", "black"),  
       pt.cex = 5,                        
       bty = "n",                          
       title = "")
legend(c(centroids[,1][3]*.1, centroids[,2][2]*.2), legend = "Hubproteins", 
       pch = 21,                    
       pt.bg = "#FF4",              # Background color of the circle
       col = "black",             
       pt.cex = 5,                 
       bty = "n")
dev.off()


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#adjacency
graph[]
get.adjacency(net, attr="weight", sparse=F)
#all vertices
V(graph)$name
#all edges
E(graph)
#get neighbours
n= neighborhood(graph, order=1, V(graph)$name %in% c( "APCS", "SERPINA3","A1BG" ))

names(unlist(n))    

#get edges between certain node and another node, and their count/length         
E(graph)[ V(graph)[name=="SERPINA3" ] %--% V(graph)[name=="A1BG"] ] %>% length()   
