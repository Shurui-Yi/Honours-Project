library(scran)
library(igraph)
library(BiocNeighbors)
library(SingleCellExperiment)
library(reshape2)
library(dplyr)
library(patchwork)

## Meta-Cell Aggregation
## Method 1: Overlapped
## data aggregation
get_mc_overlapped <- function(sc_L1_log_103g,k=10){
  # find k nearest neighbours
  cells_knn10_idx <- findKNN(t(sc_L1_log_103g), k) # to cells by genes
  cells_knn10_idx <- as.data.frame(cells_knn10_idx)
  
  #randomly choose 1000 cols as the center of meta cell
  set.seed(2022)
  sample_col <- sample(ncol(sc_L1_log_103g), 1000)
  
  #get average 
  
  mc <- get_avg(cells_knn10_idx[sample_col,],sc_L1_log_103g,sc_L1_log_103g[,sample_col],k)
  
  return(mc)# g by mc
}

get_avg <- function(idx_1000,mat_log,temp_log,k){
  for (g in 1:102){
    for (c in 1:1000){
      
      cells_idx <-  as.numeric(idx_1000[c,1:k])
      #get data of selected cells
      values <- as.numeric(mat_log[g,cells_idx])
      avg <- mean(values)
      
      temp_log[g,c] <- avg
    }
    
  } 
  return(temp_log)
}

# get the meta cells (disjoint groups)
get_mc_disjoint <- function(sc_L1_log_103g,k=3){
  #build graph and clustering
  g <- buildKNNGraph(sc_L1_log_103g,k)
  clusters_knn_w_3 <- igraph::cluster_walktrap(g)$membership #clustering of 4091 cells
  
  temp <- as.data.frame(t(sc_L1_log_103g)) #r:cell, c:gene
  temp$cluster <- clusters_knn_w_3 
  sc_L1_log_103g_avg <- aggregate(temp,list(temp$cluster),FUN= mean)
  sc_L1_log_103g_avg <- subset(sc_L1_log_103g_avg, select=-c(cluster))
  sc_L1_log_103g_avg <- sc_L1_log_103g_avg[,-1]
  return(t(sc_L1_log_103g_avg))#gene by mc
}

get_mc_disjoint_Louvain <- function(sc_L1_log_103g,k=3){
  #build graph and clustering
  g <- buildKNNGraph(sc_L1_log_103g,k)
  clusters_knn_w_3 <- igraph::cluster_louvain(g)$membership #clustering of 4091 cells
  
  temp <- as.data.frame(t(sc_L1_log_103g)) #r:cell, c:gene
  temp$cluster <- clusters_knn_w_3 
  sc_L1_log_103g_avg <- aggregate(temp,list(temp$cluster),FUN= mean)
  sc_L1_log_103g_avg <- subset(sc_L1_log_103g_avg, select=-c(cluster))
  sc_L1_log_103g_avg <- sc_L1_log_103g_avg[,-1]
  return(t(sc_L1_log_103g_avg))#gene by mc
}
