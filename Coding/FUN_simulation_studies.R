# install.packages("R.methodsS3")
# install.packages("nnls")
# install.packages("lcmix", repos="http://R-Forge.R-project.org")
library(lcmix)#rmvgamma
library(dplyr)
library(ggplot2)
library(ggstatsplot)


##data simulation
simulate_logcounts <- function(sc_L_log_g,distribution="mv_gamma",n_cells,mean_type="all"){
  # calculate the orginal correlation
  org_Pcor <- cor(t(sc_L_log_g), use="pairwise.complete.obs",method = "pearson")
  
  if (mean_type=="all"){
    meanbar <- rowMeans(sc_L_log_g)
    variance <- rowVars(sc_L_log_g)
  }else if (mean_type=="expressed"){
    meanbar <- apply(sc_L_log_g,1, function(x) mean(x[x!=0]))
    variance <-apply(sc_L_log_g,1, function(x) var(x[x!=0]))
  }
  
  if (distribution=="mv_gamma"){
    
    shape <-  meanbar^2/variance # shape is alpha
    scale <- variance/meanbar #scale is beta
    rate <- 1/scale
    
    set.seed(123)
    simulated_gamma_data <- rmvgamma(n=n_cells, shape = shape, rate = 1/scale , corr =org_Pcor)
    
    #name the columns (genes)
    colnames(simulated_gamma_data) <- colnames(org_Pcor) #correlation is always cell by gene
    
    # change back to genes by cells
    simulated_gamma_data <- t(simulated_gamma_data)
  }
  
  return(simulated_gamma_data)
}
## add sparsity
add_zeros<- function(gene_perc,cell_perc,simulated_gamma_data){
  
  total_genes <- dim(simulated_gamma_data)[1]
  total_cells <- dim(simulated_gamma_data)[2]
  
  #select genes
  set.seed(2)
  gene_index <- sample(total_genes, total_genes*gene_perc) 
  
  set.seed(2)
  for (gene in gene_index){
    cell <- sample(total_cells, total_cells*cell_perc) 
    simulated_gamma_data[gene,cell] <- 0
  }
  return(list(gene_index=gene_index,simulated_data=simulated_gamma_data)) #return selected gene idx
}

##correlation calculation for gene by cell matrix
get_cor <- function(sc_L1_log_103g,method){
  cor <- cor(t(sc_L1_log_103g), use="pairwise.complete.obs",method=method)
  return(cor)
}

##correlation of gene pairs
get_cor_uppertri <- function(org_log_cor_P,cor_name){
  
  rowCol <- expand.grid(rownames(org_log_cor_P), colnames(org_log_cor_P))
  labs <- rowCol[as.vector(upper.tri(org_log_cor_P,diag=F)),]
  df <- cbind(labs, org_log_cor_P[upper.tri(org_log_cor_P,diag=F)])
  colnames(df) <- c("gene1","gene2",cor_name)
  
  return(df)
}

## add groups to logcounts matrix
get_cor_uppertri_type <- function(org_log_cor_P,cor_name,sparse_genes){
  
  rowCol <- expand.grid(rownames(org_log_cor_P), colnames(org_log_cor_P))
  labs <- rowCol[as.vector(upper.tri(org_log_cor_P,diag=F)),]
  df <- cbind(labs, org_log_cor_P[upper.tri(org_log_cor_P,diag=F)])
  colnames(df) <- c("gene1","gene2",cor_name)
  
  df <- df%>% mutate(type=case_when(
    gene1 %in% sparse_genes & gene2 %in% sparse_genes ~"both_added",
    gene1 %in% sparse_genes | gene2 %in% sparse_genes ~"either_added",
    TRUE ~"neither_added"
  ))
  
  return(df)
}

add_cor_gene_groups <- function(merged_cor_df,sparse_genes){
  
  merged_cor_df<- merged_cor_df%>% mutate(type=case_when(
    gene1 %in% sparse_genes & gene2 %in% sparse_genes ~"both_added",
    gene1 %in% sparse_genes | gene2 %in% sparse_genes ~"either_added",
    TRUE ~"neither_added"
  ))
  return(merged_cor_df)
}

## 
calculate_rms <- function(x_true,y_obs){
  return(sqrt(mean((y_obs-x_true)^2)))
}
