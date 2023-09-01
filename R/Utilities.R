#' The function is to show the DE genes in a dotplot. The color of the dots represent the gene expression level. The size of the dots
#' represents percentage of cells expressing a given transcript for the clusters. 
#' @param rna_seurat A Seurat object of scRNA-seq data after normalization.
#' @param gene_list A list of genes users want to show in the dotplot. 
#' @return A data frame to prepare for dotplot visualization.
do_dotplot = function(rna_seurat,gene_list){
  # Load libraries
  library(tidyr)
  library(dplyr)
  library(Seurat)
  library(plyr)
  
  feature_rna_plot = rna_seurat@meta.data
  gene_df = as.data.frame(t(as.matrix(rna_seurat@assays$RNA@data[(rownames(rna_seurat@assays$RNA@data) %in% gene_list),])))
  
  feature_rna_plot = cbind(feature_rna_plot,gene_df)
  
  # Calculate ave_expressing values for each cluster separately and pct for each cluster separately
  cluster_gene_df = matrix(NA,nrow = length(unique(feature_rna_plot$seurat_clusters)),ncol = 2*ncol(gene_df))
  clusters = c(1:22)
  
  for(i in 1:length(clusters)){
    subcluster_df = subset(feature_rna_plot,feature_rna_plot$seurat_clusters == clusters[i])
    for(j in 7:ncol(feature_rna_plot)){
      cluster_gene_df[i,j-6] = sum(subcluster_df[,j])/nrow(subcluster_df)
      cluster_gene_df[i,(j-6+ncol(gene_df))] = length(which(subcluster_df[,j] != 0))/nrow(subcluster_df) * 100
    }
  }
  
  cluster_gene_df = as.data.frame(cluster_gene_df)
  rownames(cluster_gene_df) = as.character(clusters)
  colnames(cluster_gene_df)[1:ncol(gene_df)] = colnames(feature_rna_plot)[7:(7+ncol(gene_df)-1)]
  colnames(cluster_gene_df)[((ncol(gene_df))+1): (2*ncol(gene_df))] = paste("pct",colnames(feature_rna_plot)[7:(7+ncol(gene_df)-1)],sep = "_")
  
  # Convert the wide format to long format
  ave_expression_df = cluster_gene_df[,1:ncol(gene_df)]
  ave_expression_df$cluster = clusters
  
  pct_df = cluster_gene_df[,(ncol(gene_df)+1): (2*ncol(gene_df))]
  pct_df$cluster = clusters
  colnames(pct_df) = colnames(ave_expression_df)
  
  return_df1 = gather(ave_expression_df,gene,ave_expression,(colnames(ave_expression_df)[1]):(colnames(ave_expression_df)[ncol(gene_df)]))
  return_df2 = gather(pct_df,gene,pct_gene,(colnames(ave_expression_df)[1]):(colnames(ave_expression_df)[ncol(gene_df)]))
  return_df = left_join(return_df1,return_df2,by = c("cluster","gene"))
  
  # return a data frame
  return(return_df)
}

#' The function is to show the DE genes in a dotplot. The color of the dots represent the gene expression level. The size of the dots
#' represents percentage of cells expressing a given transcript for the clusters. 
#' @param data A dataframe returned from the `do_dotplot` function.
#' @param gene_list A list of genes users want to show in the dotplot. 
#' @return A dotplot.
make_dotplot = function(data,gene_list,reverse = FALSE,cluster_anno = NULL,label_angle = 90, font_size = 20, aspect_ratio = 1){
  library(ggplot2)
  
  mid = mean(data$ave_expression)
  data$gene = factor(data$gene,levels = gene_list)
  
  data$cluster = factor(data$cluster,levels = c("3","4","14","5","15","7","8","10","13","16","21","2","19","11","6","1","9","17","22","18","12","20"))
  g1 = ggplot(data) + geom_point(aes(x = gene,y = as.factor(cluster),color = ave_expression,size = pct_gene)) +
    scale_color_gradient2(low = "white",high = "#E65948",mid = "#FEF6B1",midpoint = mid) +
    theme_classic() + 
    theme(axis.text = element_text(size = font_size),axis.text.x = element_text(angle = label_angle)) +
    theme(axis.line.x  = element_line(size = 0.8),axis.line.y  = element_line(size = 0.8),axis.ticks = element_line(size = 0.8)) +
    theme(axis.title = element_blank()) +
    scale_y_discrete(limits = rev(levels(data$cluster))) +
    theme(legend.title = element_blank(),legend.text = element_text(size = 20), legend.position = "right") +
    labs(color = "Norm.ave.exp",size = "Pct.expression")
  return(g1)
}

#' cosine_distance() function is to calculate cosine distance between two vectors
#' @param X1 Matrix1 with n obersavation on rows where columns represent features
#' @param X2 Matrix2 with n obersavation on rows where columns represent features
#' @param index Indicator values to locate the number of rows for both `X1` and `X2`. The `index`
#' is a vector containing two values.
#' @return The function returns 1 element, which is a similarity matrix calculated by 
#' 1-cosine_distance. The higher value, the more similarity between two cells. 
cosine_distance = function(X1, X2, index){
  if(index[1] == 1){
    print(index[2])
  }
  
  A = X1[index[1],]
  B = X2[index[2],]
  
  # cosine distance is 1 - similarity
  cos_dis = 1 - sum(A*B)/sqrt(sum(A^2)*sum(B^2))
  
  return(cos_dis)
}

#' calculate_cell_similarity() function is to calculate pairwise cell similarity
#' @param data1 Data frame 1 or matrix1 with n obersavation on rows where columns represent features
#' @param data2 Data frame 2 or matrix2 with n obersavation on rows where columns represent features
#' @param method The method to calculate pairwise cell similarity. There are a few options to choose
#' `cosine` represents cosine distance and `Euclidean` indicates Euclidean distance. 
#' @return The function returns 1 element, which is a matrix of pairwise cell similarity. 
calculate_cell_similarity = function(data1,data2,method = c("cosine","Euclidean")){
  library(pdist)
  
  c1 = nrow(data1)
  c2 = nrow(data2)
  
  # Calculate cosine distance
  if(method == "cosine"){
    cmb = expand.grid(i = 1:c1,j = 1:c2)
    result_matrix = matrix(apply(cmb,1, function(x) cosine_distance(data1, data2, x)), c1, c2)
  }
  # Calculate Eucledian distance
  if(method == "Euclidean"){
    euclidean_matrix = pdist::pdist(data1, data2)
    result_matrix = as.matrix(euclidean_matrix)
  }
  
  return(result_matrix)
}


#' distance_to_similarity() function is to convert
#' @param distance_matrix Matrix calculated from calculate_cell_similarity() function. Each element is the distance_
#' matrix is the distance between two datasets
#' @param method The method to convert distance matrix to similarity matrix. `inverse` and `gaussian_kernel` are for
#' Euclidean distance. `inverse` take the inverse of distance with 1 to avoid divided by 0. `gaussian_kernel` has 
#' an additional parameter `gaussian_sigma` to control the band width. `cosine` similarity if for cosine distance, 
#' which is basically 1-cosine_distance.
#' @param gaussian_sigma The band width of gaussian kernel. By default, `gaussian_sigma` is 1
#' @return The function returns 1 element as matrix format, which is the similarity matrix between two datasets
distance_to_similarity = function(distance_matrix, method = c("inverse","cosine","gaussian_kernel"),gaussian_sigma = 1){
  if(method == "cosine"){
    similarity_matrix = 1 - distance_matrix
  }
  else if(method == "inverse"){
    similarity_matrix = 1/(1 + distance_matrix)
  }
  else if(method == "gaussian_kernel"){
    similarity_matrix = exp(-distance_matrix^2/(2*gaussian_sigma^2))
  }
  
  return(similarity_matrix)
}

#' generate_meta_cell() is to generate meta cells to smooth the dataset
#' @param data A dataframe or a matrix with n obersavation on rows where columns represent features
#' @param order_index A vector of order index of cells to group the cells
#' @param chunk The number of chunk to divide for meta cells. By default, `chunk` is NULL. 
#' @param ncell The number of cells to generate meta cell. By default, `ncells` is 10.
#' @return The function returns 1 element, which is an average ordered matrix/data frame.
generate_meta_cell = function(data,order_index,chunk = NULL,ncells = 10){
  ordered_data = data[order_index,]
  nr = nrow(ordered_data)
  
  # Using the same number of cell for each chunk
  if(is.null(chunk) == TRUE){
    split_list = split(ordered_data, rep(1:ceiling(nr/ncells), each=ncells, length.out = nr))
    
    # For loop to calculate mean gene expression value for each meta cell
    ave_ordered_data = c()
    for(i in 1:length(split_list)){
      chunk_df = split_list[[i]]
      ave_chunk_df = base::colMeans(chunk_df)
      ave_ordered_data = rbind(ave_ordered_data,ave_chunk_df)
    }
  }
  # Use the a fixed chunk to generate meta cells
  else{
    ncells = nr %/% (chunk-1)
    split_list = split(ordered_data, rep(1:chunk, each = ncells,length.out = nr))
    
    # For loop to calculate mean gene expression value for each meta cell
    ave_ordered_data = c()
    for(i in 1:length(split_list)){
      chunk_df = split_list[[i]]
      ave_chunk_df = base::colMeans(chunk_df)
      ave_ordered_data = rbind(ave_ordered_data,ave_chunk_df)
    }
  }
  
  return(ave_ordered_data)
}

#' heatmap_visual() is to visualize the similarity matrix using heatmap
#' @param similarity_matrix A matrix from distance_to_similarity() function
#' @param index1 A vector includes the order of samples/cells. The length of index1 should be the same as the number 
#' of rows of `similarity_matrix`. By default, `index1` is NULL. 
#' @param index2 A vector includes the order of samples/cells. The length of index2 should be the same as the number 
#' of columns of `similarity_matrix`. By default, `index2` is NULL.
#' @param inverse A Boolean indicator to whether transpose the `similarity_matrix`. By default, `inverse` is TRUE.
#' @param scale Whether to scale the heatmap. It is a parameter from pheatmap() function. There are three options for
#' `scale`, which are row, column, and none. By default, `scale` is none, which indicates no scaling for the heatmap
#' @param annotation_rcol A dataframe that specifies the annotation for column
#' @param annotation_row A dataframe that specifies the annotation for row
#' @return The function returns 1 element, which is a heatmap plot object. 
heatmap_visual = function(similarity_matrix,index1=NULL,index2=NULL,inverse = T,scale = "none",normalize = "none", annotation_row=NULL, annotation_col=NULL, annotation_colors=NULL, annotation_legend, legend, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)){
  library(pheatmap)
  
  # Scale both rows and column
  if(normalize == "both"){
    similarity_matrix = scale(similarity_matrix)
    similarity_matrix = t(scale(t(similarity_matrix)))
  }
  else if(normalize == "none"){
    similarity_matrix = similarity_matrix
  }
  
  # Visualization based on indices
  if(is.null(index1) == TRUE & is.null(index2) == TRUE){
    if(inverse == TRUE){
      p1 = pheatmap(t(similarity_matrix),cluster_rows = F,cluster_cols = F,scale = scale, annotation_row=categories_row, annotation_col=categories_col, annotation_colors = annotation_colors, annotation_names_row = F, annotation_names_col = F, show_rownames = F, show_colnames = F, annotation_legend = annotation_legend, legend = legend, color = color)
    }
    else{
      p1 = pheatmap(similarity_matrix,cluster_rows = F,cluster_cols = F,scale = scale, annotation_row=categories_row,  annotation_col=categories_col, annotation_colors = annotation_colors, annotation_names_row = F, annotation_names_col = F, show_rownames = F, show_colnames = F, annotation_legend = annotation_legend, legend = legend, color = color)
    }
  }
  else if(is.null(index1) == FALSE & is.null(index2) == FALSE){
    if(inverse == TRUE){
      p1 = pheatmap(t(similarity_matrix[index1,index2]),cluster_rows = F,cluster_cols = F,scale = scale, annotation_row=categories_row, annotation_col=categories_col, annotation_colors = annotation_colors, annotation_names_row = F, annotation_names_col = F, show_rownames = F, show_colnames = F, annotation_legend = annotation_legend, legend = legend, color = color)
    }
    else{
      p1 = pheatmap(similarity_matrix[index1,index2],cluster_rows = F,cluster_cols = F,scale = scale, annotation_row=categories_row, annotation_col=categories_col, annotation_colors = annotation_colors, annotation_names_row = F, annotation_names_col = F, show_rownames = F, show_colnames = F, annotation_legend = annotation_legend, legend = legend, color = color)
    }
  }
  else if(is.null(index1) == TRUE & is.null(index2) == FALSE){
    if(inverse == TRUE){
      p1 = pheatmap(t(similarity_matrix[,index2]),cluster_rows = F,cluster_cols = F,scale = scale, annotation_row=categories_row, annotation_col=categories_col, annotation_colors = annotation_colors, annotation_names_row = F, annotation_names_col = F, show_rownames = F, show_colnames = F, annotation_legend = annotation_legend, legend = legend, color = color)
    }
    else{
      p1 = pheatmap(similarity_matrix[,index2],cluster_rows = F,cluster_cols = F,scale = scale, annotation_row=categories_row, annotation_col=categories_col, annotation_colors = annotation_colors, annotation_names_row = F, annotation_names_col = F, show_rownames = F, show_colnames = F, annotation_legend = annotation_legend, legend = legend, color = color)
    }
  }
  else if(is.null(index1) == FALSE & is.null(index2) == TRUE){
    if(inverse == TRUE){
      p1 = pheatmap(t(similarity_matrix[index1,]),cluster_rows = F,cluster_cols = F,scale = scale, annotation_row=categories_row, annotation_col=categories_col, annotation_colors = annotation_colors, annotation_names_row = F, annotation_names_col = F, show_rownames = F, show_colnames = F, annotation_legend = annotation_legend, legend = legend, color = color)
    }
    else{
      p1 = pheatmap(similarity_matrix[index1,],cluster_rows = F,cluster_cols = F,scale = scale, annotation_row=categories_row, annotation_col=categories_col, annotation_colors = annotation_colors, annotation_names_row = F, annotation_names_col = F, show_rownames = F, show_colnames = F, annotation_legend = annotation_legend, legend = legend, color = color)
    }
  }
  
  return(p1)
}


#' linear_mapping_function() is to learn matrix A to map two datasets
#' @param X Matrix1 with n rows (samples/cells) and p columns (features/genes)
#' @param Y Matrix1 with n rows (samples/cells) and p columns (features/genes)
#' @return The function returns 1 element, which is the mapping matrix with (p+1)-by-p dimension
linear_mapping_function = function(X,Y){
  library(pracma)
  n = nrow(X)
  p = ncol(X)
  
  X = as.matrix(X)
  Y = as.matrix(Y)
  
  # Trick
  Xtil = cbind(X,rep(1,n))
  
  # Closed form for least square problem in matrix form
  Atil = pinv(t(Xtil) %*% Xtil) %*% (t(Xtil) %*% Y)
  
  return(Atil)
}

#' add_label() function is to generate metacell state information 
#' @param data A vector that contains state information for each cell
#' @param ncells The number of cells to generate meta cell. By default, `ncells` is 10.
add_label = function(data,ncells){
  split_list = split(data, ceiling(seq_along(data)/ncells))
  
  #for loop to decide each metacell's state
  metacell_data = c()
  for(i in 1:length(split_list)){
    chunk_vec = split_list[[i]]
    metacell_data = c(metacell_data,names(sort(table(chunk_vec), decreasing = TRUE))[1])
  }
  return(metacell_data)
}

