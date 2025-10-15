#### Input paths ####
# A.limnaeuse expression data in escape trajectory
escape_exr_path <- ""
# A.limnaeuse expression data in diapause trajectory
diapause_exr_path <- ""
# Fitted model's parameters feom escape expression
escape_param_path <- ""
# Fitted model's parameters feom diapause expression
diapause_param_path <- ""

#### Upload and modify data ####
# Fitted model values. based on linear, sigmoid or dauble sigmoide models fit to the original expression data
df_e <- read.csv(escape_exr_path)
df_d <- read.csv(diapause_exr_path)
# The parameters of the fitted models
df_e_param <- read.csv(escape_param_path)
df_d_param <- read.csv(diapause_param_path)

#### Identify the pattern of the model data: linear, increaseor decrease ####
add_status <- function(df_param, df_data) {
  status_vec <- c()  
  for (i in 1:nrow(df_param)) {
    model_type <- df_param[i, 10]
    # If the model type is linear the pattern is linear
    if (model_type == 0) {
      status_vec <- c(status_vec, 0)
    # If the model type is sigmoid the pattern is increase or decrease determined by the first and last points
    } else if (model_type == 1) {
      start  <- df_data[i,2]
      end <- df_data[i,8] 
      if (start < end) {
        status_vec <- c(status_vec, 1)   # increasing
      } else {
        status_vec <- c(status_vec, -1)  # decreasing
      }
    # If the model type is dauble sigmoid the pattern is increase or decrease determined by the order of the minimal and maximal point
    } else if (model_type == 2) {
      max_i <- which.max(df_data[i, -1])
      min_i <- which.min(df_data[i, -1])
      if (max_i > min_i) {
        status_vec <- c(status_vec, 1)
      } else {
        status_vec <- c(status_vec, -1)
      }
    }
  }
  df_data$status <- status_vec
  return(df_data)
}
# Execute the add status function
df_e <- add_status(df_e_param, df_e)
df_d <- add_status(df_d_param, df_d)

#### Naming the clusters ####
clusters <- c("Hot up_Cold up","Hot up_Cold constant","Hot up_Cold down",
              "Hot constant_Cold up","Hot constant_Cold constant","Hot constant_Cold down",
              "Hot down_Cold up","Hot down_Cold constant","Hot down_Cold down")

#### Plotting the clusters ####
# Open file
pdf(paste("ALM pattern clusters split by max gap",".pdf"),width = 28,height = 10)
# Get latout (2x6)
layout(
  matrix(c(
    1,  3,  4,  6, 8, 11,
    2,  5,  9, 7,10,  12
  ), nrow = 2, byrow = TRUE),)
counter <- 1
# Single cluster plot function
plot_cluster <- function(df_d,df_e,with_polygon,normalize_method,title){
  # Mean vectors
  d_vec <- colMeans(df_d[, c(2:8)])
  e_vec <- colMeans(df_e[, c(2:8)])
  # Standard deviation vectors
  d_sd <- apply(df_d[, c(2:8)], 2, sd)
  e_sd <- apply(df_e[, c(2:8)], 2, sd)  
  # adjust starting point
  if (normalize_method == "dif starts"){
    common_mean <-  mean(c(d_vec[1],e_vec[1]))
    d_vec <- d_vec - common_mean
    e_vec <- e_vec - common_mean
  } else{
    d_vec <- d_vec - d_vec[1]
    e_vec <- e_vec - e_vec[1]
  }
  x_labels <- c("CD", "NK", "6S", "10S", "16S", "20S", "24S")
  x_vals <- 1:length(x_labels)
  p <- plot(
    x = x_vals,
    y = d_vec,
    type = "l",
    col = "blue",
    xaxt = "n",        
    ylab = "Mean expression",
    ylim = range(c(-4,5)),
    main = title,
    cex.main = 2
  )
  axis(1, at = x_vals, labels = x_labels)
  # add polygon fir sd
  if (with_polygon){
  # Add shaded polygon for standard deviation (blue - diapause)
  polygon(
    c(x_vals, rev(x_vals)),
    c(d_vec + d_sd, rev(d_vec - d_sd)),
    col = rgb(0, 0, 1, alpha = 0.2),
    border = NA
  )
  # Add shaded polygon for standard deviation (red - escape)
  polygon(
    c(x_vals, rev(x_vals)),
    c(e_vec + e_sd, rev(e_vec - e_sd)),
    col = rgb(1, 0, 0, alpha = 0.2),
    border = NA
  )
  }
  # Add lines
  lines(x_vals, d_vec, type = "l",lwd=2, col = "blue")
  lines(x_vals, e_vec, type = "l", lwd=2,col = "red")
  return(p)  
}
# Execute plot_cluster function several times
# Run over 3 conditions for the escape data
for (hot_status in c(1,0,-1)){
  # Run over 3 conditions for the diapause data
  for (cold_status in c(1,0,-1)){
    # Split clusters in case the pattern is similar in both trajectories
    if (cold_status == hot_status){
      # Set subclusters by max gap
      max_gap_idx <- apply(abs(df_e[, paste0("X", 1:7)] - df_d[, paste0("X", 1:7)]), 1, which.max)
      valid_genes <- df_e$status == hot_status & df_d$status == cold_status
      rows <- which(valid_genes)
      cluster_df_e1 <- df_e[rows[df_e[rows, paste0("X", 1:7)][cbind(seq_along(rows), max_gap_idx[rows])] > 
                                 df_d[rows, paste0("X", 1:7)][cbind(seq_along(rows), max_gap_idx[rows])]], ]
      cluster_df_d1 <- df_d[rownames(cluster_df_e1), ]
      cluster_df_e2 <- df_e[rows[df_e[rows, paste0("X", 1:7)][cbind(seq_along(rows), max_gap_idx[rows])] < 
                                   df_d[rows, paste0("X", 1:7)][cbind(seq_along(rows), max_gap_idx[rows])]], ]
      cluster_df_d2 <- df_d[rownames(cluster_df_e2), ]
      # Get the genes of the cluster
      clusters_genes1 <- cluster_df_e1$gene
      clusters_genes2 <- cluster_df_e2$gene
      # Output with the genes of the cluster - higher in escape
      write.csv(clusters_genes1,paste(clusters[counter],".a (higher in escape).csv"))
      title <- paste(clusters[counter],".a ","(",length(clusters_genes1),")")
      plot_cluster(cluster_df_d1,cluster_df_e1,F,"dif starts",title)
      # Output with the genes of the cluster - higher in diapause
      write.csv(clusters_genes2,paste(clusters[counter],".b (higher in diapause).csv"))
      title <- paste(clusters[counter],".b ","(",length(clusters_genes2),")")
      plot_cluster(cluster_df_d2,cluster_df_e2,F,"dif starts",title)    
  # Get genes of cluster with distinct pattern  
    }else{
      cluster_df_e <-df_e[df_e$status == hot_status & df_d$status == cold_status,]
      cluster_df_d <- df_d[df_e$status == hot_status & df_d$status == cold_status,]
      # Get the genes of the cluster
      clusters_genes <- cluster_df_e$gene
      # Output with the genes of the cluster
      write.csv(clusters_genes,paste(clusters[counter],".csv"))
      title <- paste(clusters[counter],"(",length(clusters_genes),")")
      plot_cluster(cluster_df_d,cluster_df_e,F,"dif starts",title)
    }
    counter <- counter + 1

  }
}
dev.off()



