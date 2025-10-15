# input
df_e <- read.csv("G:/My Drive/Private folders/Omer Ben Dor/RNA seq project/Tools/Panels/RNA dynamic data/ALM/Michal's models/240401 ALM singe&double_sigmoid_fitted_data_E.csv")
df_d <- read.csv("G:/My Drive/Private folders/Omer Ben Dor/RNA seq project/Tools/Panels/RNA dynamic data/ALM/Michal's models/240401 ALM singe&double_sigmoid_fitted_data_D.csv")


df_e_param <- read.csv("G:/My Drive/Private folders/Omer Ben Dor/RNA seq project/Tools/Panels/RNA dynamic data/ALM/Michal's models/240327 ALM model_fit_param_D.csv")
df_d_param <- read.csv("G:/My Drive/Private folders/Omer Ben Dor/RNA seq project/Tools/Panels/RNA dynamic data/ALM/Michal's models/240327 ALM model_fit_param_E.csv")

# calc status
add_status <- function(df_param, df_data) {
  status_vec <- c()
  
  for (i in 1:nrow(df_param)) {
    model_type <- df_param[i, 10]
    
    if (model_type == 0) {
      status_vec <- c(status_vec, 0)
      
    } else if (model_type == 1) {
      start  <- df_data[i,2]
      end <- df_data[i,8]
      
      if (start < end) {
        status_vec <- c(status_vec, 1)   # increasing
      } else {
        status_vec <- c(status_vec, -1)  # decreasing
      }
      
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
# Run the add status function
df_e <- add_status(df_e_param, df_e)
df_d <- add_status(df_d_param, df_d)

# spesific clusters
clusters <- c("Hot up_Cold up","Hot up_Cold constant","Hot up_Cold down",
              "Hot constant_Cold up","Hot constant_Cold constant","Hot constant_Cold down",
              "Hot down_Cold up","Hot down_Cold constant","Hot down_Cold down")

# Open file
pdf(paste("ALM pattern clusters split by max gap",".pdf"),width = 28,height = 10)

layout(
  matrix(c(
    1,  2,  3,  4,
    5,  6,  7,  8,
    9,  10,  11,  12
  ), nrow = 3, byrow = TRUE),)

layout(
  matrix(c(
    1,  3,  4,  6, 8, 11,
    2,  5,  9, 7,10,  12
  ), nrow = 2, byrow = TRUE),)

counter <- 1

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

for (hot_status in c(1,0,-1)){
  for (cold_status in c(1,0,-1)){
    if (cold_status == hot_status){
      # Sety subclusters by max gap
      max_gap_idx <- apply(abs(df_e[, paste0("X", 1:7)] - df_d[, paste0("X", 1:7)]), 1, which.max)
      valid_genes <- df_e$status == hot_status & df_d$status == cold_status
      rows <- which(valid_genes)
      cluster_df_e1 <- df_e[rows[df_e[rows, paste0("X", 1:7)][cbind(seq_along(rows), max_gap_idx[rows])] > 
                                 df_d[rows, paste0("X", 1:7)][cbind(seq_along(rows), max_gap_idx[rows])]], ]
      cluster_df_d1 <- df_d[rownames(cluster_df_e1), ]
      cluster_df_e2 <- df_e[rows[df_e[rows, paste0("X", 1:7)][cbind(seq_along(rows), max_gap_idx[rows])] < 
                                   df_d[rows, paste0("X", 1:7)][cbind(seq_along(rows), max_gap_idx[rows])]], ]
      cluster_df_d2 <- df_d[rownames(cluster_df_e2), ]
      
      
      # set sub-clusters by one point
      #cluster_df_e1 <-df_e[df_e$status == hot_status & df_d$status == cold_status & df_e$X7 > df_d$X7,]
      #cluster_df_d1 <- df_d[df_e$status == hot_status & df_d$status == cold_status & df_e$X7 > df_d$X7,]
      #cluster_df_e2 <-df_e[df_e$status == hot_status & df_d$status == cold_status & df_e$X7 < df_d$X7,]
      #cluster_df_d2 <- df_d[df_e$status == hot_status & df_d$status == cold_status & df_e$X7 < df_d$X7,]
      
      clusters_genes1 <- cluster_df_e1$gene
      clusters_genes2 <- cluster_df_e2$gene
      
      write.csv(clusters_genes1,paste(clusters[counter],".a (higher in escape).csv"))
      title <- paste(clusters[counter],".a ","(",length(clusters_genes1),")")
      plot_cluster(cluster_df_d1,cluster_df_e1,F,"dif starts",title)
      
      write.csv(clusters_genes2,paste(clusters[counter],".b (higher in diapause).csv"))
      title <- paste(clusters[counter],".b ","(",length(clusters_genes2),")")
      plot_cluster(cluster_df_d2,cluster_df_e2,F,"dif starts",title)
      
    
    
    }else{
      cluster_df_e <-df_e[df_e$status == hot_status & df_d$status == cold_status,]
      cluster_df_d <- df_d[df_e$status == hot_status & df_d$status == cold_status,]
      clusters_genes <- cluster_df_e$gene
      print(length(df_e_param$model[df_e_param$id %in% clusters_genes]))
      write.csv(clusters_genes,paste(clusters[counter],".csv"))
      title <- paste(clusters[counter],"(",length(clusters_genes),")")
      plot_cluster(cluster_df_d,cluster_df_e,F,"dif starts",title)
    }
    counter <- counter + 1

  }
}
dev.off()


