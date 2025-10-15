library(ggplot2)
library(reshape2)

#### Input data ####
NFZ_genes2clusters <- read.csv("G:/My Drive/Private folders/Omer Ben Dor/RNA seq project/Tools/Panels/RNA dynamic data/NFZ/Michal's models/240327 corrent version/Kmean & GO/240103correct version/240103 All genes in clusters.csv")[,c(1,3)]
NFZ_genes2clusters <- NFZ_genes2clusters[!duplicated(NFZ_genes2clusters), ]
ALM_genes2clusters <- read.csv("G:/My Drive/Private folders/Omer Ben Dor/RNA seq project/Tools/Panels/RNA dynamic data/ALM/Michal's models/clustering by pattern/Split cluster 1 and 9 by max gap point/All genes in all clusters.csv")[,c(1,3)]
constant_ALM <- read.csv("G:/My Drive/Private folders/Omer Ben Dor/RNA seq project/Tools/Panels/RNA dynamic data/ALM/Michal's models/clustering by pattern/Single model genes/All single model genes.csv")[,c(1,3)]

# Add constant genes data
#constant_ALM$Cluster[constant_ALM$Cluster == 1] <- 10
#constant_ALM$Cluster[constant_ALM$Cluster == 0] <- 11
#constant_ALM$Cluster[constant_ALM$Cluster == -1] <- 12
#ALM_genes2clusters <- rbind(ALM_genes2clusters,constant_ALM)
ALM_genes2clusters <- ALM_genes2clusters[!duplicated(ALM_genes2clusters), ]

# Set cluster's names and order
ALM_clusters_names <- c("ALM1.1","ALM1.2","ALM2","ALM3","ALM4","ALM5.1","ALM5.2","ALM6","ALM7","ALM8","ALM9.1","ALM9.2")
# 2V2
#ALM_order <-c("ALM1.1","ALM1.2","ALM2","ALM4","ALM3","ALM7","ALM5.1","ALM5.2","ALM6","ALM8","ALM9.1","ALM9.2")
# 6v6
ALM_order <-c("ALM1.1","ALM2","ALM3","ALM5.1","ALM6","ALM9.1","ALM1.2","ALM4","ALM7","ALM5.2","ALM8","ALM9.2")
ALM_n <- length(ALM_clusters_names)

NFZ_clusters_names <- c("NFZ1","NFZ2","NFZ3","NFZ4","NFZ5","NFZ6","NFZ7","NFZ8","NFZ9")
NFZ_order <- c("NFZ8","NFZ9","NFZ1","NFZ4","NFZ3","NFZ2","NFZ5","NFZ6","NFZ7")
NFZ_n <- length(NFZ_order)

#### Get mutual genes ####
NFZ_genes2clusters <- NFZ_genes2clusters[NFZ_genes2clusters$Hum %in% ALM_genes2clusters$Hum,]
ALM_genes2clusters <- ALM_genes2clusters[ALM_genes2clusters$Hum %in% NFZ_genes2clusters$Hum,]
# Output mutual genes
write.csv(x=data.frame("Hum" = union(ALM_genes2clusters$Hum,NFZ_genes2clusters$Hum)), file = "250719 All mutual genes in compere.csv", row.names = F)

# Generate empy matrix for compere
cluster_matrix <- matrix(0, nrow = ALM_n, ncol = NFZ_n)
rownames(cluster_matrix) <- ALM_clusters_names
colnames(cluster_matrix) <- NFZ_clusters_names

# Generate empy matrix for ratio
cluster_ratio_matrix <- matrix(0, nrow = ALM_n, ncol = NFZ_n)
rownames(cluster_ratio_matrix) <- ALM_clusters_names
colnames(cluster_ratio_matrix) <- NFZ_clusters_names

# Clusters vec
NFZ_clusters <- seq(1,ALM_n)
ALM_clusters <- c(1.1,1.2,2,3,4,5.1,5.2,6,7,8,9.1,9.2)

#### True data matrix ####
# Loop over the clusters from ALM and NFZ to count mutual genes
for (ALM_cluster_i in seq(1,ALM_n)) {
  for (NFZ_cluster_i in seq(1,NFZ_n)) {
    genes_ALM <- ALM_genes2clusters$Hum[ALM_genes2clusters$Cluster == ALM_clusters[ALM_cluster_i]]
    genes_NFZ <- NFZ_genes2clusters$Hum[NFZ_genes2clusters$Cluster == NFZ_clusters[NFZ_cluster_i]]
    mutual_genes <- intersect(unique(genes_NFZ), unique(genes_ALM))
    all_genes <- union(unique(genes_NFZ), unique(genes_ALM))
    # Count the number of mutual genes and store it in the matrix
    cluster_matrix[ALM_cluster_i, NFZ_cluster_i] <- length(mutual_genes)
    cluster_ratio_matrix[ALM_cluster_i, NFZ_cluster_i] <- length(mutual_genes)/length(all_genes)
    
  }
}
cluster_ratio_matrix<-round(cluster_ratio_matrix,2)


#### Permutation test - shuffle ALM clusters ####
set.seed(123)  # for reproducibility

n_perm <- 1000  # number of permutations
perm_results <- array(0, dim = c(ALM_n, NFZ_n, n_perm))  # store permuted ratios

# Cluster labels (your original order)
ALM_clusters <- c(1.1,1.2,2,3,4,5.1,5.2,6,7,8,9.1,9.2)
NFZ_clusters <- seq(1,NFZ_n)

# Run permutations
for (p in 1:n_perm) {
  shuffled_ALM <- ALM_genes2clusters
  shuffled_ALM$Cluster <- sample(shuffled_ALM$Cluster)  # shuffle clusters
  
  for (i in 1:ALM_n) {
    for (j in 1:NFZ_n) {
      genes_ALM <- shuffled_ALM$Hum[shuffled_ALM$Cluster == ALM_clusters[i]]
      genes_NFZ <- NFZ_genes2clusters$Hum[NFZ_genes2clusters$Cluster == NFZ_clusters[j]]
      
      mutual_genes <- intersect(unique(genes_NFZ), unique(genes_ALM))
      all_genes <- union(unique(genes_NFZ), unique(genes_ALM))
      write.csv(data.frame(Hum=mutual_genes),paste0("Mutual ALM",i," NFZ",j,".csv"))
      perm_mutual <- mutual_genes
      perm_results[i, j, p] <- length(perm_mutual)
    }
  }
}


# get p values
p_matrix <- matrix(0, nrow = ALM_n, ncol = NFZ_n)
rownames(p_matrix) <- rownames(cluster_matrix)
colnames(p_matrix) <- colnames(cluster_matrix)


for (i in 1:ALM_n) {
  for (j in 1:NFZ_n) {
    observed <- cluster_matrix[i, j]
    null_dist <- perm_results[i, j, ]
    p_value <- sum(null_dist >= observed) / n_perm
    p_matrix[i, j] <- p_value
  }
}

# correction
p_matrix_adj <- matrix(p.adjust(as.vector(p_matrix), method = "BH"), nrow = ALM_n)
rownames(p_matrix_adj) <- rownames(p_matrix)
colnames(p_matrix_adj) <- colnames(p_matrix)

# Create star labels for plotting
stars <- ifelse(p_matrix < 0.01, "**",
                ifelse(p_matrix < 0.05, "*", ""))
pval_labels <- p_matrix
pval_labels[pval_labels >= 0.05] <- NA
pval_labels[pval_labels == 0] <- 0.001

pval_labels <- -log10(pval_labels)

#### Re-order the data ####
desired_row_order <- factor(ALM_order)
desired_col_order <- factor(NFZ_order)
ordered_cluster_ratio_matrix <- cluster_ratio_matrix[desired_row_order,desired_col_order]
ordered_cluster_matrix <-cluster_matrix[desired_row_order,desired_col_order]
ordered_p_matrix_adj <- p_matrix_adj[desired_row_order, desired_col_order]
ordered_stars <- stars[match(desired_row_order, rownames(stars)),
                       match(desired_col_order, colnames(stars))]# Melt ratio matrix for plotting
ordered_pvalue <- pval_labels[match(desired_row_order, rownames(pval_labels)),
                        match(desired_col_order, colnames(pval_labels))] 

#### plotting ####
df <- melt(ordered_cluster_ratio_matrix)
df["Label"] <- as.vector(ordered_stars)
df["Adj_p_value"] <- as.numeric(ordered_pvalue)
colnames(df) <- c("Row", "Column", "Value","Label","Adj_p_value")

df_number <- melt(ordered_cluster_matrix)
colnames(df_number) <- c("Row", "Column", "Number")
df$Number <- df_number$Number
df$LabelFull <- ifelse(df$Label == "",
                       "",                     # No star → no label at all
                       paste0(df$Number))  # Show number + star

df$LabelFull <- ifelse(df$Label == "" | df$Number < 15,
                       "",
                       paste0(df$Number))
df$Adj_p_value[df$Number < 15] <- NA

# White-out cells with no stars
df[df$LabelFull=="", "Value"] <- NA

# Count genes per cluster
ALM_cluster_counts <- table(ALM_genes2clusters$Cluster)
NFZ_cluster_counts <- table(NFZ_genes2clusters$Cluster)

# Named vectors of new cluster names with counts
named_row_labels <- sapply(desired_row_order, function(cluster) {
  cluster_num <- as.numeric(gsub("[^0-9.]", "", cluster))
  count <- ALM_cluster_counts[as.character(cluster_num)]
  paste0(cluster, " (", count, ")")
})
names(named_row_labels) <- desired_row_order

named_col_labels <- sapply(desired_col_order, function(cluster) {
  cluster_num <- as.numeric(gsub("[^0-9.]", "", cluster))
  count <- NFZ_cluster_counts[as.character(cluster_num)]
  paste0(cluster, " (", count, ")")
})
names(named_col_labels) <- desired_col_order

# Update factor levels in df$Row and df$Column
df$Row <- factor(df$Row, levels = desired_row_order, labels = named_row_labels)
df$Column <- factor(df$Column, levels = desired_col_order, labels = named_col_labels)


# Plot
g <- ggplot(df, aes(x = Column, y = Row, fill = Adj_p_value)) +
  geom_tile(color = "gray") +
  geom_text(aes(label = LabelFull), color = "black", size = 4.5) +
  # NFZ annotation
  geom_vline(xintercept = c(1.5,4.5,7.5), color = "black", linetype = "dashed", lwd = 1.5) +
  annotate("text", x = 1, y = 13.5, label = "Immediate", size = 5, color = "black",angle=45)+
  annotate("text", x = 3, y = 13.7, label = "Intermediate", size = 5, color = "black",angle=45)+
  annotate("text", x = 6, y = 13.7, label = "Later", size = 5, color = "black",angle=45)+
  annotate("text", x = 8.5, y = 13.7, label = "Decrease", size = 5, color = "black",angle=45)+
  annotate("text", x = 8.5, y = 15, label = "", size = 5, color = "black")+
  
  # ALM annotation 2v2 
  #geom_hline(yintercept = c(2.5,4.5,6.5,8.5,10.5), color = "black", linetype = "dashed", lwd = 0.5) +
  #annotate("text", x = 10.3, y = 1.5, label = "^^", size = 5, color = "black")+
  #annotate("text", x = 10.3, y = 3.5, label = "^>", size = 5, color = "black")+
  #annotate("text", x = 10.3, y = 5.5, label = "^v", size = 5, color = "black")+
  #annotate("text", x = 10.3, y = 7.5, label = ">>", size = 5, color = "black")+
  #annotate("text", x = 10.3, y = 9.5, label = ">v", size = 5, color = "black")+
  #annotate("text", x = 10.3, y = 11.5, label = "vv", size = 5, color = "black")+
  #annotate("text", x = 10.6, y = 9.5, label = "", size = 5, color = "black")+
  
  # ALM annotation 6v6 
  geom_hline(yintercept = c(6.5), color = "black", linetype = "dashed", lwd = 1.5) +
  annotate("text", x = 10.3, y = 3.5, label = "Hot > Cold", size = 5, color = "black")+
  annotate("text", x = 10.3, y = 9.5, label = "Cold > Hot", size = 5, color = "black")+
  #annotate("text", x = 10.5, y = 14, label = "Single model", size = 5, color = "black")+
  annotate("text", x = 11, y = 9.5, label = "", size = 5, color = "black")+
  
  scale_fill_gradient(name= "adj.p-value\n[-log10]",breaks = c(3,2.75,2.5,2.25,2),labels = c(">3","2.75","2.5","2.25","2"),na.value = "white", low ="#FFFDD0" , high = "red") +
  theme_minimal() +
  ggtitle("Permutation test results (ALM shuffle)")

g

ggsave("250804 ALM vs NFZ compere 6by6.pdf", plot = g, width = 11, height = 8)
