#### Packages ####
library(clusterProfiler)
library(ggplot2)
library(msigdbr)
library(readxl)
library(dplyr)

#### Input paths ####
# Get ontology
subC <- c('CP:KEGG', 'CP:BIOCARTA', 'CP:REACTOME', 'GO:BP', 'GO:MF')
#subC <- c('CP:KEGG','CP:REACTOME')
all_gene_sets = msigdbr(species = "human")
msigdbr_t2g  = all_gene_sets[all_gene_sets$gs_subcollection %in% subC, c('gs_name', 'gene_symbol')]
msigdbr_t2g  <- msigdbr_t2g[!duplicated(msigdbr_t2g), ]

# Background with all the genes that expressed in the RNA samples - genes that didn't fit to the model and constant genes
bg <- read.csv("G:/My Drive/Private folders/Omer Ben Dor/RNA seq project/Tools/Panels/RNA dynamic data/ALM/Michal's models/kmean/corrent/250626 ALM all genes in all clusters.xlsx")$Hum

# Input data
genes_of_interest <- read.csv("G:/My Drive/Private folders/Omer Ben Dor/RNA seq project/Tools/Panels/RNA dynamic data/ALM/Michal's models/kmean/corrent/Cluster 9 genes.csv")$Hum

# GO enrichment analysis 
msigdb_enrich <- enricher(
                      gene         = genes_of_interest,
                      TERM2GENE = msigdbr_t2g[, c("gs_name", "gene_symbol")],
                      universe     = bg, 
                      pAdjustMethod = "BH",   
                      qvalueCutoff  = 0.1)       
msigdb_enrich_df <- as.data.frame(msigdb_enrich)
rownames(msigdb_enrich_df) <- NULL
msigdb_enrich_df$ID <- sub("_.*", "", msigdb_enrich_df$ID)
msigdb_enrich_df$Description <- sub("^[^_]*_", "", msigdb_enrich_df$Description)
msigdb_enrich_df <- msigdb_enrich_df[order(msigdb_enrich_df$p.adjust), ]

# Get results
write.csv(msigdb_enrich_df,"ORA combined datasets cluster 9.csv")




# Plot the top enriched
if (nrow(msigdb_enrich_df) !=0){
p <- dotplot(msigdb_enrich, showCategory = 10) + 
  labs(title= "GO",x = "Gene Ratio\n[Genes in the data/Genes in onthology]") +
  theme(plot.title = element_text(hjust = 0.5))
}

