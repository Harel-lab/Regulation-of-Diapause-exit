# Drug enrichment analysis
library(enrichR)
library(webchem)
library(readxl)

##### 1.Input #####
# Load gene lists ALM
#gene_lists_folder <- "G:/My Drive/Private folders/Omer Ben Dor/RNA seq project/Tools/Panels/RNA dynamic data/ALM/Michal's models/clustering by pattern/Split cluster 1 and 9 by max gap point/"
#gene_list_files <- list.files(gene_lists_folder, pattern = "\\.csv$", full.names = TRUE)
#gene_list_files <- gene_list_files[gene_list_files != "G:/My Drive/Private folders/Omer Ben Dor/RNA seq project/Tools/Panels/RNA dynamic data/ALM/Michal's models/clustering by pattern/Split cluster 1 and 9 by max gap point/All genes in all clusters.csv"]

# Load gene list NFZ
gene_lists_folder <- "G:/My Drive/Private folders/Omer Ben Dor/RNA seq project/Tools/Panels/RNA dynamic data/NFZ/Michal's models/240327 corrent version/Kmean & GO/240103correct version/250728 Genes list for drug enrichment/"
gene_list_files <- list.files(gene_lists_folder, pattern = "\\.xlsx$", full.names = TRUE)

# Load tested drugs
provided_drugs <- read_excel("G:/My Drive/Private folders/Omer Ben Dor/RNA seq project/Tools/Panels/Chemical screen/250707 summary.xlsx", skip = 1)[,c(2,3)]
colnames(provided_drugs)[1] <- "Drug" 
# Load drug to CAS DB
drug2CAS <- read.csv("G:/My Drive/Private folders/Omer Ben Dor/RNA seq project/Tools/Panels/Chemical screen/Drugs and CAS.csv")

# Prepare a data frame to store combined enrichment results
all_enrichments <- data.frame()



#### 2.Loop over each gene list file ####
# Check CAS fit
get_cas <- function(drug_name) {
  # Check if the drug_name is already a CAS
  if (grepl("^[0-9\\-]+$", drug_name) && grepl("-", drug_name)) {
    return(drug_name)
  }
  
  # Remove drugs with special RE
  special_chars_pattern <- "[.()\\^{}+$*?\\[\\]|\\-]"
  if (grepl(special_chars_pattern, drug_name)) {
    return(NA)
  }
  
  # Check drug common name
  direct_matches <- which(tolower(drug2CAS$Common.name) == tolower(drug_name))
  if (length(direct_matches) > 0) {
    return(drug2CAS$CAS[direct_matches[1]])
  }
  
  # Check within synonymus
  pattern <- paste0("(^|\\|\\s*)", tolower(drug_name), "(\\s*\\||$)")
  synonym_matches <- which(grepl(pattern, tolower(drug2CAS$Synonyms)))
  if (length(synonym_matches) > 0) {
    return(drug2CAS$CAS[synonym_matches[1]])
  }
  
  # Fail to match
  return(NA)
}
for (file in gene_list_files) {
  gene_list <- read_excel(file)$Hum
  # Enrichr drug enrichment for this gene list
  enrich_results <- enrichr(gene_list, databases = "DSigDB")
  significant_hits <- enrich_results[["DSigDB"]][
    enrich_results[["DSigDB"]]$Adjusted.P.value < 0.1, ] # set p-value cut-off
  if (nrow(significant_hits) == 0) next
  # Clean drug names from Term field
  significant_hits$Term <- sub(" CTD.*","",significant_hits$Term) 
  significant_hits$Term <- sub(" BOSS.*","",significant_hits$Term) 
  significant_hits$Term <- sub(" PC3.*","",significant_hits$Term) 
  significant_hits$Term <- sub(" TTD.*","",significant_hits$Term) 
  significant_hits$Term <- sub(" DOWN.*","",significant_hits$Term) 
  significant_hits$Term <- sub(" MCF7.*","",significant_hits$Term) 
  significant_hits$Term <- sub(" ssMCF7.*","",significant_hits$Term)
  significant_hits$Term <- sub(" HL60.*","",significant_hits$Term) 
  
  # Add CAS number
  significant_hits$CAS <- sapply(significant_hits$Term, get_cas)
  # Add gene list source column for traceability
  significant_hits$GeneListFile <- basename(file)
  # Append to combined data frame
  all_enrichments <- rbind(all_enrichments, significant_hits)
  print(paste("Done....",file))
}
# Remove NAs and duplicates
all_enrichments <- all_enrichments[!is.na(all_enrichments$CAS), ]

#### 3.Check which provided drugs are enriched ####
provided_drugs$IsEnriched <- provided_drugs$CAS %in% all_enrichments$CAS
print(paste("True:",sum(provided_drugs$IsEnriched),length(provided_drugs$Drug),round(sum(provided_drugs$IsEnriched)/length(provided_drugs$Drug)*100,2),"%"))
#### 4.Save output files ####
# Tested drug
write.csv(provided_drugs, "provided_drugs_with_enrichment_flag NFZ.csv", row.names = FALSE)
# Enrichment results
write.csv(all_enrichments,"NFZ Drug enrichment results.csv",row.names = F)
