#### Packages ####
library(enrichR)
library(webchem)
library(readxl)

#### input paths ####
# Folder with all the gene lists that need to be analysed
gene_lists_folder <-
# Excel format - Second column = drug name. Theird column = CAS number
Tested_drugs_path <-
# Table that fit drug to cas number
Drug2CAS_path <- 
Output_title <-

##### Load and modify data #####
gene_list_files <- list.files(gene_lists_folder, pattern = "\\.xlsx$", full.names = TRUE)
provided_drugs <- read_excel(Tested_drugs_path, skip = 1)[,c(2,3)]
colnames(provided_drugs)[1] <- "Drug" 
drug2CAS <- read.csv(Drug2CAS_path)

#### Loop over each gene list file ####
# Find the CAS number of given drug name
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
# Prepare a data frame to store combined enrichment results
all_enrichments <- data.frame()
# Execute drug enrichment analysis per provided file
for (file in gene_list_files) {
  gene_list <- read_excel(file)$Hum
  # Enrichr drug enrichment 
  enrich_results <- enrichr(gene_list, databases = "DSigDB")
  significant_hits <- enrich_results[["DSigDB"]][
    enrich_results[["DSigDB"]]$Adjusted.P.value < 0.1, ] # set p-value cut-off
  if (nrow(significant_hits) == 0) next
  # Clean DB names from Term field
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
  # Add the source file 
  significant_hits$GeneListFile <- basename(file)
  # Append to combined data frame
  all_enrichments <- rbind(all_enrichments, significant_hits)
}
# Remove NAs and duplicates
all_enrichments <- all_enrichments[!is.na(all_enrichments$CAS), ]

#### Check which provided drugs are enriched ####
provided_drugs$IsEnriched <- provided_drugs$CAS %in% all_enrichments$CAS
print(paste("True:",sum(provided_drugs$IsEnriched),length(provided_drugs$Drug),round(sum(provided_drugs$IsEnriched)/length(provided_drugs$Drug)*100,2),"%"))

#### Output files ####
# Tested drug
write.csv(provided_drugs, "provided_drugs_with_enrichment_flag.csv", row.names = FALSE)
# Enrichment results
write.csv(all_enrichments,paste0(Output_title,".csv"),row.names = F)

