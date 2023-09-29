#' Quantitative matrix to FragPipe input
#' conversion using local UniProt human proteome database
#'
#' \code{expr_to_frag_input} creates a FragPipe-type output file
#' based on a protein intensity matrix. 
#' 
#' @param quant_matrix Data.frame,
#' Protein intensity matrix.
#' 
#' \code{expr_manifest_to_frag_ano} creates a FragPipe-type annotation file
#' based on a Spectronaut Condition experimental design.
#' 
#' @param ano Data.frame,
#' Annotation file 
#' (output from Spectronaut).
#' 
#' @return A data.frame
#' with the necessary columns for the LFQ protein
#' level analysis


expr_to_frag_input <- function(quant_matrix){
  
  sample_names <- colnames(quant_matrix)[2:ncol(quant_matrix)]
  colnames(quant_matrix)[-1] <- paste(colnames(quant_matrix[-1]), "Intensity")
  #colnames(quant_matrix)[1] <- "ProteinName"
  quant_matrix[-1] <- apply(quant_matrix[-1], c(1, 2), function(x) 10^x)
  
  info_proteins <- read.csv("local_database/uniprotkb_proteome_UP000005640_2023_09_06.tsv", sep='\t')
  
  frag_input <- info_proteins %>%
    mutate(Protein=paste(Entry, Entry.Name, sep='|'), `Combined Total Peptides`=NA) %>%
    rename(., "ProteinName"="Entry") %>%
    merge(., quant_matrix, by="ProteinName") %>%
    rename(., "Protein ID"="ProteinName", "Entry Name"="Entry.Name", 
           "Gene"="Gene.Names..primary.", "Description"="Protein.names")
  
  return(frag_input)
}

expr_manifest_to_frag_ano <- function(ano){
  ano <- ano[c("Run.Label", "Condition", "Replicate", "File.Name")] %>%
    mutate(sample=paste(ano$Run.Label,ano$Replicate, sep="_")) %>%
    mutate(sample_name=sample) %>%
    rename(., "file"="File.Name", "condition"="Condition", "replicate"="Replicate") %>%
    select(., -Run.Label) 

  ano$condition[which(ano$condition==""|ano$condition=="#N/A")] = "?"
  return(ano)
}






