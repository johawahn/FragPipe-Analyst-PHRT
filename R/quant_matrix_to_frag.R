


expr_to_frag_input <- function(quant_matrix){
  sample_names <- colnames(quant_matrix)[2:ncol(quant_matrix)]
  colnames(quant_matrix)[-1] <- paste(colnames(quant_matrix[-1]), "Intensity") 
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

expr_manifest_to_frag_ano <- function(df, quant_matrix){
  sample_names <- colnames(quant_matrix)[2:ncol(quant_matrix)]
  df <- df[c("Run.Label", "Condition", "Replicate", "File.Name")] %>%
    mutate(sample=paste(df$Run.Label,df$Replicate, sep="_")) %>%
    mutate(sample_name=sample) %>%
    rename(., "file"="File.Name", "condition"="Condition", "replicate"="Replicate") %>%
    select(., -Run.Label) %>%
    filter(sample %in% sample_names) 
  
  df$condition[which(df$condition==""|df$condition=="#N/A")] = "?"
  return(df)
}






