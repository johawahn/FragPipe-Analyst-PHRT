


expr_to_frag_input <- function(quant_matrix){
  sample_names <- colnames(quant_matrix)[2:ncol(quant_matrix)]
  colnames(quant_matrix)[-1] <- paste(colnames(quant_matrix[-1]), "Intensity") 
  quant_matrix[-1] <- apply(quant_matrix[-1], c(1, 2), function(x) 10^x)

  cores <- detectCores()-1
  
  batch_size <- 250
  num_batches <- ceiling(length(quant_matrix$ProteinName) / batch_size)
  batch_indices <- split(seq_along(quant_matrix$ProteinName), 
                         rep(1:num_batches, each = batch_size, length.out = length(quant_matrix$ProteinName)))
  
  # Define a function to process a batch of data
  process_batch <- function(indices) {
    batch_data <- quant_matrix$ProteinName[indices]
    batch_result <- lapply(batch_data, function(x) GetNamesTaxa(x))
    return(batch_result)
  }
  
  # Apply the function to each batch in parallel
  result_batches <- mclapply(batch_indices, process_batch, mc.cores = cores)
  
  # Combine the results from all batches
  result <- unlist(result_batches, recursive = FALSE)
  result_df <- do.call(rbind, result)[c("Entry", "Entry.Name", "Gene.Names..primary.", 
                                        "Protein.names", "Organism")]
  frag_input <- result_df %>%
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













