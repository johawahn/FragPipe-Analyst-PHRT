

DGIAPI_R <- function(genes=c(), interaction_sources=c(), interaction_types=c(), 
                     gene_categories=c(), source_trust_levels=c(), antineoplastic_only=NULL){
  
  
  variables = list(
    'genes' = genes,
    'interaction_sources'= interaction_sources,
    'gene_categories'= gene_categories,
    'interaction_types'= interaction_types,
    'source_trust_levels'= source_trust_levels
  )  
  
  for (idx in 1:length(variables)){
    if (is.null(variables[idx])){
      break 
    } 
    variables[idx] <- paste(unname(unlist(variables[idx])), collapse=",")
    
    if (grepl("\\.", variables[idx])){
      gsub("\\.", " ", variables[idx])
    }
  }
  
  if (!is.null(antineoplastic_only)){
    variables['drug_types'] <- 'antineoplastic'
  }
  
  print(variables)
  server <- 'https://dgidb.org/api/v1/interactions.json'
  response <- GET(url=server, query=variables)
  
  content <- httr::content(response, "parse")
  matches <- content$matchedTerms
  df <- data.frame()
  
  if (!is.null(matches)){
    for (match in 1:length(matches)){
      if (length(matches[[match]]$interactions) == 0){
        break
      }
      gene <- unlist(matches[[match]]['geneName'], use.names = F)
      categories <- matches[[match]]['geneCategories']
      categories <- sort(unlist(categories, use.names = F))
      joined_categories <- paste(categories, collapse=",")
      interactions_list <- matches[[match]]$interactions
      for (idx in 1:length(interactions_list)){
        interaction <- interactions_list[[idx]]
        source <- interaction$source
        drug <- interaction$drugName
        interaction_type <-  interaction$interactionType
        df <- rbind(df, c(gene, drug, interaction_type, source, joined_categories))
      }
      
    }
  }
  
  if (length(df) == 0){
    return(NULL)
  } 
  colnames(df) <- c('gene_name', 'drug_name','interaction_type', 'source','gene_categories') 
  
  return(df)
}

#example = DGIAPI_R(genes = 'FLT1')
#example = DGIAPI_R(genes = c("DKK1", "MRPS21", "CALD1", "OS9", "SLBP", "DPH2", "USP36", "SNX12", "RABGAP1"))

#example = DGIAPI_R(genes = c("DKK1", "MRPS21", "CALD1", "OS9", "SLBP", "DPH2", "USP36", "SNX12", "RABGAP1"), 
#                  interaction_sources = NULL,
#                   interaction_types = NULL,
#                   gene_categories = NULL,
#                   source_trust_levels = NULL,
#                   antineoplastic_only = NULL)






