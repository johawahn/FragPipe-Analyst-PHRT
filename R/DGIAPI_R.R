#' Obtain drug interactions of a list of proteins using the DGIdb 
#' \url{https://www.dgidb.org/}
#' 
#' \code{DGIAPI_R} creates a Data.Frame of the interaction
#' \parameters List of available options are under {'FragPipe-Analyst-PHRT/local_database/drug_prediction_DGIdb'}
#' 
#' @param genes vector (REQUIRED),
#' List of genes 
#'  
#' @param interaction_sources vector,
#' List of source names to include in the result set. 
#' If this field is omitted, all sources will be included. 
#' 
#' @param interaction_types vector,
#' List of interaction types to include in the result set. 
#' If this field is omitted, all interaction types will be included. 
#' 
#' @param source_trust_levels vector,
#' List of source trust levels to include in the result set. 
#' If this field is omitted, all trust levels will be included. 
#' 
#' @param antineoplastic_only boolean,
#' A flag denoting whether or not to limit interactions to only 
#' the ones involving antineoplastic drugs.
#' 
#' @return A data.frame
#' containing the interactions found in the query with the columns:
#' "gene_name", "drug_name", "interaction_type", "source", "gene_categories" 
#' 
#' @examples:
#' example = DGIAPI_R(genes = c("FLT1","STK1","FAKE1"))
#' example = DGIAPI_R(genes = c("FLT1","STK1","FAKE1"), 
#'                  interaction_sources = c("TTD","DrugBank"),
#'                  interaction_types = c("inhibitor","activator"),
#'                   gene_categories = c("kinase","tumor suppressor"),
#'                   source_trust_levels = c("Expert curated"),
#'                   antineoplastic_only = TRUE)

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







