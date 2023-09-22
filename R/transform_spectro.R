#' Spectronaut protein report to FragPipe protein level report
#'
#' \code{create_annotation} creates a TMT-type output file
#'
#' @param spectro_cond Data.frame,
#' Protein level Spectronaut report
#' 
#' \code{create_quant} creates a TMT-type annotation file
#' 
#' @param spec_quant Data.frame,
#' Annotation file 
#' (output from Spectronaut).
#' 
#' @return A data.frame
#' with the necessary columns for the FragPipe
#' LFQ protein-level analysis

create_annotation <- function(spectro_cond){
  #############################
  # CREATE ANNOTATION FILE
  #############################
  #Select only useful columns
  colnames(spectro_cond) <- tolower(colnames(spectro_cond))
  spectro_cond <- spectro_cond[c('run.label', 'condition')]
  
  #Order by condition 
  spectro_cond <- spectro_cond[order(spectro_cond$condition),]
  row.names(spectro_cond) <- NULL
  
  #Add 'replicate' column
  spectro_cond$replicate <- rep(NA,nrow(spectro_cond))
  
  #Obtain list of conditions
  conditions <-  unique(spectro_cond$condition)
  
  #Add replicate numbers to the condition groups
  for (cond in conditions){
    idx <- which(spectro_cond$condition==cond)
    n_reps <- c(1:length(spectro_cond$replicate[idx]))
    spectro_cond$replicate[idx] <- n_reps
  }
  
  #Rename columns to FragPipe annotation file
  colnames(spectro_cond) <- c("sample", "condition", "replicate")
  
  #Remove duplicated samples
  spectro_cond <-  subset(spectro_cond, 
                              !(rownames(spectro_cond) %in% rownames(spectro_cond)[duplicated(spectro_cond$sample)]))
  
  spectro_cond['sample_name'] = spectro_cond['sample']
  return(spectro_cond)
}

create_quant <- function(spec_quant){
  #############################
  # TRANSFORM PROTEIN FILE
  #############################
  
  #Remove contaminants
  remove_contam_spectro <- function(df){
    return(df[-grep("CONT|iRT", df$`Protein ID`),]) 
  }
  
  #Dictionary of Spectronaut to FragPipe equivalences (not samples)
  spec_to_frag <- c("PG.ProteinAccessions" = "Protein ID",
                    "PG.Genes" = "Gene",
                    "PG.ProteinDescriptions" = "Description",
                    "PG.ProteinNames" = "Entry Name")
  
  #Rename the columns
  for (cols in names(spec_to_frag)){
    colnames(spec_quant)[which(colnames(spec_quant)==cols)] <-  spec_to_frag[cols]
  }
  
  
  # Remove QValue or any other column we don't need
  spec_quant <- spec_quant[c(spec_to_frag, colnames(spec_quant)[grep("Quantity",colnames(spec_quant))])] 
  spec_cols <- colnames(spec_quant)
  
  # Add Protein and Combined Total Peptides columns
  spec_quant["Protein"] <-  paste("sp",spec_quant$`Protein ID`,spec_quant$`Entry Name`,sep = '|')
  spec_quant["Combined Total Peptides"] <- rep(0,nrow(spec_quant))
  
  # Obtain sample names
  spec_sample_names <- spec_cols[grep("Quantity",colnames(spec_quant))]
  
  # Rename the quantity columns to sample intensity 
  spec_sample_names_short <-  c()
  for (spec_sample in spec_sample_names){
    sample_short <-  substr(spec_sample, unlist(gregexpr("\\.\\.", spec_sample))+2, 
                          unlist(gregexpr(".PG.Quantity", spec_sample))-1)
    spec_sample_names_short <-  c(spec_sample_names_short, sample_short)
    colnames(spec_quant)[which(colnames(spec_quant)==spec_sample)] <-paste(sample_short, "Intensity", sep=" ")  
  }
  
  # Remove rows that contain more than one protein NEED TO CORRECT THIS
  spec_quant <- spec_quant[!grepl(";", spec_quant$`Protein ID`), ]
  
  #Remove duplicated samples
  spectro_quant <-subset(spec_quant, select=!(duplicated(colnames(spec_quant))))

  return(spectro_quant)
}





