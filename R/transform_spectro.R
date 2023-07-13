#setwd("~/230707_Spectronaut18_trans")
#Sys.setlocale("LC_MESSAGES","en_US.UTF-8")

####################################################
#                 FUNCTIONS
####################################################

create_annotation <- function(spectro_cond){
  #############################
  # CREATE ANNOTATION FILE
  #############################
  
  #Order by condition 
  spectro_cond <- spectro_cond[order(spectro_cond$Condition),]
  row.names(spectro_cond) <- NULL
  
  #Add 'replicate' column
  spectro_cond$replicate <- rep(NA,nrow(spectro_cond))
  
  #Obtain list of conditions
  conditions <-  unique(spectro_cond$Condition)
  
  #Add replicate numbers to the condition groups
  for (cond in conditions){
    idx <- which(spectro_cond$Condition==cond)
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

  return(list(spectro_quant, spec_sample_names_short))
}

spectronaut_to_fragpipe <-  function(spectro_quant, spectro_cond){
  #############################
  # ENSURE THAT QUANTITY AND
  # ANNOTATION FILE HAVE THE
  # SAME SAMPLES
  #############################

  spec_list  <-  create_quant(spectro_quant)
  spec = as.data.frame(spec_list[[1]])
  spec_sample_names_short = unlist(spec_list[[2]])
  
  spec_annotations <- create_annotation(spectro_cond)

  
  correct_samples = intersect(spec_sample_names_short, spec_annotations$sample)
  
  row.names(spec_annotations) <-  spec_annotations$sample
  spec_annotations = spec_annotations[correct_samples,]
  
  spec = spec[,c("Protein ID", "Gene", "Description", "Entry Name", "Protein", 
                "Combined Total Peptides", paste(correct_samples, "Intensity", 
                                                 sep=" "))]
  
  #write.table(spec, file='spec_quantities.tsv', quote=FALSE, sep='\t', 
  #           col.names = TRUE, row.names = FALSE)
  
  
  #write.table(spec_annotations, file='spec_annotations.tsv', quote=FALSE, sep='\t', 
  #            col.names = TRUE, row.names = FALSE)
  
  return(list(spec, spec_annotations))
}


####################################################
#                 MAIN
####################################################

# Read files
#conditions <- read.csv("data/202302_PHRT-15_ES903_ConditionSetup.csv", header = T, sep=",")[c("Run.label","Condition")]
#spec_output <- read.csv("data/20230614_090727_202204_PHRT-15_cohort-1_all_Protein-Report.csv", header = T)


# Transform them
#files <-  spectronaut_to_fragpipe(spec_output, conditions)

#spec_quant <-  as.data.frame(files[[1]])
#spec_annotations <-  as.data.frame(files[[2]])

#write.table(spec_quant, file='spec_quantities_1.tsv', quote=FALSE, sep='\t', 
#           col.names = TRUE, row.names = FALSE)


#write.table(spec_annotations, file='spec_annotations_1.tsv', quote=FALSE, sep='\t', 
#           col.names = TRUE, row.names = FALSE)
