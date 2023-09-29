#' FragPipe LFQ peptide report to TMT-peptide level report
#'
#' \code{quant_lfq_to_tmt} creates a TMT-type peptide-level file
#'
#' @param lfq_df_path Data.frame,
#' Peptide level FragPipe report
#' @param lfq_type string,
#' Intensity type
#' 
#' \code{anot_lfq_to_tmt} creates a TMT-type annotation file
#' 
#' @param ano Data.frame,
#' Annotation file 
#' (output from FragPipe).
#' 
#' \code{quant_spectro_to_tmt} creates a TMT-type peptide-level file
#' 
#' @param spectro_df Data.Frame,
#' Peptide level Spectronaut report
#' 
#' \code{anot_spectro_to_tmt} creates a TMT-type annotation file
#' 
#' @param ano Data.frame,
#' Annotation file 
#' (output from Spectronaut).
#' 
#' @return A data.frame
#' with the necessary columns for the TMT peptide
#' level analysis

quant_lfq_to_tmt <- function(df, lfq_type){
  if (lfq_type == "Intensity") {
    lfq_columns <- setdiff(grep("Intensity", colnames(df)),
                           grep("MaxLFQ", colnames(df)))
    lfq_columns <- setdiff(lfq_columns, grep("Total Intensity", colnames(df)))
    lfq_columns <- setdiff(lfq_columns, grep("Unique Intensity", colnames(df)))
  } else if (lfq_type == "MaxLFQ") {
    lfq_columns<-grep("MaxLFQ", colnames(df))
    if (length(lfq_columns) == 0) {
      stop(safeError("No MaxLFQ column available. Please make sure your files have MaxLFQ intensity columns."))
    }
  } else if (lfq_type == "Spectral Count") {
    lfq_columns <- grep("Spectral", colnames(df))
    lfq_columns <- setdiff(lfq_columns, grep("Total Spectral Count", colnames(df)))
    lfq_columns <- setdiff(lfq_columns, grep("Unique Spectral Count", colnames(df)))
  }
  
  samples_og <- colnames(df)[lfq_columns]
  interesting_cols <- c("Peptide Sequence", "Modified Sequence", "Charges", 
                        "Protein ID", "Gene", samples_og)
  
  df_sub <-  subset(df, select=interesting_cols)
  
  #Add index and reference intensity
  df_sub <- cbind("Index"=paste(df$`Protein ID`, df$`Peptide Sequence`, sep="_"), 
                  "ReferenceIntensity"=as.numeric(rowMeans(df_sub[colnames(df_sub) %in% samples_og])),
                  "MaxPepProb"=rep(1,nrow(df)), 
                  df_sub) %>%
    rename("Peptide"="Peptide Sequence",
           "ModPeptides"="Modified Sequence",
           "ProteinID"="Protein ID")
  
  
  #Rename columns
  
  samples <- unlist(strsplit(samples_og, " Intensity"))
  colnames(df_sub)[colnames(df_sub) %in% samples_og]  <- samples
  
  df_sub[is.na(df_sub)] <- 0
  
  return(df_sub[c("Index", "Gene", "ProteinID",	
                  "Peptide", "MaxPepProb", "ReferenceIntensity",
                  "ModPeptides", samples)])
}

anot_lfq_to_tmt <- function(df){
  
  
  df <- df[!colnames(df) %in% c("file")]
  
  df <- cbind("Index"=c(1:nrow(df)), 
              "plex"=rep("TMT_01s", nrow(df)), 
              "channel"=rep(0,nrow(df)), 
              df[c("sample", "sample_name", "replicate", "condition")])
  return(df)
}

quant_spectro_to_tmt <- function(spectro_df){
  info_proteins <- read.csv("local_database/uniprotkb_proteome_UP000005640_2023_09_06.tsv", 
                            sep='\t') %>%
    column_to_rownames('Entry')
  
  ptms_df <- read.csv('local_database/230929_spectro_ptms.csv',
                      header=F, row.names = 1)
  
  ptms <- unlist(as.vector(ptms_df))
  names(ptms) <- row.names(ptms_df)
  
  spectro_df <- spectro_df[-grep("CONT", spectro_df$PG.ProteinAccessions),] # Remove contaminants
  
  process <- spectro_df %>%
    select(PG.ProteinAccessions, EG.PrecursorId, 
           colnames(.)[grep("EG.TotalQuantity",colnames(.))]) %>%
    rename('Protein ID'=PG.ProteinAccessions, 'Modified Sequence'=EG.PrecursorId) %>%
    mutate('Modified Sequence'=unlist(lapply(spectro_df$EG.PrecursorId, function(x){unlist(strsplit(x,"\\_"))[2]})),
           .after='Protein ID') %>%
    mutate('Peptide Sequence'=unlist(lapply(`Modified Sequence`, function(x){str_replace_all(x, "\\[[^\\]]*\\]", "")})),
           .after='Modified Sequence') %>%
    mutate('Charges'=unlist(lapply(spectro_df$EG.PrecursorId, function(x){unlist(strsplit(x,"\\_."))[3]})),
           .after='Peptide Sequence') %>%
    mutate('Gene'=info_proteins[spectro_df$PG.ProteinAccessions,]$Gene.Names..primary.,
           .after='Charges') %>%
    mutate("Modified Sequence" = str_replace_all(`Modified Sequence`, ptms))
  
  sample_idx <- grep("Quantity",colnames(process))
  process[,sample_idx] <- apply(process[,sample_idx], 2, function(x){str_replace_all(x, ',' , '.')})
  process[,sample_idx] <- apply(process[,sample_idx], 2, as.numeric)
  
  colnames(process)[sample_idx] <- paste0(unlist(lapply(colnames(process)[sample_idx], 
                                                        function(x){str_match(x, "\\.{2}(.*?)\\.EG")[,2]})),
                                          " Intensity")
  
  if(sum(duplicated(colnames(process)))!=0){
    stop(safeError("You have duplicated samples in the Peptide report.\n 
                   Please correct this in both annotation and experimental reports"))
  }
  
  process_tmt <- quant_lfq_to_tmt(process, 'Intensity')
  return(process_tmt)
}

anot_spectro_to_tmt <- function(spectro_ano){
  # First transform to FragPipe anotation then TMT
  return(anot_lfq_to_tmt(create_annotation(spectro_ano))) 
}
