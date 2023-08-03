

quant_lfq_to_tmt <- function(lfq_df_path, lfq_type){
  df <- read.table(lfq_df_path,
                   header = TRUE,
                   fill= TRUE, # to fill any missing data
                   sep = "\t",
                   quote = "",
                   comment.char = "",
                   blank.lines.skip = F,
                   check.names = F)
  
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
                  df_sub)
  

  #Rename columns
  rename <- c("Peptide Sequence"="Peptide",
              "Modified Sequence"="ModPeptides",
              "Protein ID"="ProteinID")
  samples <- unlist(strsplit(samples_og, " Intensity"))
  colnames(df_sub)[colnames(df_sub) %in% names(rename)] <- as.character(rename)
  colnames(df_sub)[colnames(df_sub) %in% samples_og]  <- samples
  
  df_sub[is.na(df_sub)] <- 0
  
  return(df_sub[c("Index", "Gene", "ProteinID",	
                       "Peptide", "MaxPepProb", "ReferenceIntensity",
                       "ModPeptides", samples)])
}

anot_lfq_to_tmt <- function(lfq_ano_path){
  df <- read.table(lfq_ano_path, header = TRUE, sep='\t')
  df <- df[!colnames(df) %in% c("file")]
  
  df <- cbind("Index"=c(1:nrow(df)), 
              "plex"=rep("TMT_01s", nrow(df)), 
              "channel"=rep(0,nrow(df)), 
              df[c("sample", "sample_name", "replicate", "condition")])
  return(df)
}
  








