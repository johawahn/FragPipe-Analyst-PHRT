

transform <-  function(df, ano){
  # This is a function that processes a FragPipe protein
  # file into an intensity matrix for use in RFE++
  # Parameters:
  #   - df: FragPipe protein report
  #   - ano: FragPipe annotation file
  # Returns:
  #   Intensity matrix with 'class' column
  load_frag_matrix <- function(X){
    # This is a function that converts the df from a wide
    # format to a long format to obtain intensities of each 
    # protein per sample
    # Parameters:
    #   - X: protein report
    # Returns:
    #   Long formatted df
    X <- X[c("Protein.ID", colnames(X)[7:ncol(X)])] 
    X_long <-  gather(X, sample, intensity, 
                      colnames(X)[2]:colnames(X)[ncol(X)], 
                      factor_key=T)
    
    return(X_long)
  }
  intensity_matrix <-  function(X){
    # This is a function that converts the df from a long
    # format to a wide format again to have samples as rows,
    # and proteins as columns. It also log normalizes intensity
    # Parameters:
    #   - X: protein report
    # Returns:
    #   Wide formatted intensity df
    X["intensity"] <- log10(X["intensity"])
    
    X <- pivot_wider(X, id_cols = sample, 
                     names_from = Protein.ID, 
                     values_from = intensity)
    X$sample <- sub(".Intensity", "", X$sample)

    return(X)
  }
  add_class <-  function(X,ano){
    # This is a function adds class column necessary for
    #RFE++ classification
    # Parameters:
    #   - X: Intensity matrix
    #   - ano: annotation file
    # Returns:
    #   Intensity matrix with class column
    ano_df <- data.frame(keys = ano$sample,
                         values = ano$condition)
    
    ano_dict <-  setNames(as.character(ano_df$values), ano_df$keys)
    X$class <- unname(ano_dict[X$sample])
    
    return(X)
  }
  
  df <- load_frag_matrix(df)
  df <- intensity_matrix(df)
  df <- add_class(df, ano)

  return(df)
}

                           





































