#' Obtain protein structure and sequence annotation using Protter 
#' \url{wlab.ethz.ch/protter/}
#' 
#' \code{request_protter} Image of the structure and annotation of a requested protein
#' 
#' @param prot_ID string (REQUIRED),
#' UniProt protein accession
#'  
#' @param peptides string,
#' Comma separated list of peptides
#' 
#' @param annotations vector,
#' List of wanted annotations, the ones available are:
#' PTMs, variants, disulfide bonds, signal peptide
#' 
#' @return A magick_image
#' vectorized magick image can be visualized using
#' image_ggplot()
#' 
#' @examples:
#' example <-  request_protter(prot_ID="Q10589")
#' example  <-  request_protter(prot_ID="Q10589", 
#'                  peptides="HLLQQELT,KKYYPSSQ",
#'                  annotations=c('signal_pep','variants','ptms','disulf_bonds'))


request_protter <- function(prot_ID, peptides="", annotations=c()){
  param <- list("up"=prot_ID,
                "tm"="auto",
                "format"="png",
                "legend"="")
  
  if (!nchar(peptides) == 0){
    param["n:Identified Peptides,fc:pink,bc:pink"] <- peptides
  }

  if ('signal_pep' %in% annotations){
    param["n:Signal peptide,fc:red,bc:red"] <- "UP.SIGNAL"
  }

  if ('variants' %in% annotations){
    param["n:Variants,s:diamond,fc:orange,bc:orange"] <- "UP.VARIANT"
  }

  if('ptms' %in% annotations){
    param["n:PTMs,s:box,fc:forestgreen,bc:forestgreen"] <- "UP.CARBOHYD,UP.MOD_RES"
  }

  if('disulf_bonds' %in% annotations){
    param["n:Disulfide bonds,s:box,fc:greenyellow,bc:greenyellow"] <- "UP.DISULFID"
  }
  
  if('lux_mods' %in% annotations){
    param["n:Lux Modifications,s:circ,bc:moccasin"] <- "H,W"
  }
  
  server <- "http://wlab.ethz.ch/protter/create?"
  
  response <- GET(server, query=param)
  image_content <- httr::content(response, "raw")
  magick_image <- image_read(image_content)
  return(magick_image)
}

