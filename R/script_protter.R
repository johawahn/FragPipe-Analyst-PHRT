library(httr)
library(cowplot)

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

#example <- request_protter("A0A0B4J2D5", annotations = c('disulf_bonds'))

