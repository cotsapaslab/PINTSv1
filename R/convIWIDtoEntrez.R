#' convert InWeb ID to Entrez ID
##' @title convert InWeb ID to Entrez ID
##' @param IWID InWeb ID
##' @param nodeAnnotPath Node annotation file
##' @param coexAnnot Entrez ID annoation
##' @param coexDB Coexpression DB 
##' @return InWeb ID to Entrez ID conversion table 

convIWIDtoEntrez <- function(IWID, nodeAnnotPath, coexAnnot, coexDB){
  
  fileA <- nodeAnnotPath
  fileB <- coexAnnot
  dir <- coexDB
  IW_to_Entrez <- matrix(nrow = length(IWID), ncol=3)    
  for(i in 1:length(IWID) ){
    
    IWid <- IWID[i]
    command1 = paste("awk "," \'{ if($1 == ", IWid,") print $3; }\' ",fileA, sep="")
    ENS <- try(system(command1, intern = TRUE))
    command2 = paste("awk "," \'{ if($1 == \"", ENS ,"\" ) print $2; }\' ",fileB, sep="")
    Entrez  <- try(system(command2,intern = TRUE))
    Entrez <- Entrez[1]
    if (file.exists(paste(dir, Entrez, sep=""))){
      Entrez <- Entrez[1]
      IW_to_Entrez[i,] <- cbind(IWid, ENS, Entrez)
    }
    else{
      IW_to_Entrez[i,] <- NA
    }
  }
  return(na.omit(IW_to_Entrez))
}
