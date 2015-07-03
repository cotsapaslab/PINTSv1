#' get coexpression coefficient 
##' @title get coexpression coefficient 
##' @param id1 a Entrez ID
##' @param id2 a Entrez ID
##' @param coexDB Coexpression DB 
##' @return a coexpression coefficient 

getCoexCoeff <- function(id1, id2, coexDB){
  
  dir <- coexDB
  fileA <- paste(dir, id1,sep="")
  fileB <- paste(dir, id2,sep="")
  if (file.exists(fileA) && file.exists(fileB)){
    command = paste("awk "," \'{ if($1 == ", id2,") print $3; }\' ",fileA, sep="")
    result <- try(system(command,intern = TRUE))
    return(as.numeric(result[1]))
  }
  else{
    return("NA")
  }
  
}
