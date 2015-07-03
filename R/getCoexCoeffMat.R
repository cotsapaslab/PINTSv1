#' get coexpression coefficient matrix 
##' @title get coexpression coefficient matrix 
##' @param entrezID a set of Entrez ID
##' @param coexDB Coexpression DB 
##' @return a coexpression coefficient matrix

getCoexCoeffMat <- function(entrezID, coexDB){
  
  msize <- length(entrezID)
  coexpress_matrix <- matrix(nrow = msize, ncol = msize)
  for (i in 1:(msize-1) ) {
    id1 <- entrezID[i]
    for (j in (i+1):msize){
      id2 <- entrezID[j]
      coexpress_matrix[i,j] <- getCoexCoeff(id1, id2, coexDB)
      coexpress_matrix[j,i] <- coexpress_matrix[i,j]
    }
  }
  diag(coexpress_matrix) <- 1
  return(coexpress_matrix)
}
