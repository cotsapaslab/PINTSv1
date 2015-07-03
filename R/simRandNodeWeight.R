#' Simulate per-gene node weights (p values) from signal + null distribution for PPI network nodes.
#' For example, the signal distribution can be sampled from a beta distribution and
#' the null distribution from a uniform distribution.
##' @title Simulate signal/null nodes weights
##' @param nNodes the number of PPI network nodes to simulate
##' @param nSigNodes the number of signal nodes
##' @param sigFun the signal distribution
##' @param nullFun the null distribution
##' @param sigThres the signal threshold
##' @param ... additional parameters to pass
##' @return a numeric vector of p values
##' @export 

simRandNodeWeight <- function(nNodes, nSigNodes, sigFun, nullFun, sigThres, ... ) {

  cat("Simulate the nodes p values from sigFun (Signal) + nullFun (Null) ...\n")
  nTot <- nNodes
  nSig <- nSigNodes
  nNull <- nTot - nSig
  p_values <- numeric(length = nTot)
  N <- 1000 #Number of avearged random Signal/Null distribution

  sig_pval_matrix <- matrix(nrow = N, ncol = nSig)
  null_pval_matrix <- matrix(nrow = N, ncol = nNull)

  for(i in 1:N){
    sig_pval <- sort(sigFun)
    null_pval <- sort(nullFun)
    #sig_pval <- sort(rbeta(nSig, shape1= Beta, shape2=1, ncp = 0))
    #null_pval <- sort(nullFun(nNull, 0, 1))

    sig_pval_matrix[i,] <- sig_pval
    null_pval_matrix[i,] <- null_pval
  }
  avg_sig_pval <- colMeans(sig_pval_matrix)
  avg_null_pval <- colMeans(null_pval_matrix)

  p_values[1:nSig] <- avg_sig_pval
  p_values[(nSig+1): nTot] <- avg_null_pval

  cat("number of p values below signal Threshold (True positive + False positive): ", length(which(p_values < sigThres)), "\n")
  nTP = length(which(p_values[1:nSig] <= sigThres))
  nFN = length(which(p_values[1:nSig] > sigThres))
  nFP = length(which(p_values[(nSig+1): nTot] <= sigThres))
  nTN = length(which(p_values[(nSig+1): nTot] > sigThres))

  cat("number of True positive : ", nTP, "\n")
  cat("number of False positive : ", nFP, "\n")
  cat("number of False negative : ", nFN, "\n")
  cat("number of True negative : ", nTN, "\n")

  return(p_values)
}
