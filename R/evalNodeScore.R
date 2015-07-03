#' Evaluete node scores
##' @title Evaluete node scores
##' @param pvalues Per-gene association score for a given phenotype
##' @param sigThres Genome-wide association threshold. e.g., sigThres <- 5e-6
##' @param nIterate It allows iterative disease subnetwork search. Use the default value "nIterate <- 0" for the top subnetwork search.'Please keep the default value until the option is completed.
##' @return calibrated association scores


evalNodeScore <- function(pvalues, sigThres, nIterate=0){
  if (nIterate == 0){
    scores <- -log10(as.numeric(pvalues)) + log10(as.numeric(sigThres))
  }else{
    cat("need to provide code\n")
    #adj_scores <- find_nodes_to_adjscore(DiseaseNetwork, sigThres, nIterate)
  }
  return(scores)
}
