#' Randomize disease-specific PPI network
#'
##' @title randomize disease-specific PPI network
##' @param diseaseNetwork A a phenotype-specific PPI network
##' @param nPermute number of permutations
##' @param permuteOptions There are two options. "Non-Degree Preserving" - Node degree is not constrained in randomization of nodes in PPI. "Degree Preserving" - Node degree is preserved in randomization of nodes in PPI.
##' @return randomized disease network node assignment table (data frame)
##' @export

randDiseaseNetwork <- function(diseaseNetwork, nPermute, permuteOptions="Non-Degree-preserving"){

  nNodes  <- length(V(diseaseNetwork)$name)
  pvalues <- V(diseaseNetwork)$weight
  if( permuteOptions == "Degree-preserving"){
    stop("Degree preserving random permutation\n")
    #Add code.
  }else if ( permuteOptions == "Non-Degree-preserving") {
          cat("Non Degree preserving random permutation\n")
          names = paste ("null", 1:nPermute, sep= "")
          randDiseaseNetwork <- replicate(nPermute, sample(seq(1:nNodes)) )
          randDiseaseNetwork <- cbind(V(diseaseNetwork)$name, V(diseaseNetwork)$weight, randDiseaseNetwork)
          colnames(randDiseaseNetwork) <- c("Ensemble","pvalues", names)
  }

  return(randDiseaseNetwork)
}
