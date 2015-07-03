#' run disease-associated subnetwork search
##' @title run disease-associated subnetwork search
##' @param diseaseNetwork A phenotype-specific PPI network
##' @param nodeAnnot Node annotation (data frame)
##' @param sigThres Genome-wide association threshold. e.g., sigThres <- 5e-6
##' @param nIterate It allows iterative disease subnetwork search. Use the default value "nIterate <- 0" for the top subnetwork search.
##' @param nullDiseaseNetwork null association disease networks (randomized disease-association network)
##' @param showResult show the top subnetwork search result
##' @return top/null scoring subnetworks SubnetList[[1]] corresponds to top scoring subnet. SubnetList[[2:]] the null top scoring subnets.
##' @export
##' @import igraph
##' @import BioNet


runDiseaseSubnetSearch <- function(diseaseNetwork, nodeAnnot, sigThres, nIterate = 0, nullDiseaseNetwork=NULL, showResult=TRUE){

  cat(" run disease-associated subnetwork search \n ")
  subnetList <- list()
  subnetID <- 0

  print(diseaseNetwork$weight)
  #Evaluate node scorese considering the signal threshold
  pvalues <- V(diseaseNetwork)$weight
  scores <- evalNodeScore(pvalues, sigThres, nIterate=0)
  names(scores) <- V(diseaseNetwork)$name

  #Run the heuristic algorithm
  bioNetDiseaseNetwork <- igraph.to.graphNEL(diseaseNetwork)
  module <- runFastHeinz(bioNetDiseaseNetwork, scores)

  #print the results on the fly
  if (showResult) showSubNetSearchResult(bioNetDiseaseNetwork, nodeAnnot, module, scores)
  subnetNodes <- attributes(module)[[1]]
  subnetID <- subnetID + 1
  subnetList[[subnetID]] <- subnetNodes

  if(!is.null(nullDiseaseNetwork)){

    All.nullNetworks <- colnames(nullDiseaseNetwork)[-c(1,2)]
    for(nullNetwork in All.nullNetworks){

      cat("run ", nullNetwork, "\n")
      randNodeIndex <- as.numeric(nullDiseaseNetwork[,nullNetwork])
      pvalues <- V(diseaseNetwork)$weight[randNodeIndex]
      scores <- evalNodeScore(pvalues, sigThres, nIterate=0)
      names(scores) <- V(diseaseNetwork)$name

      #Run the heuristic algorithm
      bioNetdiseaseNetwork <- igraph.to.graphNEL(diseaseNetwork)
      module <- runFastHeinz(bioNetdiseaseNetwork, scores)

      #print the results on the fly
      if (showResult) showSubNetSearchResult(bioNetdiseaseNetwork, nodeAnnot, module, scores)
      subnetNodes <- attributes(module)[[1]]
      subnetID <- subnetID + 1
      subnetList[[subnetID]] <- subnetNodes
    }
  }

  cat("It returns the list of the detected subnets: \n
      subnetList[[1]] corresponds to top scoring subnet \n
      subnetList[[2:]] the null top scoring subnets. \n ")
  return(subnetList)

}


