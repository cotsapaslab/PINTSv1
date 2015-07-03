#'Provides a phenotype-specific PPI network: The nodes in the PPI network are assigned weights reflect#'ing the association with a phenotype. The weights are given by -log(p) + log(p_sigThres).
#'
#' @title Build disease-specific PPI network
##' @param network A PPI network (Pajek format graph)
##' @param nodeWeight Node weight (per-gene association p value) and Ensemble gene symbol data frame
##' @param nodeAnnot Node annotation (data frame)
##' @param nodeNameBy Gene Symbol option. It takes either "Ensemble" or "HUGO".
##' @param nodeSelected By default, we use the entire network nodes (nodeSelected = NULL). The set of nodes can be limited to a subset of the entire network.
##' @return A disease-specific PPI network (Pajek format graph)
##' @export
##' @import igraph


buildDiseaseNetwork <- function(network, nodeWeight, nodeAnnot, nodeNameBy, nodeSelected = NULL){

  #Convert and Match the nodeWeights entered with PPI nodeAnnot (Ensemble or HUGO)
  nodeNameDB <- c("Ensemble", "HUGO")
  if (nodeNameBy == "HUGO"){
    cat("HUGO names are converted to Ensemble name. \n
        As InWeb PPI network use Ensemble symbol as a cannonical name. \n ")
    nodeWeight <- nodeWeight[match(nodeAnnot[,"HUGO"], nodeWeight[,nodeNameBy]), ]
  }
  if (nodeNameBy == "Ensemble"){
    cat("Ensemble names are sorted in the order of InWeb gene. \n")
    nodeWeight <- nodeWeight[match(nodeAnnot[,"Ensemble"], nodeWeight[,nodeNameBy]), ]
  }
  if (!(nodeNameBy %in% nodeNameDB)){
    cat("Plaese check the nodeNames associated node weights are properly entered: \n
        Currently, we take either Ensemble and HUGO symbol.\n
        Make sure that the column names are set either Ensemble or HUGO.\n")
  }
  if ("pvalues" %in% colnames(nodeWeight)){
    V(network)$weight <- nodeWeight[,"pvalues"]
  }else{
    cat("Please provide the pvalues.\n")
  }

  #V(network)$weight <- nodeWeight[,"pvalues"]
  #Consider how to deal with missing nodeinfo (e.g., p values)
  #or resticting the network within the expreseed genes across all tissues or a set of tissues.
  #is.missing.info(nodeWeights)
  if (!is.null(nodeSelected)){
    network <- induced.subgraph(network, nodeSelected, "copy_and_delete")
  }

  return(network)
}
