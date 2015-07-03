#' show subNet search result on the fly
##' @title show subNet search result on the fly x
##' @param diseaseNetwork a phenotype-specific PPI network (Pajek format graph)
##' @param nodeAnnot node annotation (data frame)
##' @param module the detected subnetwork
##' @param scores the node scores
##' @return summary
##' @import BioNet

showSubNetSearchResult <- function(diseaseNetwork, nodeAnnot, module, scores){

  #Print out basic information of the retrieved subnetwork
  cat("The subnet nodes", nodes(module), "\n")
  cat("The score of nodes", scores[nodes(module)], "\n")
  cat("The score sum of all nodes", sum(scores[nodes(module)]), "\n")
  cat("The score sum of signal nodes", sum(scores[nodes(module)] > 0), "\n")
  cat("The score sum of null nodes", sum(scores[nodes(module)] < 0), "\n")

  #Need further work.
  #Plot the subnetwork
  #subnetNodes <- attributes(module)[[1]]
  #subnetNodesHUGO <- nodeAnnot[ which(nodeAnnot[,"Ensemble"] %in% subnetNodes), "HUGO"]
  #sigSubnet <- induced.subgraph(diseaseNetwork, subnetNodes)
  #plot(sigSubnet, layout = layout.fruchterman.reingold, main = "signal", vertex.color = "skyblue",
  #     vertex.size = 10, vertex.label = subnetNodesHUGO, vertex.label.font = 3, vertex.label.cex = 0.5,
  #     edge.width = 2)

}
