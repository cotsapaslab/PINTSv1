#' Load protein protein interaction network and annotate each node by Ensemble gene symbol (default). 
#' We provide InWeb3 protein interaction network (Pajek graph format) and the annotation file (text table). 
#' They are located under /Data. User can provide their own graph and the annotation file.   
#'
#' @title Load protein protein interaction network 
#' @param networkfile a pajek format graph file  
#' @param nodeAnnotfile A file that contains the network node symbols and annotations. i.e, Ensemble/HUGO gene symbol 
#' @return a graph object with node names annotated by "Ensemble" gene symbol 
#' @return a data frame with node annotation (Ensemble and HUGO gene symbol etc.)   
#' @export 
#' @import igraph

loadPPINetwork <- function(networkfile, nodeAnnotfile){
  
  #Check file formats and show warning
  nodeAnnot <- read.table(nodeAnnotfile, header = TRUE, stringsAsFactors=FALSE)
  g <- read.graph(networkfile, format="pajek")
  g <- simplify(g)
  #Use Ensemble as a cannonical gene symbol
  if ("Ensemble" %in% colnames(nodeAnnot)) V(g)$name <- nodeAnnot[,"Ensemble"] 
  
  #Retures network in "pajek' and 'graphNEL' format
  cat("PPI Network loaded with igraph pajek format with nodeAnnot:\n 
      g : igraph (pajek) \n
      nodeAnnot : geneAnnotation (Ensemble, HUGO, Genomic position etc.) \n  ")
  
  list(g = g, nodeAnnot = nodeAnnot)  
}
