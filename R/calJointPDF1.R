#' calculate joint probability distribution (PDF) of a UGM model
##' @title calculate joint PDFs of a UGM model
##' @param diseaseNetwork a phenotype-specific PPI network (Pajek format graph)
##' @param topSubnet top subnetworks
##' @param edgePotList the edge potentials for the top subnetworks
##' @param nodeStatesList the on/off states of nodes fro the top subnetworks
##' @param logZvals the logarithmic value of the partition function
##' @param nNull the number of null top subnetworks
##' @return the joint PDFs of the top subnetworks
##' @import igraph

calJointPDF1 <- function(diseaseNetwork, topSubnet, edgePotList, nodeStatesList, logZvals, nNull){

  All.subnets <- topSubnet[1:nNull]
  All.cells <- colnames(nodeStatesList[[1]])[-c(1)]

  factorProdM <- data.frame( matrix(nrow = nNull, ncol = length(All.cells)) )
  jointPDF <- data.frame( matrix(nrow = nNull, ncol = length(All.cells)) )
  colnames(factorProdM) <- All.cells
  colnames(jointPDF) <- All.cells

  k <- 0
  for(subnet in All.subnets){
    k <- k + 1
    g2 <- induced.subgraph(diseaseNetwork, subnet)
    g2.adj <- get.adjacency(g2)

    edgePot <- edgePotList[[k]]
    nodeStates <- as.data.frame(nodeStatesList[[k]], stringsAsFactors=F)

    for(cell in All.cells) {
      edgeStates <- calEdgeStates(g2, g2.adj, edgePot, nodeStates, cell)
      factorProdM[k,cell] <- sum(log(edgeStates))
      jointPDF[k,cell] <- exp( factorProdM[k,cell] - logZvals[k] )
    }
  }
  jointPDFSig <-  apply(jointPDF, 2, function(x)  trunc(rank(x))/length(x) )

  cat("It returns a list: the joint PDF and the significance against null top subnetworks. \n")

  list(jointPDF, jointPDFSig)
}
