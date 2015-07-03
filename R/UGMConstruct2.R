#' UGM construct with null top subnetwork
##' @title UGM construct with null top subnetwork
##' @param diseaseNetwork A phenotype-specific PPI network
##' @param geneExpression a data frame for the tissue specificity
##' @param topSubnets a list of detected top subnetworks
##' @param nNull the number of null top subnetworks
##' @return a list of UGM construct outputs
##' @import igraph

UGMConstruct2 <- function(diseaseNetwork, geneExpression, topSubnets, nNull){

  subnet <- topSubnets
  gE <- geneExpression
  All.cells <- colnames(gE)[-c(1)]

  adjMatList <- list()
  nodePotList <- list()
  edgePotList <- list()
  nodeStatesMList <- list()

  #subnet geneExpression (tissue specific)
  subnet.gE <- gE[gE[,"Ensemble"] %in% subnet, ]
  colnames(subnet.gE)[1] <- c("gene")

  #adjacency matrix, nodePost, nodeStates, edgePotentials of a subnet
  sigSubnet <- induced.subgraph(diseaseNetwork, subnet)
  sigSubnet.adj <- get.adjacency(sigSubnet, sparse=FALSE)
  node.pot <- as.data.frame(matrix(1, ncol = 2, nrow = length(V(sigSubnet)$name) ))
  nodeStates <- subnet.gE[ match(V(sigSubnet)$name, subnet.gE$gene), ]
  All.nodes <- nodeStates$gene

  for(k in 1:nNull){

    edge.pot <- data.frame()
    if (k==1){
      colnames(nodeStates) <- c("gene", All.cells)
      edge.pot <- calEdgePot(sigSubnet, sigSubnet.adj, nodeStates, All.cells)
      nodeStatesM <- nodeStates
    }else{
      nodeStatesM <- as.matrix(nodeStatesM[,-c(1)])
      nodeStatesM <-  t(apply(nodeStatesM, 1, function(x) sample(x) ))
      nodeStatesM <- as.data.frame(cbind(All.nodes, nodeStatesM), stringsAsFactors=F)
      colnames(nodeStatesM) <- c("gene", All.cells)
      edge.pot <- calEdgePot(sigSubnet, sigSubnet.adj, nodeStatesM, All.cells)
    }

    adjMatList[[k]] <- sigSubnet.adj
    nodePotList[[k]] <- node.pot
    edgePotList[[k]] <- edge.pot
    nodeStatesMList[[k]] <-  nodeStatesM
  }

  list(adjMatList = adjMatList, nodePotList = nodePotList, edgePotList = edgePotList, nodeStatesList = nodeStatesMList)
}
