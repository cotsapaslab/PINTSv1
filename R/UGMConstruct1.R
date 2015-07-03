#' UGM construct with null top subnetwork
##' @title UGM construct with null top subnetwork
##' @param diseaseNetwork A phenotype-specific PPI network
##' @param geneExpression a data frame for the tissue specificity
##' @param topSubnets a list of detected top subnetworks
##' @param nNull the number of null top subnetworks
##' @return a list of UGM construct outputs
##' @import igraph

UGMConstruct1 <- function(diseaseNetwork, geneExpression, topSubnets, nNull){

  All.subnets <- topSubnets # top subnetwork
  gE <- geneExpression # tissue specificity matrix
  All.cells <- colnames(gE)[-c(1)] #All cell types in gene Exprssion data

  subnet.id <- 0
  adjMatList <- list()
  nodePotList <- list()
  edgePotList <- list()
  nodeStatesList <- list()

  for (subnet in All.subnets){

    subnet.id <- subnet.id + 1
    edge.pot <- data.frame()
    node.pot <- data.frame()
    node.states <- data.frame()

    #subnet geneExprssion (tissue specific)
    subnet.gE <- gE[gE[,"Ensemble"] %in% subnet, ]
    colnames(subnet.gE)[1] <- c("gene")

    #adjacency matrix, nodePost, nodeStates, edgePotentials of a subnet
    sigSubnet <- induced.subgraph(diseaseNetwork, subnet)
    sigSubnet.adj <- get.adjacency(sigSubnet, sparse=FALSE)
    node.pot <- as.data.frame(matrix(1, ncol = 2, nrow = length(V(sigSubnet)$name) ))
    node.states <- subnet.gE[ match(V(sigSubnet)$name, subnet.gE$gene),]
    edge.pot <- calEdgePot(sigSubnet, sigSubnet.adj, subnet.gE, All.cells)

    adjMatList[[subnet.id]] <- sigSubnet.adj
    nodePotList[[subnet.id]] <- node.pot
    edgePotList[[subnet.id]] <- edge.pot
    nodeStatesList[[subnet.id]] <-  node.states
  }

  list(adjMatList = adjMatList, nodePotList = nodePotList, edgePotList = edgePotList, nodeStatesList = nodeStatesList)
}
