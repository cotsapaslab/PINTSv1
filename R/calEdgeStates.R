#' calculate edge states
##' @title calculate edge states
##' @param g2 the top subnetwork
##' @param g2.adj the adjacency matrix of top subnetwork
##' @param edgePot the edge potentails of the top subnetwork edges
##' @param nodeStates on/off states of the top subnetwork nodes
##' @param cell a cell/tissue
##' @return edge states a vector of all the edge potentials
##' @import igraph

calEdgeStates <- function(g2, g2.adj, edgePot, nodeStates, cell){

  edgeStates <- c()
  activeEdge <- 0
  cnt <- 0
  for(i in 1:dim(g2.adj)[1]){
    for(j in i:dim(g2.adj)[2]){
      if( g2.adj[i,j] == 1){
        cnt <- cnt + 1
        n1 <- V(g2)$name[i]
        n2 <- V(g2)$name[j]
        n1.state <- as.numeric(nodeStates[nodeStates[,"gene"] == n1, cell])
        n2.state <- as.numeric(nodeStates[nodeStates[,"gene"] == n2, cell])

        n1.id <- which(V(g2)$name %in% n1)
        n2.id <- which(V(g2)$name %in% n2)

        if( n1.state == 1 && n2.state == 1){ state <- 4 ; }
        if( n1.state == 1 && n2.state == 0){ state <- 3 ; }
        if( n1.state == 0 && n2.state == 1){ state <- 2 ; }
        if( n1.state == 0 && n2.state == 0){ state <- 1 ; }

        edgeStates <- c(edgeStates, edgePot[cnt,state])
      }
    }
  }
  return(edgeStates)
}
