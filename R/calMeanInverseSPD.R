#' calculate mean Inverse Shortest path distance among activeNodes
##' @title calculate mean inverse shortest path distance
##' @param sigSubnet the top subnetwork
##' @param nodeStates the on/off states of nodes
##' @param cell a cell/tissue
##' @return mean inverse shortest path distance
##' @export
##' @import igraph 

calMeanInverseSPD <- function(sigSubnet, nodeStates, cell){

  nodes.On <- nodeStates[which(nodeStates[,cell] == 1), "gene"]
  if(length(nodes.On) > 1){
      SPDs <- calSPD(sigSubnet, nodes.On)
      SPDs <- SPDs$SPDs
      #AllPathNodes <- SPDs$AllPathNodes
      #print(SPDs)
      MISPD <- sum(1/SPDs)/length(SPDs)
  }else{
      MISPD <- 0
  }

  return(MISPD)
}
