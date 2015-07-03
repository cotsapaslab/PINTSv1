#' calculate the edge potential of UGM
##' @title calculate the edge potential
##' @param sigSubnet the top subnetwork graph 
##' @param sigSubnet.adj the adjcency matrix of the top subnetwork
##' @param nodeStates the tissue specificity matrix for the top subnetwork genes
##' @param All.cells all the cells in the tissue specificity matrix
##' @return the edge potentials of a UGM (data frame)
##' @import igraph

calEdgePot <- function(sigSubnet, sigSubnet.adj, nodeStates, All.cells){

  edge.pot <- data.frame()
  epsilon <- 0.000001
  cnt <- 0

  for(i in 1:dim(sigSubnet.adj)[1]){
    for(j in i:dim(sigSubnet.adj)[2]){
      if(sigSubnet.adj[i,j] == 1){
        cnt <- cnt + 1
        n1 <- V(sigSubnet)$name[i]
        n2 <- V(sigSubnet)$name[j]

        if ( n1 %in% nodeStates[,"gene"]  && n2 %in% nodeStates[,"gene"] ){
          #n1.vec <- as.numeric( subset(nodeStates, gene == n1, select = All.cells ) )
          #n2.vec <- as.numeric( subset(nodeStates, gene == n2, select = All.cells ) )
          n1.vec <- as.numeric( nodeStates[which(nodeStates$gene == n1), All.cells] )
          n2.vec <- as.numeric( nodeStates[which(nodeStates$gene == n2), All.cells] )
          n1.on <- which(n1.vec == 1)
          n1.off <- which(n1.vec == 0)
          n2.on <- which(n2.vec == 1)
          n2.off <- which(n2.vec == 0)
          psi11 <- length( intersect(n1.on, n2.on) )
          psi10 <- length( intersect(n1.on, n2.off) )
          psi01 <- length( intersect(n1.off, n2.on) )
          psi00 <- length( intersect(n1.off, n2.off) )
          edge.pot[cnt,1] <- (psi00 + epsilon)
          edge.pot[cnt,2] <- (psi01 + epsilon)
          edge.pot[cnt,3] <- (psi10 + epsilon)
          edge.pot[cnt,4] <- (psi11 + epsilon)
          #cat(edge.pot[i,1], "\t", edge.pot[i,2], "\t", edge.pot[i,3], "\t", edge.pot[i,4],"\t", sum(edge.pot[i,1:4]), "\n")
        }
        else{
          stop("Edge potential is not defined! \n")
          #edge.pot[i,4] <- length(All.cells)/4
          #edge.pot[i,3] <- length(All.cells)/4
          #edge.pot[i,2] <- length(All.cells)/4
          #edge.pot[i,1] <- length(All.cells)/4
        }
      }
    }
  }
  return(edge.pot)
}
