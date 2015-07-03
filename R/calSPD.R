#' get Shortest path distance among a set of nodes
##' @title get Shortest path distance among a set of nodes
##' @param graph a network/graph
##' @param nodelist a set of nodes
##' @return shortest path distance
##' @return all the nodes in shortest path distance
##' @export
##' @import igraph

calSPD <- function(graph, nodelist){

  SPDs <- c()
  AllSPDs <- c()
  nodelist <- nodelist
  for(i in 1:(length(nodelist)-1)){
    for(j in (i+1):length(nodelist)){

      node1.id <- which(V(graph)$name %in% nodelist[i])
      node2.id <- which(V(graph)$name %in% nodelist[j])
      SPD <- get.shortest.paths(graph, node1.id , node2.id)
      SPD <- SPD$vpath[[1]]
      SPD <- length(SPD) - 1
      AllSPD <- get.all.shortest.paths(graph, node1.id , node2.id)
      for (AllSPD in AllSPD$res){
        AllSPDs <- c(AllSPDs, AllSPD)
       }
      SPDs <- c(SPDs, SPD)
    }
  }

  AllPathNodes <- V(graph)$name[unique(AllSPDs)]

  list(SPDs = SPDs,  AllPathNodes = AllPathNodes)
}
