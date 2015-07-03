#' Provides a random node assignment with a signal subnetwork. The signal subnetwork refers to a
#' connected subnetwork where all the nodes are assumed to be associated with a given phenotype in con#' sideration.
#'
#' @title Simulate random subnetwork from PPI network
##' @param nNodes the number of PPI network nodes
##' @param nSigNodes the number of signal nodes
##' @param nodeAnnot Node annotation (data frame)
##' @param subnetMatrix A sampled signal subnetwork. ex. /Data/randomSubnet_ex1.txt
##' @param nSubnet the number of sub network : Please use a default value for now, ie.,  nSubnet <- 1
##' @param nPermute number of permutations with the Subnet : Please use a default value for now, ie.,  nPermute <- 1
##' @return node_table - A vector (matrix) of sampled nodes.
##' @return node_table[1:nSigNodes] - signal node index of the network
##' @return node_table[(nSigNodes+1):nNodes] - null node index of the network
##' @export

simRandSubnet <- function(nNodes, nSigNodes, nodeAnnot, subnetMatrix, nSubnet, nPermute){

  cat("Node_sampling_from_network postive control simulation  ......\n")
  nTot <- nNodes
  nSig <- nSigNodes
  nNull <- nTot - nSig
  subnet_size <- nSig + 1
  node_table <- matrix(nrow = nTot, ncol= nSubnet*nPermute)
  #cat("The number of nodes:", nTot, "The number of signal nodes", nSig , "\n")

  for( i in 1:nSubnet ){
    CC <- subnetMatrix[1,i]
    nodelist <- subnetMatrix[2:subnet_size,i]
    #cat("\n",CC ,"\t", nodelist[1:10],"\n")
    #cat("sampled gene index :", sort(data0$name[which(data0$name %in% nodelist)])[1:10], "\n" )
    #cat("sampled gene index :", sort(which(data0$name %in% nodelist))[1:10], "\n" )
    sampled_tot_index <- seq(1:nTot)
    sampled_signal_index <- which(nodeAnnot[,"HUGO"] %in% nodelist)
    sampled_null_index <- sampled_tot_index[which(! (sampled_tot_index %in% sampled_signal_index))]

    col1 <- i + (i-1)*nPermute
    if(nSig > 0) node_table[1:nSig,col1] <- sampled_signal_index
    #cat("origninal:", sort(node_table[1:nSig,col1])[1:10], "\n")
    if(nNull > 0) node_table[(nSig+1):nTot, col1] <- sampled_null_index

    if (nPermute > 1){
      for (j in 1:(nPermute - 1)){
        col2 <- col1 + j
        node_table[1:nSig, col2] <- sample(sampled_signal_index)
        #cat("permuted:", sort(node_table[1:nSig,col2])[1:10], "\n")
        #cat("intersect:", length(intersect(node_table[1:nSig,col1],node_table[1:nSig,col2])),"\n")
        node_table[(nSig+1):nTot, col2] <- sample(sampled_null_index)
      }
    }
  }

  cat("The number of signals is ", length(sampled_signal_index), "\n")
  cat("It returns a table of signal nodes (row: 1 to nSigNodes) and null nodes \n")
  return(node_table)
}
