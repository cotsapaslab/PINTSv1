#' test the significance of the connectivity of a set of expressed genes for target tissues.
##' @title test the Significance of coordination of gene set (e.g., tissue specific genes)
##' @param diseaseNetwork a phenotype-specific PPI network (Pajek format graph)
##' @param geneExpression a data frame for tissue specificity.
##' @param targetCells A set of cells/tissues under test.
##' @param topSubnets Top sub network genes
##' @param nNull number of null gene sets (randomly selected genes in a give tissue)
##' @return p values of MISPD (mean inverse shortest path distance)
##' @export
##' @import igraph

testSigTissueCoherence <- function(diseaseNetwork, geneExpression, targetCells, topSubnets, nNull){

  subnet <- topSubnets
  gE <- geneExpression
  All.cells <- colnames(gE)[-c(1)]

  nodeStatesMList <- list()

  #subnet geneExpression (tissue specific)
  subnet.gE <- gE[gE[,"Ensemble"] %in% subnet, ]
  colnames(subnet.gE)[1] <- c("gene")

  #nodeStates of a subnet
  sigSubnet <- induced.subgraph(diseaseNetwork, subnet)
  nodeStates <- subnet.gE[ match(V(sigSubnet)$name, subnet.gE$gene), ]
  All.nodes <- nodeStates$gene

  for(k in 1:nNull){

    if (k==1){
      colnames(nodeStates) <- c("gene", All.cells)
      nodeStatesM <- nodeStates
    }else{
      nodeStatesM <- as.matrix(nodeStatesM[,-c(1)])
      nodeStatesM <- t(apply(nodeStatesM, 1, function(x) sample(x) ))
      nodeStatesM <- as.data.frame(cbind(All.nodes, nodeStatesM), stringsAsFactors=F)
      colnames(nodeStatesM) <- c("gene", All.cells)
    }
    nodeStatesMList[[k]] <-  nodeStatesM
  }

  mISPD <- data.frame( matrix( ncol=length(targetCells), nrow=length(nodeStatesMList) ) )
  colnames(mISPD) <- targetCells
  nodeState.id <- 0
  for (nodeStates in nodeStatesMList){
    nodeState.id <- nodeState.id + 1
    for (cell in targetCells){
      meanInverseSPD <- calMeanInverseSPD(sigSubnet, nodeStates, cell)
      mISPD[nodeState.id, cell] <- meanInverseSPD
      #cat(cell, meanInverseSPD, "\n")
    }
  }

  #Calculate the significance
  #mISPDSig <-  apply(mISPD, 2, function(x) 1 - trunc(rank(x))/length(x))
  mISPDSig <- t( apply(mISPD, 2, function(x)  1 - ( trunc(rank(x,ties.method = c("average"))) - 1)/nNull ) )
  
  list(nodeStatesList = nodeStatesMList, meanInverseSPD = mISPD, meanInverseSPDSig = mISPDSig)
}
