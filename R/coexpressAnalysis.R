#' coexpression Analysis 
##' @title coexpression Analysis 
##' @param coexDB Coexpression DB
##' @param coexAnnot Entrez ID annoation
##' @param nodeAnnotPath Node annotation file
##' @param topSubnets Top sub network genes
##' @param nNull The number of null top subnetworks      
##' @return a coexpression coefficient matrix
##' @export

coexpressAnalysis <-function(coexDB, coexAnnot, nodeAnnotPath, topSubnets, nNull){
  
  meanCoexps <- c()
  topSubN <- 1:nNull
  nodeAnnot <- read.table(nodeAnnotPath, header=T, stringsAsFactors=F)
  for(k in topSubN){
    topSub <- topSubnets[[k]]
    Node_index <- nodeAnnot[which(nodeAnnot$Ensemble %in% topSub), "comp0_NID"]
    InWeb_to_Entrez <- convIWIDtoEntrez(Node_index, nodeAnnotPath, coexAnnot, coexDB)
    Coexpression <- getCoexCoeffMat(as.numeric(InWeb_to_Entrez[,3]), coexDB )
    Coexpression <- Coexpression[lower.tri(Coexpression)]
    meanCoexps <- c(meanCoexps, mean(Coexpression))
  }
  return(cbind(topSubN, meanCoexps))
}
