#' simulate gene expression data
##' @title simulated gene expression data
##' @param nodeAnnot node annotation data frame
##' @param allCells a set of cells/tissues
##' @param topSubnet the top subnetwork
##' @param targetCells a set of target cells/tissues
##' @param geneNameBy gene symbol options, e.g., "Ensemble" or "HUGO"
##' @return simulated gene expression
##' @export

simGeneExpression <- function(nodeAnnot, allCells, topSubnet, targetCells, geneNameBy="Ensemble"){

  All.Genes <- nodeAnnot[, geneNameBy]
  nCells <- length(allCells)
  nGenes <- length(All.Genes)

  #simulate null tissue specific gene expression
  MAT <- matrix( rnorm(nGenes*nCells, mean=0, sd=1), nGenes, nCells)
  MAT <- ifelse( MAT > 2, 1, 0)
  gExpression <- data.frame( All.Genes, MAT, stringsAsFactors=F)
  colnames(gExpression) <- c("Ensemble",allCells)

  #Set the signal gene specity to 1
  gExpression[ which(All.Genes %in% topSubnet) , targetCells] <- 1
  #print(str(gExpression))
  return(gExpression)
}
