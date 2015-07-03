#' test the significance of the overlap between gene sets from target cell/tissues.
##' @title test the significance of the overlap among target cell/tissues.
##' @param diseaseNetwork A phenotype-specific PPI network (Pajek format graph)
##' @param geneExpression A data frame for the tissue specificity.
##' @param targetCells A set of cells/tissues under test.
##' @param topSubnet Top subnetwork genes
##' @return p values
##' @export

testSigTissueOverlap <- function(diseaseNetwork, geneExpression, targetCells, topSubnet){

#Tissue specificity (binary call) matrix for the top subnetwork.
topsubGE <- geneExpression[ which( geneExpression[,"Ensemble"] %in% topSubnet), ]
allGenes <- topsubGE[,"Ensemble"]
totNumGenes <- length(allGenes)
#print(str(topsubGE))

#Significant overlap between the tissue specific genes in the target tissues
#cat("cell1","cell2","n1", "n2", "nc", "np", "pval", "\n")
TSSigOverlap1 <- as.data.frame(matrix(nrow=0, ncol=7))
TSSigOverlap1 <- as.data.frame(matrix(nrow= length(targetCells)*(length(targetCells)-1), ncol=7 ))
colnames(TSSigOverlap1) <- c("cell1","cell2","n1", "n2", "nc", "np", "pval")
cnt <- 0
for(cell1 in targetCells){
  for (cell2 in targetCells){
    if (cell1 != cell2){
      cnt <- cnt + 1
      TSNodes1 <- topsubGE[which(topsubGE[,cell1] == 1), "Ensemble"]
      TSNodes2 <- topsubGE[which(topsubGE[,cell2] == 1), "Ensemble"]
      comTSNodes <- intersect(TSNodes1, TSNodes2)

      pval <- evalSigIntersect(length(TSNodes1), length(TSNodes2), length(comTSNodes), totNumGenes)
      result <- c(cell1, cell2, length(TSNodes1), length(TSNodes2), length(comTSNodes), totNumGenes, pval)
      TSSigOverlap1[cnt,] <- result
    }
  }
}


#find common connector nodes
#TSconnNodesList <- list()
#cell.id <- 0
#for (cell in targetCells){
#  cell.id <- cell.id + 1
#  sigSubnet <- induced.subgraph(diseaseNetwork, topSubnet)
#  nodeSet <- topsubGE[which(topsubGE[,cell] == 1), "Ensemble"]
#  SPD <- calSPD(sigSubnet, nodeSet)
#  allPathNodes <- SPD$AllPathNodes
#  TSConnNodes <- setdiff(allPathNodes, nodeSet)
  #print(allPathNodes)
  #print(length(allPathNodes))
#  TSconnNodesList[[cell]] <- TSConnNodes
#}

#Significant overlap between the tissue specific genes in the target tissues
#cat("cell1","cell2","n1", "n2", "nc", "np", "pval", "\n")
TSSigOverlap2 <- as.data.frame(matrix(nrow= length(targetCells)*(length(targetCells)-1), ncol=7 ))
colnames(TSSigOverlap2) <- c("cell1","cell2","n1", "n2", "nc", "np", "pval")
#cnt <- 0
#for(cell1 in targetCells){
#  for (cell2 in targetCells){
#    if (cell1 != cell2){
#      cnt <- cnt + 1
#      TSNodes1 <- TSconnNodesList[[cell1]]
#      TSNodes2 <- TSconnNodesList[[cell2]]
#      comTSNodes <- intersect(TSNodes1, TSNodes2)
#      pval <- evalSigIntersect(length(TSNodes1), length(TSNodes2), length(comTSNodes), totNumGenes)
#      result <- c(cell1, cell2, length(TSNodes1), length(TSNodes2), length(comTSNodes), totNumGenes, pval)
#      TSSigOverlap2[cnt,] <- result
#    }
#  }
#}

list(TSOverlap = TSSigOverlap1, TSConOverlap = TSSigOverlap2)

}

