#' test the significance of the proportion
##' @title test the significance of the proportion
##' @param diseaseNetwork a phenotype-specific PPI network (Pajek format graph)
##' @param geneExpression a data frame for the tissue specificity
##' @param targetCells  A set of cells/tissues under test
##' @param topSubnet top subnetwork genes
##' @param TSSigThres the significance cutoff for the tissue specificity
##' @return p values of Fisher's exact test
##' @export

testSigProportion <- function(diseaseNetwork, geneExpression, targetCells, topSubnet, TSSigThres){

  allCells <- colnames(geneExpression)[-c(1)]
  allGenes <- geneExpression[,"Ensemble"]

  topsubGE <- geneExpression[which(geneExpression[,"Ensemble"] %in% topSubnet), ]
  topsubGenes <- topsubGE[,"Ensemble"]

  contingencyTables <- list()
  pvalueList <- list()
  for(cell in targetCells){
    #The group1 : All expressed genes
    allTSGenes <- geneExpression[which(geneExpression[ ,cell] > TSSigThres) ,"Ensemble"]

    #The group2 : The topsubnet genes
    TSGenes <- topsubGE[which(topsubGE[ ,cell] > TSSigThres) ,"Ensemble"]

    #Category 1 : TS genes
    a <- length(TSGenes) #TS genes in top subnetwork
    b <- length(allTSGenes) - length(TSGenes) # TS genes outside the top subnetwork

    #Category 2: Non TS genes
    allNonTSGenes <- length(allGenes) - length(allTSGenes)
    c <- length(topsubGenes) - a #Non-TS genes in the top subnetwork
    d <- allNonTSGenes - c       #Non-TS genes outside the top subnetwork 
    
    #Fisher's exact test on the over-representaion
    table <- rbind(c(a,b), c(c,d))
    colnames(table) <- c("Topsub","Non-Topsub")
    rownames(table) <- c("TS","NonTS")

    #print(table)
    result <- fisher.test((table) , alternative="greater")

    contingencyTables[[cell]] <- table
    pvalueList[[cell]] <- result$p.value
    #cat(cell, result$p.value, "\n")
  }

  list(contingencyTables = contingencyTables, pvalues = pvalueList)
}
