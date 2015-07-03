#' test rank enrichment of gene sets
##' @title test rank enrichment of gene sets
##' @param diseaseNetwork a phenotype-specific PPI network (Pajek format graph)
##' @param geneExpression a data frame for the tissue specificity
##' @param targetCells A set of cells/tissues under test
##' @param topSubnet top subnetwork genes
##' @param TSSigThres the significance cutoff for the tissue specificity
##' @return a list of p values: KS standard test and KS bootstrap test
##' @export

testRankEnrichment <- function(diseaseNetwork, geneExpression, targetCells, topSubnet, TSSigThres){

  allCells <- colnames(geneExpression)[-c(1)]
  allGenes <- geneExpression[,"Ensemble"]

  topsubGE <- geneExpression[which(geneExpression[,"Ensemble"] %in% topSubnet), ]
  topsubGenes <- topsubGE[,"Ensemble"]
  KStest1 <- as.data.frame(matrix(nrow=4, ncol = length(targetCells)))
  KStest2 <- as.data.frame(matrix(nrow=4, ncol = length(targetCells)))
  colnames(KStest1) <- targetCells
  colnames(KStest2) <- targetCells
  rownames(KStest1) <- c("All","TS","NonTS", "TSConn")
  rownames(KStest2) <- c("All","TS","NonTS", "TSConn")
  
  #Induce the top subnetwork
  sigSubnet <- induced.subgraph(diseaseNetwork, topsubGenes)
  
  for (cell in targetCells){
    TSGenes <- topsubGE[which(topsubGE[ ,cell] > TSSigThres) ,"Ensemble"]
    NonTSGenes <- setdiff(topsubGenes, TSGenes)
    SPDs <- calSPD(sigSubnet, TSGenes)
    TSConnGenes <- setdiff(SPDs$AllPathNodes, TSGenes)

    #print(cell)
    #cat(length(topsubGenes), length(TSGenes), length(NonTSGenes), "\n")

    TSGeneScores <- topsubGE[topsubGE[,"Ensemble"] %in% TSGenes, cell]
    NonTSGeneScores <- topsubGE[topsubGE[,"Ensemble"] %in% NonTSGenes, cell]
    allGeneScores <- topsubGE[, cell]
    TSConnGeneScores <- topsubGE[topsubGE[,"Ensemble"] %in% TSConnGenes, cell]

    TSRefScores <- geneExpression[which(!(geneExpression[,"Ensemble"] %in% TSGenes)), cell]
    NonTSRefScores <- geneExpression[which(!(geneExpression[,"Ensemble"] %in% NonTSGenes)), cell]
    allRefGeneScores <- geneExpression[which(!(geneExpression[,"Ensemble"] %in% topSubnet)), cell]
    TSConnRefScores <- geneExpression[which(!(geneExpression[,"Ensemble"] %in% TSConnGenes)), cell]

    x <- allGeneScores
    x1 <- TSGeneScores
    x2 <- NonTSGeneScores
    x3 <- TSConnGeneScores
   # print("TS Connector gene scores \n")
   # print(sort(x3))
    boxplot(x, x1, x2, x3, names=c("All", "TS", "Non-TS", "TSConn"), main=cell)
   # cat(length(x), length(x1), length(x2), length(x3), "\n")

    y <- allRefGeneScores
    y1 <- TSRefScores
    y2 <- NonTSRefScores
    y3 <- TSConnRefScores
    boxplot(y, y1, y2, y3, names=c("All", "TS", "Non-TS", "TSConn"), main=paste(cell,"_Ref",sep="") )
  #  cat(length(y), length(y1), length(y2), length(y3),  "\n")

    #standard KS test
    KStest1[1, cell] <- calKStest(x, y, "standard")
    KStest1[2, cell] <- calKStest(x1, y1, "standard")
    KStest1[3, cell] <- calKStest(x2, y2, "standard")
    KStest1[4, cell] <- calKStest(x3, y3, "standard")

    #bootstrap version of KS test
    KStest2[1, cell] <- calKStest(x, y, "bootstrap")
    KStest2[2, cell] <- calKStest(x1, y1, "bootstrap")
    KStest2[3, cell] <- calKStest(x2, y2, "bootstrap")
    KStest2[4, cell] <- calKStest(x3, y3, "bootstrap")
  }

  list(KSstandard = KStest1, KSbootstrap = KStest2)
}


