#' Load gene expression data
##' @title load gene expression data
##' @param expressionDataFile gene expression data file (tissue specificity)
##' @param nodeAnnot Node annotation (data frame)
##' @param nodeNameby Gene Symbol option. It takes either "Ensemble" or "HUGO"
##' @return gene expression data frame
##' @export 

loadGeneExpression <- function(expressionDataFile, nodeAnnot, nodeNameby="Ensemble"){

  gExpression <- read.table(expressionDataFile, stringsAsFactors = F, header = T)

  #Read gene expressoin data: Work on gene expression data types.
  nodeNameDB <- c("Ensemble", "HUGO")
  if (nodeNameby == "HUGO"){
    cat("HUGO names are converted to Ensemble name as InWeb PPI network use Ensemble symbol as a cannonical name. \n ")
    #gExpression <-  gExpression[ match(nodeAnnot[,"HUGO"], gExpression[,"HUGO"]), ]
    gExpression <- mergeWithOrder(nodeAnnot, gExpression, by = nodeNameby, keep_order = 1)
    nodeWeight <- nodeWeight[,c("Ensemble", "pvalues")]
    #gExpression <- subset(gExpression, select=-c(HUGO))
  }
  if (nodeNameby == "Ensemble"){
    cat("Ensemble names are sorted according to InWeb gene list. \n")
    gExpression <-  gExpression[ match(nodeAnnot[,"Ensemble"], gExpression[,"Ensemble"]), ]
  }
  if (!(nodeNameby %in% nodeNameDB)){
    stop("Plaese check the nodeNames associated nodeWeights are properly entered: \n
        Currently, we take either Ensemble and HUGO symbol.\n
        Make sure that the column names are set either Ensemble or HUGO.\n")
  }

  #Number of genes that matches with InWeb
  cat("The number of genes that match with InWeb gene list is ", dim(gExpression)[1]-1  ,"\n")
  cat("The number of tissues represented is ", dim(gExpression)[2]-1, "\n")

  return(gExpression)

}
