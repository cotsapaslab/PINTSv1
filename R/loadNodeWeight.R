#' Load node weights (e.g. per gene association p values) for a protein protein interaction network. 
#' The input file is required to have "Ensemble" or "HUGO" field as a gene symbol and "pvalues" field as per-gene association.    
#' We provide per gene scores from mutationally constraint study as an example. It is located under /Data.   
#'
#'@title Load node weights (association p values) 
#'@param nodeWeight Node weight file
#'@param nodeAnnot Node annotation (data frame)
#'@param nodeNameby Gene Symbol option. It takes either "Ensemble" or "HUGO".  
#'@return nodeWeigtTable (a data frame sorted according to InWeb node annotation.)  
#'@export

loadNodeWeight <- function(nodeWeight, nodeAnnot, nodeNameby){
  
  if (nodeNameby == "HUGO"){
    cat("HUGO names are to Ensemble name as InWeb PPI network use Ensemble symbol as a cannonical name. \n ")
    nodeWeight <- mergeWithOrder(nodeAnnot, nodeWeight, by = nodeNameby, keep_order = 1)
    nodeWeight <- nodeWeight[,c("Ensemble", "pvalues")]
    #nodeWeight <- subset(nodeWeight, select=-c(HUGO))
    
  }
  if (nodeNameby == "Ensemble"){ 
    cat("Ensemble names are sorted according to InWeb gene list. \n ")
    nodeWeight <- nodeWeight[match(nodeAnnot[,"Ensemble"], nodeWeight[,"Ensemble"]), ] 
  }
  
  #Needs to consider the situation where there are missing genes. 
  
  return(nodeWeight)  
}