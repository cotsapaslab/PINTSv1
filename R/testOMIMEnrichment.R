#' test OMIM record enrichment of gene set
##' @title test OMIM record enrichment of gene set
##' @param geneOMIM a data frame for gene-OMIM record
##' @param topSubnet top subnetwork genes
##' @return p values of Fisher's exact test and contingency table
##' @export

testOMIMEnrichment <- function(geneOMIM, topSubnet){
  
  #get all gene-OMIM record 
  df1 <- geneOMIM
  df1 <- unique(df1[which(df1$OMIM_ID != "NA"),])
  allOMIMGenes <- length(unique(df1$Ensemble))
  allOMIMIDs <- length(unique(df1$OMIM_ID))
  
  #get gene-OMIM record for the top subnetwork genes
  AllTS <- topSubnet
  df3 <- df1[which(df1$Ensemble %in% AllTS),]
  TSOMIM <- length(unique(df3$Ensemble))
  
  #Categrory 1: Top subnet genes
  a <- TSOMIM
  b <- length(AllTS) - a
  
  #Category 2: Non-Topsub genes
  #Total Roadmap: 9729 
  allNonTS <- 9729 - length(AllTS)
  c <- allOMIMGenes - TSOMIM
  d <- allNonTS - c
  
  table <- rbind(c(a,b), c(c,d))
  colnames(table) <- c("OMIM","NonOMIM")
  rownames(table) <- c("TS","NonTS")
  
  result <- fisher.test((table) , alternative="greater")
  #cat(class, result$p.value, "\n")    
  contingencyTable <- table
  pvalue <- result$p.value
  
  list(contingencyTable = contingencyTable, pvalue = pvalue)
}