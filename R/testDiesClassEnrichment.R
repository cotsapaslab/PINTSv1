#' test disease class enrichment of gene sets
##' @title test disease class enrichment of gene sets
##' @param geneDiseaseclass a data frame for gene-Disease class matrix
##' @param topSubnet top subnetwork genes
##' @return p values of Fisher's exact test and contingency table
##' @export

testDiesClassEnrichment <- function(geneDiseaseclass, topSubnet){
  
  dfRef <- geneDiseaseclass 
  colSumsRef <- colSums(dfRef[,-c(1)])
  nGenesRef <- dim(dfRef)[1]
  
  dfTS <- dfRef[dfRef$Ensemble %in% topSubnet, ]
  colSumsTS <- colSums(dfTS[,-c(1)])
  nGenesTS <- dim(dfTS)[1]
  
  contingencyTables <- list()
  pvalueList <- list()
  
  for (class in names(dfRef)[-c(1)]){
    
    if(as.numeric(colSumsTS[class]) > 0) {
      #Categrory 1: OMIM recorded vs in a Disease class for top subnetwork
      a <- colSumsTS[class] #Classified in TS
      b <- nGenesTS - colSumsTS[class] #Not-classified in TS
      
      nGenesNonTS <- nGenesRef - nGenesTS
      #Category 2: OMIM recorded vs in a Disease class for entire gene set. 
      #c <- nGenesRef - a #Not-classified outside TS
      c <- colSumsRef[class] - a
      d <- nGenesNonTS - (colSumsRef[class] - a) #Not-classified outside TS
      
      table <- rbind(c(a,b), c(c,d))
      colnames(table) <- c("Class","NonClass")
      rownames(table) <- c("TS","NonTS")
      #print(table)
      result <- fisher.test((table) , alternative="greater")
      #cat(class, result$p.value, "\n")
      
      contingencyTables[[class]] <- table
      pvalueList[[class]] <- result$p.value
      
    }
  }
  
  list(contingencyTables = contingencyTables, pvalues = pvalueList)
}