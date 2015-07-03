#' Keep gene annotation up to date.
##' @title gene symbol conversion
##' @param geneNames Ensemble ID
##' @return geneNames HUGO gene symbol
##' @export
##' @import biomaRt

geneAnnotation <- function(geneNames){
  
  cat("This is yet to complete.\n")
  ensembl <- useMart("ensembl")
  ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")

  Ensemble <- geneNames
  HUGO <- c()
  for (gene in Ensemble){
    result <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"), filters="ensembl_gene_id", values = gene, mart = ensembl)
    if (length(result$hgnc_symbol) != 0 ){
      HUGO <- append(HUGO, result$hgnc_symbol)
    }else{
      HUGO <- append(HUGO, "NA")
    }
  }

  return()
}
