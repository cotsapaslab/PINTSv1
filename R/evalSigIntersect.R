#' evaluate the significance of the intersect of the two sampled sets from a total population.
#'
##' @title evaluate intersect significance
##' @param set1TotalSize The size of the first set
##' @param set2TotalSize The size of the second set
##' @param intersectSize The sized of intersect between the two sets
##' @param TotalPopulation The size of total population where the two sets are sampled
##' @return p value of hypergeometic test


evalSigIntersect <- function(set1TotalSize, set2TotalSize, intersectSize, TotalPopulation) {

  #Hypergeometric test
  return(phyper(q=intersectSize, m=set1TotalSize, k=set2TotalSize, n=TotalPopulation-set1TotalSize, lower.tail=F))

}
