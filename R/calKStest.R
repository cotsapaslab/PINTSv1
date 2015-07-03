#' perform two sample KS test either the statdard KS test or bootstrap version
##' @title two sample KS test
##' @param x A test set
##' @param y A reference set
##' @param option either "standard" or "bootstrap"
##' @param nboots number of bootstrap
##' @param alternative "less" or "greater"
##' @return KS test p values
##' @import Matching

calKStest <- function(x, y, option, nboots=NULL, alternative){

  if (option == "standard"){
    ks1 <- ks.test(x, y, alternative="less")
    ksPval <- ks1$p.value
  }

  if (option == "bootstrap"){
    ks2 <- ks.boot(x, y, nboots=5000, alternative="less")
    ksPval <- ks2$ks.boot.pvalue
  }

  return(ksPval)
}
