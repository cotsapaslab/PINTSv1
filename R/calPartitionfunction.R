#' calculate partition function of a UGM (Undirected Graphical Model)
##' @title calculate partiontion function
##' @param nNull the number of times to calculate the partition function
##' @param resultDir a dirctory where the result is saved
##' @param UGMPATH a path for the UGM software
##' @param MATLABPATH a path for MATLAB software
##' @return a set of calculated partition function

calPartitionfunction <- function(nNull, resultDir, UGMPATH, MATLABPATH){

  #Matlab path and UGM path
  #UGMpath <- "/home/jc2437/Complex_disease_network/Algorithm/UGM"

  #Parameters
  nNull <- nNull
  UGMInDir <- paste(resultDir, sep="")
  UGMOutDir <- paste(resultDir, sep="")

  #runUGM
  UGMbash <- paste(UGMPATH, 'runUGM.sh', sep="")
  command=paste(UGMbash, nNull, UGMInDir, UGMOutDir, UGMPATH, MATLABPATH, sep=" ")
  system(command)

}
