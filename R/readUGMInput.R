#' read UGM inputs
##' @title read UGM inputs
##' @param UGMresultDir UGM result directory
##' @param nNull number of times to calcualte the partition function
##' @return a list of UGM constructs

readUGMInput <- function(UGMresultDir, nNull){

  #Read partionfucntion calculation
  filename <- paste("Partionfunction.txt",sep="")
  file <- paste(UGMresultDir, filename, sep="")
  df <- read.table(file, header = FALSE, stringsAsFactors=FALSE)
  Partitionfunction <- df[,1]

  #Read edgePot and nodeStates files
  edgePotList <- list()
  nodeStatesList <- list()
  for(i in 1:nNull){

    file1 <- paste(UGMresultDir, "edgePot",i,".txt", sep="")
    file2 <- paste(UGMresultDir, "nodeStates",i,".txt", sep="")
    df1 <- read.table(file1, header = FALSE, stringsAsFactors=FALSE)
    df2 <- read.table(file2, header = TRUE, stringsAsFactors=FALSE)
    edgePotList[[i]] <- df1
    nodeStatesList[[i]] <- df2
  }

  list(Partitionfunction = Partitionfunction, edgePotList = edgePotList, nodeStatesList = nodeStatesList)
}
