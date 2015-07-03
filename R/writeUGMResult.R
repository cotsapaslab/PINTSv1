#' writeUGMResult
##' @title writeUGMResult
##' @param UGM A list of UGM constructs
##' @param resultDir A directory where the UGM constructs are written to

writeUGMResult <- function(UGM, resultDir){

  cat("The results are written under", resultDir, ".\n")
  for(i in 1:length(UGM$adjMatList)){
    write.table(UGM$adjMatList[[i]], paste(resultDir, "/adjMatrix", i, ".txt", sep=""), sep="\t", quote=F, row.names=F, col.names=F)
    write.table(UGM$edgePotList[[i]], paste(resultDir, "/edgePot", i, ".txt", sep=""), sep="\t", quote=F, row.names=F, col.names=F)
    write.table(UGM$nodePotList[[i]], paste(resultDir, "/nodePot", i, ".txt", sep=""), sep="\t", quote=F, row.names=F, col.names=F)
    write.table(UGM$nodeStatesList[[i]], paste(resultDir, "/nodeStates", i, ".txt", sep=""), sep="\t", quote=F, row.names=F, col.names=T)
  }

}
