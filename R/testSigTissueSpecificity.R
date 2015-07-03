#' test the significance of tissue specificity
#'
##' @title Tissue specificity analysis
##' @param diseaseNetwork a phenotype-specific PPI network (Pajek format graph)
##' @param geneExpression a data frame tissue specificity.
##' @param topSubnets a list of the detected subnetworks
##' @param nNull the number of null top subnetworks
##' @param permuteOption  "nullTopsubnet" or "randGeneExpression"
##' @param algorithmOption UGMPath (UGM software installation directory)
##' @param MATLABPATH MATLABPATH (Matlab installation directory)
##' @return p values
##' @references UGM software
##' @export

testSigTissueSpecificity <- function(diseaseNetwork, geneExpression, topSubnets, nNull, permuteOption, algorithmOption=NULL, MATLABPATH=NULL){

  nNull <- nNull + 1
  cat("run tissue specificity test\n ")
  workingDir <- getwd()
  if(permuteOption == "nullTopsubnet"){
    dtime <- format(Sys.time(), "%H%M%S-%b%d%Y")
    resultDir <- paste(workingDir,"/UGMresultOP1_",dtime,"/",sep="")
    cat("Null model: Null Topsubnetwork\n")
    UGM <- UGMConstruct1(diseaseNetwork, geneExpression, topSubnets, nNull)

    if (file.exists(resultDir)){
      writeUGMResult(UGM, resultDir)
    } else {
      dir.create(resultDir)
      writeUGMResult(UGM, resultDir)
    }

    UGMPATH <- algorithmOption
    if ( !is.null(MATLABPATH) & !is.null(UGMPATH) ){
        cat("calculate the partitionfunction with UGM software \n")
        calPartitionfunction(nNull, resultDir, UGMPATH, MATLABPATH)
        UGMinputs <- readUGMInput(resultDir, nNull)
    }else{
      #Need to check if the file exist and provide an instruction.
      cat("If MATLAB/UGM are not available in run time,
          you can provide a partition function calculation and place it
          under", resultDir, ".\n")
      UGMinputs <- readUGMInput(resultDir, nNull)
    }

    UGMJointPDF <- calJointPDF1(diseaseNetwork, topSubnets, UGMinputs$edgePotList, UGMinputs$nodeStatesList, UGMinputs$Partitionfunction, nNull)

  }else if(permuteOption == "randGeneExpression" ) {
    cat("Null model: randomized tissue specifcity matrix. \n")
    #Randomize the tissue specificity matrix while preserving the tissue specificity of each gene.
    #For example, gene1 is tissue specific in a set of tissues {T_i}
    #In randomization, gene 1 is randomly reassigned to another set of tissues {T_j} 
    #while the number of specific tissues is kept the same (#{T_i} = #{T_j}).  
    
    dtime <- format(Sys.time(), "%H%M%S-%b%d%Y")
    resultDir <- paste(workingDir,"/UGMresultOP2_",dtime,"/",sep="")
    UGM <- UGMConstruct2(diseaseNetwork, geneExpression, topSubnets, nNull)

    if (file.exists(resultDir)){
      writeUGMResult(UGM, resultDir)
    } else {
      dir.create(resultDir)
      writeUGMResult(UGM, resultDir)
    }

    UGMPATH <- algorithmOption
    if (!is.null(MATLABPATH) & !is.null(UGMPATH)){
      cat("calculate the partitionfunction with UGM software \n")
      calPartitionfunction(nNull, resultDir, UGMPATH, MATLABPATH)
      UGMinputs <- readUGMInput(resultDir, nNull)
    }else{
      #Need to check if the necessary files exist and then provide an instruction.
      cat("If MATLAB/UGM are not available in run time,
         you can provide a partition function calculation and place it
          under", resultDir, ".\n")
      UGMinputs <- readUGMInput(resultDir, nNull)
    }

    UGMJointPDF <- calJointPDF2(diseaseNetwork, topSubnets, UGMinputs$edgePotList, UGMinputs$nodeStatesList, UGMinputs$Partitionfunction, nNull)
  }

  return(UGMJointPDF)
}
