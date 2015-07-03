#' Test the significance of association signal clustering (top subnetwork) against null association signal clustering (null top subnetworks)test the significance of association signal clustering
##' @title Test the significance of association signal clustering
##' @param diseaseNetwork a phenotype-specific PPI network (Pajek format graph)
##' @param nodeAnnot Node annotation (data frame)
##' @param detectedSubnets A list of detected nodes for the top subnetworks
##' @param nullNetwork null association disease networks (randomized disease-association network)
##' @param nNull the number of null networks
##' @param sigThres genome-wide association threshold. e.g., sigThres <- 5e-6
##' @return p values of the permutation test for node #, edge #, clustering coefficient, and genetic score
##' @export
##' @import igraph

testSigClustering <- function(diseaseNetwork, nodeAnnot, detectedSubnets, nullNetwork, nNull, sigThres){
  subnetID <- 0
  sigSubnetSummary <- data.frame()

  for (subnet in detectedSubnets[1:(nNull+1)]){
    subnetID <- subnetID + 1
    #print(subnetID)

    #Retreive p values for Subnet
    if (subnetID == 1){
      #Top subnetwork
      nodeNamesIndex <- which( V(diseaseNetwork)$name %in% subnet)
      subnetPvalues <- as.numeric(V(diseaseNetwork)$weight[nodeNamesIndex])
      subnetNames <- V(diseaseNetwork)$name[nodeNamesIndex]
      dftemp <- as.data.frame(cbind(subnetNames, subnetPvalues), stringsAsFactors=F)
      #print(str(dftemp))
      steinerNodes <- dftemp$subnetNames[which( as.numeric(dftemp$subnetPvalues) > sigThres)]
      #print(steinerNodes)
    }else{
      #Null top subnetworks
      nullNetworkID <- paste("null", subnetID-1, sep="")
      randNodeIndex <- as.numeric(nullNetwork[,nullNetworkID])
      pvalues <- V(diseaseNetwork)$weight[randNodeIndex]
      names(pvalues) <- V(diseaseNetwork)$name
      subnetPvalues <- pvalues[subnet]

    }
    subnetScore <- sum(-log10(as.numeric(subnetPvalues)))

    #Induce graph and calculate the network clustering coefficents
    sigSubnet <- induced.subgraph(diseaseNetwork, subnet)
    subnetSize <- length(subnet)
    subnetEdge <- ecount(sigSubnet)
    if (subnetSize == 1) { subnetCC <- 0 }
    else{ subnetCC <- transitivity(sigSubnet) }
    subnetMeanDegree <- mean(degree(diseaseNetwork, subnet))
    subnetDiameter <- diameter(sigSubnet)

    #cat("Node #", subnetSize, "Edge #", subnetEdge, "subnetCC :", subnetCC, "\n")
    #cat("mean degree: ", subnetMeanDegree, "subnetDiameter", subnetDiameter, "\n")
    sigSubnetInfo <- c(subnetSize, subnetEdge, subnetCC, subnetScore, subnetMeanDegree)
    sigSubnetSummary <- rbind(sigSubnetSummary, sigSubnetInfo)
  }

  colnames(sigSubnetSummary) <- c("subnetSize", "subnetEdge", "subnetCC", "subnetScore", "subnetMeanDegree")
  #sigSubnetSummaryPvals <- t(apply(sigSubnetSummary, 2, function(x) 1 - trunc(rank(x))/length(x)))
  sigSubnetSummaryPvals <- t(apply(sigSubnetSummary, 2, function(x)  1 - ( trunc(rank(x,ties.method = c("average"))) - 1)/nNull ) ) 
  #print(sigSubnetSummary)
  #print(sigSubnetSummaryPvals)
  list(sigSubnetSummary, sigSubnetSummaryPvals)
}
