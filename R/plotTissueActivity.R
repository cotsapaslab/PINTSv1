#' plot the tissue activity
##' @title plot the tissue activity
##' @param diseaseNetwork a phenotype-specific PPI network (Pajek format graph)
##' @param geneExpression a data frame for the tissue specificity
##' @param targetCells A set of cells/tissues under test
##' @param topSubnets top subnetwork genes 
##' @return A network plot showing the tissue specific genes in a top subnetwork for target tissues
##' @export
##' @import igraph

plotTissueActivity <- function(diseaseNetwork, geneExpression, targetCells, topSubnets){

  set.seed(123456)
  subnet <- topSubnets
  gE <- geneExpression
  All.cells <- colnames(gE)[-c(1)]

  #subnet geneExpression (tissue specific)
  subnet.gE <- gE[gE[,"Ensemble"] %in% subnet, ]
  colnames(subnet.gE)[1] <- c("gene")

  #nodeStates of a subnet
  sigSubnet <- induced.subgraph(diseaseNetwork, subnet)
  nodeStates <- subnet.gE[ match(V(sigSubnet)$name, subnet.gE$gene), ]
  All.nodes <- nodeStates$gene

  connNodes <- list()
  sigSubnet.layout <- layout.fruchterman.reingold(sigSubnet)
  
  #pdf("tissueActivityPathFigures.pdf")
  cell.id <- 0
  for (cell in targetCells){

    cell.id <- 1
    #Identify active nodes in a cell
    nodes.On <- nodeStates[which(nodeStates[,cell] == 1), "gene"]
    nodes.On.index <- which(V(sigSubnet)$name %in% nodes.On)

    #find all the SPD nodes that connects the active nodes
    SPDs <- calSPD(sigSubnet, nodes.On)
    AllPathNodes <- SPDs$AllPathNodes

    #gSubComp <- induced.subgraph(sigSubnet, AllPathNodes)
    AllPathNodes.index <- which(V(sigSubnet)$name %in% AllPathNodes)

    #common connector genes
    #commonConn <- c("ENSG00000187555","ENSG00000150991","ENSG00000196924","ENSG00000109971", "ENSG00000168036")
    #commonConn <- c("ENSG00000087460", "ENSG00000196924", "ENSG00000150991", "ENSG00000187555", "ENSG00000109971", "ENSG00000168036", "ENSG00000181222", "ENSG00000072501")
    #commonConn.index <- which(V(sigSubnet)$name %in% commonConn)

    #steiner nodes
    #steinerNodes <-getSteinerNodes(nodeWeight,  subnet, sigThres)
    steinerNodes<- c("ENSG00000150991", "ENSG00000145675", "ENSG00000109971", "ENSG00000171608", "ENSG00000168036")
    steinerNodes.index <- which(V(sigSubnet)$name %in% steinerNodes)

    #Color/shape coding
    V(sigSubnet)$color <- c("white")
    V(sigSubnet)$frame.color <- c("black")
    V(sigSubnet)$shape <- "circle"

    #Steiner node
    #V(sigSubnet)$frame.color[steinerNodes.index] <- c("blue")
    V(sigSubnet)$shape[steinerNodes.index] <- "square"

    #Connector nodes
    V(sigSubnet)$color[setdiff(AllPathNodes.index, nodes.On.index) ] <- c("grey")

    #Common connectors in two tissues
    #V(sigSubnet)$color[commonConn.index] <- c("gray48")
    #V(sigSubnet)$frame.color[commonConn.index] <- c("red")
    V(sigSubnet)$color[nodes.On.index] <- c("black")

    #Edge coloring
    sigSubnetEdges <- get.edgelist(sigSubnet, names=TRUE)
    for ( edgeCount in 1:dim(sigSubnetEdges)[1] ){
      node1 <- sigSubnetEdges[edgeCount, 1]
      node2 <- sigSubnetEdges[edgeCount, 2]
      if( (node1 %in% AllPathNodes) && (node2 %in% AllPathNodes) ){
          E(sigSubnet)$color[edgeCount] <- c("black")
          E(sigSubnet)$lty[edgeCount] <- c(1)
     }else{
        E(sigSubnet)$color[edgeCount] <- c("gray")
        E(sigSubnet)$lty[edgeCount] <- c(2)
     }
    }

    connNodes[[cell.id]] <- V(sigSubnet)$name[setdiff(AllPathNodes.index, nodes.On.index)]
    #lapply(connNodes, cat, "\n", file = "./ConnectorNodes.txt", append=TRUE)

    plot(sigSubnet, layout = sigSubnet.layout, main= cell, vertex.color=V(sigSubnet)$color, vertex.shape=V(sigSubnet)$shape, vertex.size= 7, vertex.label="", vertex.label.font=3, vertex.label.cex=0.5, edge.color=E(sigSubnet)$color, edge.lty=E(sigSubnet)$lty, edge.width=2)
    #plot(gSubComp, vertex.color="black", main= cell, vertex.size=5, vertex.label="", vertex.label.cex=0.5, edge.width=2 )
  }

  #gSubComp.layout <- layout.fruchterman.reingold(gSubComp)
  #clusters(gSubComp)
  #plot(gSubComp, layout = gSubComp.layout)
  #dev.off()
  return()
}

