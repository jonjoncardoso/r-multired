#' @import Matrix
#' @import plyr
library("Matrix")
library("plyr")

# Read and print multilayer network ---------------------------

#' Read a multilayer network from a text file
#' 
#' The read.multilayer function reads only text files.
#' 
#' The file must contain three or four columns, each line with the following format:
#' 
#'  <layerID> <source> <target> (<weight>)
#'  
#'  If the fourth column is not informed, all edges will be assigned weight 1.
#'  
#' 
#' @param filename The path to the file
#' @param n The total number of nodes in the multilayer network. 
#' If not given, multired will infer the number of nodes from the edge list automatically, 
#' but depending on the size of the network, this might take a while.
#' @param networkName The name of the network (optional).
#' @return A multilayer object
#' @export
read_multilayer <- function(filename, n, networkName=NULL) {
    multilayer <- list()
    multilayer$edgelist <- read.table(filename, header = F)
    if (ncol(multilayer$edgelist) == 4) {
        names(multilayer$edgelist) = c("layerID", "source", "target", "weight")
    } else if (ncol(multilayer$edgelist) == 3) {
        multilayer$edgelist[, 4] <- 1
        names(multilayer$edgelist) = c("layerID", "source", "target", "weight")
    } else {
        stop("Multilayer file must be in the format:\n", "\t<layerID> <source> <target> (<weight>).\n", "")
    }
    multilayer$noLayers <- length(unique(multilayer$edgelist$layerID))
    if (missing(n)) {
        warning("Number of nodes were not informed, inferring from the edge list.",call. = F)
        multilayer$n <- c("noNodes" = length(union(multilayer$edgelist$source, multilayer$edgelist$target)))
    } else {
        multilayer$n <- c("noNodes"=n)
    }
    multilayer$name <- c("name"=networkName)
    class(multilayer) <- "multilayer"
    
    #Store adjacency matrices
    multilayer$adjacencies <- lapply(1:multilayer$noLayers,get_adjacency_of_layer,multilayer)
    
    multilayer
}

#' Converts a List of Matrix to Data frame
#' 
#' Each position in the list must be a n x n matrix
#' 
#' @export
list_to_multilayer <- function(l, n, networkName=NULL){
  multilayer <- list()
  multilayer$edgelist <- NA
  multilayer$noLayers <- length(l)
  if (missing(n)) {
    warning("Number of nodes were not informed, inferring from the list of adjacencies.",call. = F)
    multilayer$n <- c("noNodes" = dim(l[[1]])[1])
  } else {
    multilayer$n <- c("noNodes"=n)
  }
  multilayer$adjacencies <- l
  multilayer$name <- c("name"=networkName)
  class(multilayer) <- "multilayer"
  multilayer
}

#' @export
print <- function(x,...){
  UseMethod("print",x)
}

#' @export
print.multilayer <- function(m,...) {
    invisible(cat(sprintf("Multilayer network %s\n Nodes = %d\n Layers = %d\n", 
                          ifelse(is.null(m$name), ":", m$name), m$n, m$noLayers)))
}

# Functions for multilayer objects ---------------------------

#' Adjacency matrix of one or more layers
#' 
#' The function will return the adjacency matrix of one or more layers of the multilayer network, 
#'  if more than one layerID is passed, the function will aggregate the layers, 
#'  summing the weights of the links.
#'  
#' 
#' @param layerIds An integer or a vector, identifying one or more layerID of the network. 
#' If layerIds is a vector, the layers will be aggregated.
#' @param m A multilayer object
#' @return A \code{\link[Matrix]{Matrix}} object
#' @export
get_adjacency_of_layer <- function(layerId,m){
  if(!inherits(m,"multilayer")){
    stop("m must be a multilayer object.")
  }
    A <- Matrix(0, nrow = m$n, ncol = m$n, sparse = TRUE)
    if(layerId <= m$noLayers & layerId > 0){
      el <- subset(m$edgelist, layerID == layerId)
      # Summing weights in all layers
      el <- ddply(el, .(source, target), summarize, weight = sum(weight))
      el <- as.matrix(el)
      tryCatch({
        A[el[, 1:2]] <- el[, "weight"]
        # Sum existing value with the weight from the other link (make it bi-directional)
        A[el[, c(2, 1)]] <- A[el[, c(2, 1)]] + A[el[, c(1, 2)]]
      }, error = function(e) {
        stop("Could not create adjacency matrix. Indices of nodes >= number of nodes informed.")
      })
      A
    }else{
     stop("Invalid layerId.") 
    }  
} 

#' Adjacency matrix of one or more layers
#' 
#' The function will return the adjacency matrix of one or more layers of the multilayer network, 
#'  if more than one layerID is passed, the function will aggregate the layers, 
#'  summing the weights of the links.
#'  
#' 
#' @param A adjacency matrix.
#' @return The entropy of the adjacency matrix.
#' @export
get_layer_entropy <- function(A){
  D <- Diagonal(nrow(A),rowSums(A)) 
  L <- (D-A)
  L <- L/sum(diag(L))
  lambda <- eigen(L,only.values=T)$values
  - Reduce("+",
             sapply(lambda,function(l){if(l > 10e-20){l * log(l)}else{0}}))
  
}

#' @export
get_entropy <- function(m){
  if(!inherits(m,"multilayer")){
    stop("m must be a multilayer object.")
  }
  Reduce("+",Map(get_layer_entropy,m$adjacencies))
}

#' Entropy of the aggregation of all layers
#' 
#' @export
get_agg_entropy <- function(m){
  if(!inherits(m,"multilayer")){
    stop("m must be a multilayer object.")
  }
  agg <- Reduce("+",m$adjacencies)
  get_layer_entropy(agg)
}

#' Computes Jensen-Shannon matrix
#'  @export
compute_JSD_matrix <- function(m){
  if(!inherits(m,"multilayer")){
    stop("m must be a multilayer object.")
  }
  allComb <- t(combn(m$noLayers,2))
  JSD <- Matrix(0,m$noLayers,m$noLayers)
  for(idx in 1:nrow(allComb)){
    i <- allComb[idx,1]
    j <- allComb[idx,2]
    aggr <- (m$adjacencies[[i]] + m$adjacencies[[j]])/2
    JSD[i,j] <- sqrt(get_layer_entropy(aggr) - 0.5*(get_layer_entropy(m$adjacencies[[i]]) + get_layer_entropy(m$adjacencies[[j]])))
    JSD[j,i] <- JSD[i,j]
  }
  JSD
}

#' Calculates relative entropy (q) for reduced multilayer network
#' 
#' Calculates the relative entropy of a reduced multilayer, 
#' relative to the aggregated entropy of the original multilayer (not the reduced). 
#' 
#' @param reducedM the multilayer network with reduced number of layers
#' @param aggEntropy the aggregated entropy of the original multilayer. See \code{\link{get_agg_entropy}}.
#' @return The relative entropy (q)
#' @export
get_relative_entropy <- function(reducedM,aggEntropy){
  if(!inherits(reducedM,"multilayer")){
    stop("reducedM must be a multilayer object.")
  }
  HC <- get_entropy(reducedM)/reducedM$noLayers
  1 - HC/aggEntropy
}

#' Find the reduced multilayer network
#' 
#' \code{find_reduced} is the main function of this package. 
#' It reduces the multilayer network passed as an argument (\code{m}) and  
#' 
#' @export
#' @examples
#' \dontrun{
#' m <- read_multilayer("test.txt")
#' find_reduced(m)
#' }
find_reduced <- function(m){
  if(!inherits(m,"multilayer")){
    stop("m must be a multilayer object.")
  }
  aggEntropy <- get_agg_entropy(m)
  
  merges <- hclust(as.dist(compute_JSD_matrix(m)),method="ward.D2")$merge
  noMerges <- nrow(merges)
  
  qVals <- vector(mode="numeric",length=noMerges)
  adjMerges <- vector(mode="list",length=noMerges-1)
  currMultilayer <- m$adjacencies
  merged <- vector(mode="logical",length=m$noLayers)
  
  #Evaluate all merges
  for(i in 1:noMerges){
    layer1 <- merges[i,1]
    layer2 <- merges[i,2]
    
    if(layer1 < 0){
      merged[abs(layer1)] <- TRUE
    }
    if(layer2 < 0){
      merged[abs(layer2)] <- TRUE
    }
    
    adj1 <- if(layer1 < 0){m$adjacencies[[abs(layer1)]]}else{adjMerges[layer1][[1]]}
    adj2 <- if(layer2 < 0){m$adjacencies[[abs(layer2)]]}else{adjMerges[layer2][[1]]}
    adjMerges[[i]] <- adj1 + adj2
    currMultilayer <- Reduce(append,m$adjacencies[which(!merged)],adjMerges[1:i])
    
    reducedM <- list_to_multilayer(currMultilayer,m$n,m$name)
    qVals[i] <- get_relative_entropy(reducedM,aggEntropy)
  }
  qVals
}

