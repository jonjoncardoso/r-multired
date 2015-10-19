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
    
    multilayer$adjacencies <- lapply(1:multilayer$noLayers,get_adjacency,multilayer)
    
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
get_adjacency <- function(layerIds,m){
  UseMethod("get_adjacency",m)
}

#' @export
get_adjacency.multilayer <- function(layerIds, m) {
    A <- Matrix(0, nrow = m$n, ncol = m$n, sparse = TRUE)
    el <- subset(m$edgelist, layerID %in% layerIds)
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
  UseMethod("get_entropy",m)
}

get_entropy.multilayer <- function(m){
  Reduce("+",Map(get_layer_entropy,m$adjacencies))
}

