## This file is part of triversity.
##
## triversity is an R package for the computation of diversity measures on
## tripartite graphs. It has been developed by researchers of the Complex
## Networks team, within the Computer Science Laboratory of Paris 6 (LIP6),
## for the ALGODIV project, founded by the French National Agency of
## Research (ANR) under grant ANR-15-CE38-0001.
## 
## Copyright © 2017 Robin Lamarche-Perrin (<Robin.Lamarche-Perrin@lip6.fr>)
## 
## triversity is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by the
## Free Software Foundation, either version 3 of the License, or (at your
## option) any later version.
## 
## triversity is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
## more details.
## 
## You should have received a copy of the GNU General Public License along
## with this program. If not, see <http://www.gnu.org/licenses/>.


#' @title Compute diversity measures on multipartite graphs.
#'
#' @description
#' \code{triversity} is an R package for the computation of diversity measures
#' on multipartite graphs. First, it implements a parametrized family of such diversity
#' measures which apply on probability distributions. Sometimes called "True Diversity",
#' this family contains famous measures such as the richness, the Shannon entropy, the
#' Herfindahl-Hirschman index, and the Berger-Parker index. Second, the package allows
#' to apply these measures on probability distributions resulting from random walks between
#' the parts of multipartite graphs. By defining an initial distribution on a given part
#' of the graph and a path to follow between the different parts, the probability of the
#' walker's position within the final part is then computed, thus providing a particular
#' instance of diversity to measure.
#'
#' This package has been developed by researchers of the Complex Networks team,
#' located within the Computer Science Laboratory of Paris 6 (LIP6),
#' for the AlgoDiv project founded by the French National Agency of Research
#' (ANR) under grant ANR-15-CE38-0001.
#'
#' Links:
#' \itemize{
#' \item AlgoDiv project: \url{http://algodiv.huma-num.fr/}
#' \item Complex Networks team: \url{http://www.complexnetworks.fr/}
#' \item LIP6: \url{https://www.lip6.fr/}
#' \item ANR: \url{http://www.agence-nationale-recherche.fr/}
#' }
#' 
#' Contact:
#' \itemize{
#' \item Robin Lamarche-Perrin: \email{Robin.Lamarche-Perrin@@lip6.fr}
#' 
#' See also my webpage: \url{https://www-complexnetworks.lip6.fr/~lamarche/}
#' }
#' 
#' 
#' List of main collaborators:
#' \itemize{
#' \item Robin Lamarche-Perrin (CNRS, ISC-PIF, LIP6)
#' \item Lionel Tabourier (UPMC, LIP6)
#' \item Fabien Tarissan (CNRS, ISP, LIP6)
#' \item Raphaël Fournier S'niehotta (CNAM, CÉDRIC)
#' \item Rémy Cazabet (UdL, LIRIS)
#' }
#' 
#' Copyright 2017 © Robin Lamarche-Perrin
#' 
#' \code{triversity} is free software: you can redistribute it and/or modify it under the
#' terms of the GNU General Public License as published by the Free Software Foundation,
#' either version 3 of the License, or (at your option) any later version. It is distributed
#' in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
#' warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#' Public License for more details. You should have received a copy of the GNU General
#' Public License along with this program. If not, see \url{http://www.gnu.org/licenses/}.
#' 
#' 
#' @docType package
#' 
#' @name triversity
#' 
#' @import stats slam data.tree
NULL


triversity.env <- new.env (parent = emptyenv())
triversity.env$distribution_precision <- 1e-5
triversity.env$verbose <- TRUE


#' @title Build a properly-structured multipartite graph from raw data.
#' 
#' @description
#' \code{get_multipartite} builds a properly-structured multipartite graph
#' from a file or from a \code{data.frame} containing raw data. This object can then
#' be used by the other functions of this package. The structure
#' of the input data and of the resulting structure is detailed below.
#'
#' @param filename A character string giving the path to the file containing the raw data.
#'
#' The input file should contain at least four columns, separated by spaces. Each row gives
#' a link between two nodes belonging to two different parts of the multipartite graph.
#' The first column gives the label of the part of the first node and the second column
#' gives the label of this node (any character string). Similarly, the third and fourth
#' columns give the labels of the second part and of the second node. The file can eventually
#' contain a fifth column giving the weights of the links (any positive integer or float
#' value).
#' 
#' @param data A \code{data.frame} containing the raw data.
#'
#' This \code{data.frame} should have the same structure than the one described above
#' for the case of an input file: four columns indicating the labels of the
#' parts and of the nodes constituting the link, and an optional
#' fifth column for its weight. At least \code{filename} or \code{data} should be
#' specified, but not both at the same time.
#' 
#' @return A properly-structured multipartite graph that can be used by the other functions
#' of the \code{triversity} package.
#' 
#' The resulting object encodes the labels of the parts and of the nodes in the multipartite
#' graph, as well as the transition matrices of random walks following different paths
#' between parts (encoded as sparse matrices of floats in [\code{0},\code{1}], with all
#' rows summing to \code{1}). These transition matrices
#' are then used by functions such as \code{\link{get_distribution_from_path}} to
#' compute the probability distributions of such random walks, or such as
#' \code{\link{get_diversity_from_path}} to compute the diversity of these distributions.
#'
#' Assuming the object returned by \code{get_multipartite} is stored in the
#' \code{graph} variable, then:
#' \itemize{
#' \item \code{graph$parts} is a vector of character strings giving the labels of the
#' different parts in the multipartite graph;
#' \item \code{graph$nodes} is a list of vectors of character strings giving the labels
#' of the nodes constituting the different parts of the multipartite graph;
#' \item \code{graph$edges} is a \code{data.tree} whose nodes each contains the
#' transition matrix of the corresponding random walk.
#' }
#'
#' @seealso Package \code{slam} for the sparse encoding of transition matrices and
#' package \code{data.tree} for the tree structure storing these matrices.
#' 
#' @examples
#' data (tripartite_example)
#' graph <- get_multipartite (data=tripartite_example)
#' 
#' graph$parts
#' 
#' part1 <- graph$parts[1]
#' graph$nodes[[part1]]
#' 
#' part2 <- graph$parts[2]
#' as.matrix (graph$edges$Climb(part1,part2)$matrix)
#'
#' @export
get_multipartite <- function (filename=NULL, data=NULL)
{
    if (is.null (filename) && is.null (data)) { stop ("'filename' or 'data' has to be specified") }
    if (!is.null (filename) && !is.null (data)) { stop ("'filename' and 'data' cannot be both specified at the same time") }

    ## LOAD DATA
    if (triversity.env$verbose) { message ('Loading data') }
    time <- Sys.time()
    total_time <- time

    if (!is.null (filename)) { data <- utils::read.table (filename, stringsAsFactors=FALSE) }
    if (length(colnames(data)) < 4) { stop ("input data should have at least 4 columns") }
    
    colnames(data)[1:4] <- c ('part_from', 'node_from', 'part_to', 'node_to')
    if (ncol(data) > 4) { colnames(data)[5] <- 'weight' } else { data$weight <- 1 }

    if (triversity.env$verbose) { message ('-> done in ', round (as.numeric (Sys.time() - time, units="secs"), 3), ' sec') }

    ## BUILD MULTIPARTITE GRAPH
    graph <- list()

    ## BUILD PARTS
    if (triversity.env$verbose) { message ('Building parts') }
    time <- Sys.time()
    
    graph$parts <- unique (append (unique (as.character (data$part_from)), unique (as.character (data$part_to))))
    graph$parts <- graph$parts[order(graph$parts)]

    if (triversity.env$verbose) { message ('-> done in ', round (as.numeric (Sys.time() - time, units="secs"), 3), ' sec') }

    ## BUILD NODES
    if (triversity.env$verbose) { message ('Building nodes') }
    time <- Sys.time()
    
    graph$nodes <- sapply (graph$parts, function (part) {
        nodes <- unique (append (as.character (data[data$part_from == part, 'node_from']), as.character (data[data$part_to == part, 'node_to'])))
        return (nodes[order(nodes)])
    }, simplify=FALSE, USE.NAMES=TRUE)

    if (triversity.env$verbose) { message ('-> done in ', round (as.numeric (Sys.time() - time, units="secs"), 3), ' sec') }

    ## BUILD TREE
    if (triversity.env$verbose) { message ('Building data tree') }
    time <- Sys.time()

    graph$edges <- Node$new ('tree')

    for (part_from in graph$parts) {
        tree_from <- graph$edges$AddChild (part_from)

        nodeNb <- length(graph$nodes[[part_from]])
        i <- 1:nodeNb
        v <- rep(1,nodeNb)
        tree_from$matrix <- simple_triplet_matrix (i, i, v, dimnames=list (graph$nodes[[part_from]], graph$nodes[[part_from]]))

        for (part_to in graph$parts) { tree_from_to <- tree_from$AddChild (part_to) }
    }

    if (triversity.env$verbose) { message ('-> done in ', round (as.numeric (Sys.time() - time, units="secs"), 3), ' sec') }

    ## BUILD EDGES
    for (p_from in 1:length(graph$parts)) {
        part_from <- graph$parts[p_from]

        for (p_to in p_from:length(graph$parts)) {
            part_to <- graph$parts[p_to]
            tree_from_to <- graph$edges$Climb (c (part_from, part_to))
            tree_to_from<- graph$edges$Climb (c (part_to, part_from))
            
            if (triversity.env$verbose) { message ('Building edges from ', part_from, ' to ', part_to) }
            time <- Sys.time()

            data_from_to <- data[data$part_from == part_from & data$part_to == part_to, c(1,2,3,4,5)]
            data_to_from <- data[data$part_from == part_to & data$part_to == part_from, c(3,4,1,2,5)]
            colnames (data_to_from) <- colnames (data_from_to)

            if (nrow(data_from_to) > 0 || nrow(data_to_from) > 0) {
                data_from_to <- rbind (data_from_to, data_to_from)
                data_from_to <- aggregate (list (weight=data_from_to$weight), by=list (node_from=data_from_to$node_from, node_to=data_from_to$node_to), FUN=sum)

                i <- match (data_from_to[,'node_from'], graph$nodes[[part_from]])
                j <- match (data_from_to[,'node_to'], graph$nodes[[part_to]])
                v <- data_from_to$weight
                
                tree_from_to$matrix <- simple_triplet_matrix (i, j, v, dimnames=list (graph$nodes[[part_from]], graph$nodes[[part_to]]))
                tree_from_to$matrix <- tree_from_to$matrix / row_sums (tree_from_to$matrix)

                if (part_from != part_to) {
                    tree_to_from$matrix <- simple_triplet_matrix (j, i, v, dimnames=list (graph$nodes[[part_to]], graph$nodes[[part_from]]))
                    tree_to_from$matrix <- tree_to_from$matrix / row_sums (tree_to_from$matrix)
                }

            } else { tree_from_to$matrix <- NULL }
            
            if (triversity.env$verbose) { message ('-> done in ', round (as.numeric (Sys.time() - time, units="secs"), 3), ' sec') }
        }
    }

    if (triversity.env$verbose) { message ('-> total time ', round (as.numeric (Sys.time() - total_time, units="secs"), 3), ' sec') }

    return (graph)
}

#' @title Build a properly-structured multipartite graph from a \code{data.frame} or a list of
#' \code{data.frame}s.
#'
#' @description
#' \code{get_multipartite_from_df} is a wrapper function that properly formats a a
#' \code{data.frame} or a list of \code{data.frame}s, then calls \code{get_multipartite}.
#'
#' @param data A a \code{data.frame} or a list of \code{data.frame}s containing the raw data.
#'
#' If it is a single \code{data.frame}, it must contain all the variables listed in colnames.
#'
#' If it is a list of \code{data.frame}s, each one must contain a couple of adjacent level
#' variables: the first \code{data;frame} contains levels 1 and 2, the second level 2 and 3,
#' etc.
#'
#' @param colnames A character vector containing the name of the variables for each level of
#' the graph, in order (n^th name is level n).
#'
#' @param weights A logical vector of length equal to number of levels minus one. Specifies
#' whether edges should be weighted (TRUE) or not (FALSE). Weights are computed with the number
#' of rows. By default, the edges between level 1 and level 2 nodes are weighted, while the
#' rest are not.
#'
#' @return A properly-structured multipartite graph that can be used by the other functions
#' of the \code{triversity} package. See \code{\link{get_multipartite}} for details.
#'
#' @seealso \code{\link{get_multipartite}}.
#'
#' @examples
#' ex_data <- data.frame(lev1 = c(1, 1, 2, 2, 3, 3),
#'                       lev2 = c("a", "b", "b", "c", "a", "c"),
#'                       lev3 = c("x", "x", "x", "y", "x", "y"))
#'
#' graph <- get_multipartite_from_df(data = ex_data,
#'                                   colnames = c("lev1", "lev2", "lev3"))
#'
#' ex_data <- data.frame(lev1 = c(1, 1, 2, 2, 3, 3),
#'                       lev2 = c("a", "b", "b", "c", "a", "c"))
#' ex_dict <- data.frame(lev2 = c("a", "b", "c"),
#'                       lev3 = c("x", "x", "y"))
#'
#' graph <- get_multipartite_from_df(data = list(ex_data, ex_dict),
#'                                   colnames = c("lev1", "lev2", "lev3"))
#'
#' @export
get_multipartite_from_df <- function(data, colnames, weights = NULL){
  if(is.null(weights)){
    weights <- c(TRUE, rep(FALSE, length(colnames)-2))
  }
  r <- vector("list", length(colnames) - 1)

  if(is.data.frame(data)){
    data <- data[, colnames]
    data[] <- lapply(data, function(x) {if(!is.character(x)) as.character(x) else x})
    for(i in 1:(length(colnames) - 1)){
      r[[i]] <- aggregate(data[, colnames[i]], data[, colnames[i:(i+1)]], FUN = length)
      names(r[[i]]) <- c("node_from", "node_to", "weight")
      if(!weights[i]){
        r[[i]]$weight <- 1
      }
      r[[i]]$part_from <- i
      r[[i]]$part_to <- i + 1
      r[[i]] <- r[[i]][, c("part_from", "node_from", "part_to", "node_to", "weight")]
    }

  } else if(is.list(data)){
    for(i in 1:(length(colnames) - 1)){
      data[[i]] <- data[[i]][, colnames[i:(i+1)]]
      data[[i]][] <- lapply(data[[i]], function(x) {if(!is.character(x)) as.character(x) else x})
      r[[i]] <- aggregate(data[[i]][, colnames[i]], data[[i]][, colnames[i:(i+1)]], FUN = length)
      names(r[[i]]) <- c("node_from", "node_to", "weight")
      if(!weights[i]){
        r[[i]]$weight <- 1
      }
      r[[i]]$part_from <- i
      r[[i]]$part_to <- i + 1
      r[[i]] <- r[[i]][, c("part_from", "node_from", "part_to", "node_to", "weight")]
    }
  }

  r <- do.call("rbind", r)
  x <- get_multipartite(data = r)
  return(x)
}

#' @title Compute the transition matrix of a random walk following a path between the parts
#' of a multipartite graph.
#'
#' @description
#' \code{get_transition_from_path} computes the transition matrix of a random walk following
#' a \code{path} between the different parts of the input multipartite \code{graph}.
#'
#' @param graph A multipartite graph obtained by calling the \code{\link{get_multipartite}}
#' function.
#'
#' @param path A vector of character strings giving the path that the random walk should
#' follow between the different parts of the input multipartite \code{graph}. This path can
#' be as long as wanted, with eventual cycles, and each string it contains should refer to a
#' label in \code{graph$parts}.
#'
#' @return A matrix of floats in [\code{0},\code{1}], with all lines summing to \code{1},
#' giving the transition matrix of the random walk following the input \code{path}.
#'
#' @details
#' Note that the multipartite graph structure implemented in this package
#' stores in memory any computed transition matrix to avoid redundant computation
#' in the future. Hence, the latter execution of \code{get_transition_from_path}, or of
#' any other function that builds on it, can be much faster than the first call.
#' The transition matrices are stored within a \code{data.tree} in the input
#' \code{graph} variable (see \code{graph$edges}).
#'
#' @examples
#' data (tripartite_example)
#' graph <- get_multipartite (data=tripartite_example)
#'
#' as.matrix (get_transition_from_path (graph, path=c(2,1,2,3)))
#' 
#' @export
get_transition_from_path <- function (graph, path)
{
    path <- as.character (path)
    if (!all (path %in% graph$parts)) { stop ("'path' contains unknown part") }

    compute_transition_from_path (graph, path)    
    return (graph$edges$Climb(path)$matrix)
}

compute_transition_from_path <- function (graph, path)
{
    path <- as.character (path)
    if (!all (path %in% graph$parts)) { stop ("'path' contains unknown part") }

    current_step <- 1
    current_node <- graph$edges$Climb (path[current_step])

    while (current_step < length (path))
    {
        past_step <- current_step
        past_node <- current_node

        current_step <- current_step + 1
        current_node <- current_node$Climb (path[current_step])
        
        if (is.null (current_node))
        {
            current_node <- past_node$AddChild (path[current_step])
            if (is.null (past_node$matrix) || is.null (graph$edges$Climb (path[c(past_step,current_step)])$matrix)) { current_node$matrix <- NULL }
            else { current_node$matrix <- as.simple_triplet_matrix (tcrossprod_simple_triplet_matrix (past_node$matrix, t(graph$edges$Climb (path[c(past_step,current_step)])$matrix))) }
        }
    }
}

    


#' @title Compute the probability distribution associated to a random walk following a path
#' between the parts of a multipartite graph.
#'
#' @description
#' \code{get_distribution_from_path} computes the probability distribution of a random walk
#' following a given \code{path} between the different parts of the input multipartite
#' \code{graph}. It starts at a given part with an initial probability distribution, then
#' randomly follows the links of the graph between the different parts according to the
#' input \code{path}, then stops at the last specified part.
#' 
#' @param graph A multipartite graph obtained by calling the \code{\link{get_multipartite}}
#' function.
#'
#' @param path A vector of character strings giving the path that the random walk should
#' follow between the different parts of the input multipartite \code{graph}. This path can
#' be as long as wanted, with eventual cycles, and each string it contains should refer to a
#' label in \code{graph$parts}.
#' 
#' @param initial_distribution (optional) A vector of floats in [\code{0},\code{1}] and
#' summing to \code{1} giving the probability distribution to start with at the first part
#' of the input \code{path}. It should hence contain as many values as there are nodes in the
#' corresponding part. If not specified, this distribution is assumed uniform.
#'
#' @param initial_node (optional) A character string giving the label of a node in the first
#' part of the input \code{path}. This node is then considered to have probability one, thus
#' being equivalent to specifying an \code{initial_distribution} with only zeros except for
#' one node. If not specified, no such node is defined and the initial distribution is then
#' assumed uniform.
#'
#' @return A vector of floats in [\code{0},\code{1}] and summing to \code{1} giving the
#' probability distribution of the random walk when arriving at the last part, after having
#' followed the input \code{path} within the different parts of the graph.
#'
#' @examples
#' data (tripartite_example)
#' graph <- get_multipartite (data=tripartite_example)
#'
#' path <- c(2,1,2,3)
#' as.matrix (get_distribution_from_path (graph, path))
#' as.matrix (get_distribution_from_path (graph, path, initial_distribution=c(1/3,0,0,2/3)))
#' as.matrix (get_distribution_from_path (graph, path, initial_node='i2'))
#' 
#' @export
get_distribution_from_path <- function (graph, path, initial_distribution=NULL, initial_node=NULL)
{
    path <- as.character (path)
    if (!all (path %in% graph$parts)) { stop ("'path' contains unknown part") }
    
    if (!is.null (initial_node) && !is.null (initial_distribution)) { stop ("'initial node' and 'initial distribution' cannot be both specified at the same time") }

    ## INITIALISE DISTRIBUTION
    length <- length (graph$nodes[[path[1]]])
    
    if (is.null (initial_node) && is.null (initial_distribution))
    {
        distribution <- rep (1 / length, length)
        names (distribution) <- graph$nodes[[path[1]]]
    }

    if (!is.null (initial_node))
    {
        if (!(initial_node %in% graph$nodes[[path[1]]])) { stop ("'initial node' was not found") }
        
        distribution <- rep (0, length)
        names (distribution) <- graph$nodes[[path[1]]]
        distribution[initial_node] <- 1
    }
    
    if (!is.null (initial_distribution))
    {
        if (length (initial_distribution) != length) { stop ("'initial distribution' has not the proper length") }
        if (abs (sum (initial_distribution) - 1) > triversity.env$distribution_precision) { stop ("'initial distribution' does not sum to 1") }

        distribution <- initial_distribution
        names (distribution) <- graph$nodes[[path[1]]]
    }

    ## COMPUTE DISTRIBUTION
    if (length (path) > 1)
    {
        compute_transition_from_path (graph, path)
        if (is.null (graph$edges$Climb(path)$matrix)) { return (NULL) }
        distribution <- as.simple_triplet_matrix (crossprod_simple_triplet_matrix (distribution, graph$edges$Climb(path)$matrix))
    }
    
    return (distribution)
}


#' @title Compute the diversity of a probability distribution.
#'
#' @description
#' \code{get_diversity_from_distribution} computes diversity values associated to an input
#' probability \code{distribution}. The implemented diversity measures all belong to the
#' parametrized family of "True Diversity" measures. They can either be specified by their
#' diversity \code{order} in [\code{0},\code{Inf}[ or by their \code{measure} name when it
#' corresponds to classical instances such as the richness, the Shannon entropy, the
#' Herfindahl-Hirschman index, or the Berger-Parker index.
#' 
#' @param distribution A vector of floats in [\code{0},\code{1}] and summing to \code{1}
#' giving the probability distribution whose diversity is measured.
#' 
#' @param order A vector of positive floats (possibly including \code{Inf}) giving the orders
#' of the diversity measures to be computed. If neither \code{order} nor \code{measure} is
#' specified, a predefined list of 8 diversity measures is computed.
#' 
#' @param measure A vector of strings giving the names of the diversity measures to compute.
#' Possible values are \code{richness}, \code{entropy}, \code{herfindahl}, and
#' \code{bergerparker}.
#'
#' @return A vector of positive floats giving the diversity values of the input probability
#' distribution.
#'
#' @examples
#' distribution <- c (1/4, 1/2, 1/12, 2/12)
#'
#' get_diversity_from_distribution (distribution)
#' get_diversity_from_distribution (distribution, order=c(0,Inf), measure='entropy')
#' 
#' @export
get_diversity_from_distribution <- function (distribution, order=NULL, measure=NULL)
{
    if (triversity.env$verbose) { message ('Computing diversity') }
    time <- Sys.time()

    if (is.null (order) && is.null (measure))
    {
        order <- c (0, 1, 2, Inf)
        measure <- c ('richness', 'entropy','herfindahl','bergerparker')
    }
    
    if (!is.null (order) && any (order < 0)) { stop ("'order' should be positive") }
    if (!is.null (measure) && any (! measure %in% c ('richness', 'entropy','herfindahl','bergerparker'))) { stop ("'measure' unknown (possible measures are 'richness', 'entropy', 'herfindahl', and 'bergerparker')") }

    if (abs (sum (distribution) - 1) > triversity.env$distribution_precision) { stop ("'distribution' does not sum to 1") }

    distribution <- distribution [distribution != 0] # suppress null values in the distribution

    diversity <- c()

    for (o in order)
    {
        if (o == 1) { d <- 1 / prod (distribution ^ distribution) } # order tends towards 1
        else if (o == Inf) { d <- 1 / max (distribution) } # order tends towards Inf
        else { d <- sum (distribution ^ o) ^ (1/(1-o)) }
        
        diversity <- append(diversity,d)
    }
    
    for (m in measure)
    {
        if (m == 'richness') { d <- length (distribution) }
        else if (m == 'entropy') { d <- - sum (distribution * log2 (distribution)) }
        else if (m == 'herfindahl') { d <- sum (distribution ^ 2) }
        else if (m == 'bergerparker') { d <- max (distribution) }
        
        diversity <- append(diversity,d)
    }

    names <- c()
    if (!is.null(order)) { names <- paste ('order', order, sep='') }
    names <- append (names, measure)
    names (diversity) <- names

    if (triversity.env$verbose) { message ('-> done in ', round (as.numeric (Sys.time() - time, units="secs"), 3), ' sec') }

    return (diversity)
}





#' @title Compute the diversity of a transition matrix, that is the individual diversities of
#' its rows or their mean diversity.
#'
#' @description
#' \code{get_diversity_from_transition} either computes the individual diversity values
#' associated to the rows of the input \code{transition} matrix, or the geometric means of
#' these values, optionally weighting them according to the input \code{mean_distribution}.
#' This hence allows to compute conditional diversity values associated to a transition
#' matrix.
#' 
#' @param transition A matrix of floats in [\code{0},\code{1}], with all lines summing to
#' \code{1}, giving the transition matrix from which the individual diversity values are
#' computed.
#'
#' @param type Either \code{'individual'}, to separately compute all individual diversities,
#' or \code{'mean'}, to compute their geometric mean.
#' 
#' @param mean_distribution (optional, only when \code{type == 'mean'}) A vector of floats
#' in [\code{0},\code{1}] and summing to \code{1} giving the probability distribution that
#' is used to weight the diversity values when computing their geometric means. It should
#' hence contain as many values as there are rows in the input \code{transition}. If not
#' specified, this distribution is assumed uniform.
#' 
#' @param order A vector of positive floats (possibly including \code{Inf}) giving the orders
#' of the diversity measures to be computed. If neither \code{order} nor \code{measure} is
#' specified, a predefined list of 8 diversity measures is computed.
#' 
#' @param measure A vector of strings giving the names of the diversity measures to compute.
#' Possible values are \code{richness}, \code{entropy}, \code{herfindahl}, and
#' \code{bergerparker}.
#' 
#' @return A matrix (or a vector) of positive floats giving the individual diversity values
#' (or their geometric mean) of the row of the input \code{transition} matrix.
#'
#' @examples
#' transition <- matrix (c (1/3, 1/3, 1/3, 0.9, 0.1, 0), nrow=2, ncol=3, byrow=TRUE)
#'
#' get_diversity_from_transition (transition, type='individual', order=c(0,Inf), measure='entropy')
#' get_diversity_from_transition (transition, type='mean')
#' get_diversity_from_transition (transition, type='mean', mean_distribution=c(1/4,3/4))
#' 
#' @export
get_diversity_from_transition <- function (transition, type='individual', mean_distribution=NULL, order=NULL, measure=NULL)
{
    if (is.null (order) && is.null (measure))
    {
        order <- c (0, 1, 2, Inf)
        measure <- c ('richness', 'entropy','herfindahl','bergerparker')
    }
    
    if (!is.null (order) && any (order < 0)) { stop ("'order' should be positive") }
    if (!is.null (measure) && any (! measure %in% c ('richness', 'entropy','herfindahl','bergerparker'))) { stop ("'measure' unknown (possible measures are 'richness', 'entropy', 'herfindahl', and 'bergerparker')") }

    if (! type %in% c('individual', 'mean')) { stop ("'type' unknown") }

    if (!is.simple_triplet_matrix (transition)) { transition <- as.simple_triplet_matrix (transition) }

    if (is.null (dimnames (transition)))
    {
        dimnames(transition)[[1]] <- seq (1, dim(transition)[1])
        dimnames(transition)[[2]] <- seq (1, dim(transition)[2])
    }

    if (!is.null(mean_distribution)) {
        mean_distribution <- as.vector (mean_distribution)
        if (length (mean_distribution) != transition$nrow) { stop ("'mean distribution' has not the proper length") }
        if (abs (sum (mean_distribution) - 1) > triversity.env$distribution_precision) { stop ("'mean distribution' does not sum to 1") }
    }

    ## TODO: test transition

    if (triversity.env$verbose) { message ('Computing diversity') }
    time <- Sys.time()

    names <- c()
    if (!is.null(order)) { names <- paste ('order', order, sep='') }
    names <- append (names, measure)

    ## BUILD
    diversity <- matrix (0, nrow=transition$nrow, ncol=length(order)+length(measure), dimnames=list(transition$dimnames[[1]], names))
 
    ## COMPUTE
    for (k in 1:length(transition$i)) {
        i <- transition$i[k]
        v <- transition$v[k]

        j <- 1
        for (o in order)
        {
            if (o == 0) { diversity[i,j] <- diversity[i,j] + 1 } # order tends towards 1
            if (o == 1) { diversity[i,j] <- diversity[i,j] - v * log2(v) } # order tends towards 1
            else if (o == Inf) { diversity[i,j] <- max (diversity[i,j], v) } # order tends towards Inf
            else { diversity[i,j] <- diversity[i,j] + v ^ o }
            j <- j+1
        }
        
        for (m in measure)
        {
            if (m == 'richness') { diversity[i,j] <- diversity[i,j] + 1 }
            else if (m == 'entropy') { diversity[i,j] <- diversity[i,j] - v * log2(v) }
            else if (m == 'herfindahl') { diversity[i,j] <- diversity[i,j] + v^2 }
            else if (m == 'bergerparker') { diversity[i,j] <- max (diversity[i,j], v) }
            j <- j+1
        }
    }

    ## NORMALISE
    j <- 1
    for (o in order)
    {
        if (o == 1) { diversity[,j] <- 2 ^ diversity[,j] } # order tends towards 1
        else if (o == Inf) { diversity[,j] <- 1 / diversity[,j] } # order tends towards Inf
        else if (o != 0) { diversity[,j] <- diversity[,j] ^ (1/(1-o)) }
        j <- j+1
    }

    if (type == 'mean') {
        if (!is.null(mean_distribution)) {
            mean_diversity <- c()

            for (o in order) {
                md <- prod (diversity[,paste('order',o,sep='')] ^ mean_distribution)
                mean_diversity <- append (mean_diversity, md)
            }

            for (m in measure) {
                if (m == 'richness') { md <- prod (diversity[,m] ^ mean_distribution) }
                else if (m == 'entropy') { md <- sum (mean_distribution * diversity[,m]) }
                else if (m == 'herfindahl') { md <- prod (diversity[,m] ^ mean_distribution) }
                else if (m == 'bergerparker') { md <- prod (diversity[,m] ^ mean_distribution) }
                mean_diversity <- append (mean_diversity, md)
            }

        } else {
            weight <- 1/length(dimnames(transition)[[1]])
            mean_diversity <- c()

            for (o in order) {
                md <- prod (diversity[,paste('order',o,sep='')] ^ weight)
                mean_diversity <- append (mean_diversity, md)
            }

            for (m in measure) {
                if (m == 'richness') { md <- prod (diversity[,m] ^ weight) }
                else if (m == 'entropy') { md <- sum (weight * diversity[,m]) }
                else if (m == 'herfindahl') { md <- prod (diversity[,m] ^ weight) }
                else if (m == 'bergerparker') { md <- prod (diversity[,m] ^ weight) }
                mean_diversity <- append (mean_diversity, md)
            }
        }
        
        names (mean_diversity) <- dimnames(diversity)[[2]]
        diversity <- mean_diversity
    }

    if (triversity.env$verbose) { message ('-> done in ', round (as.numeric (Sys.time() - time, units="secs"), 3), ' sec') }

    return (diversity)
}




#' @title Compute the diversity associated to a random walk following a path between the
#' parts of a multipartite graph.
#' 
#' @description
#' \code{get_diversity_from_path} computes diversity values of the probability distribution
#' of a random walk following a \code{path} between the different parts of the input
#' multipartite \code{graph}. It starts at a given part with an initial probability
#' distribution, then randomly follows the links of the graph between the different parts
#' according to the input \code{path}, then stops at the last specified part. The
#' implemented diversity measures all belong to the parametrized family of "True Diversity"
#' measures. They can either be specified by their diversity \code{order} in
#' [\code{0},\code{Inf}[ or by their \code{measure} name when it corresponds to classical
#' instances such as the richness, the Shannon entropy, the Herfindahl-Hirschman index,
#' or the Berger-Parker index.
#'
#' @param graph A multipartite graph obtained by calling the \code{\link{get_multipartite}}
#' function.
#' 
#' @param path A vector of character strings giving the path that the random walk should
#' follow between the different parts of the input multipartite \code{graph}. This path can
#' be as long as wanted, with eventual cycles, and each string it contains should refer to a
#' label in \code{graph$parts}.
#'
#' @param type Either \code{'individual'}, to separately compute all individual diversities,
#' \code{'mean'}, to compute their geometric mean, or \code{'collective'}, to compute the
#' overall diversity.
#' 
#' @param mean_distribution (optional, only when \code{type == 'mean'}) A vector of floats
#' in [\code{0},\code{1}] and summing to \code{1} giving the probability distribution that
#' is used to weight the diversity values when computing their geometric means. It should
#' hence contain as many values as there are rows in the input \code{transition}. If not
#' specified, this distribution is assumed uniform.
#' 
#' @param initial_distribution (optional, only when \code{type == 'collective'})
#' A vector of floats in [\code{0},\code{1}] and
#' summing to \code{1} giving the probability distribution to start with at the first part
#' of the input \code{path}. It should hence contain as many values as there are nodes in the
#' corresponding part. If not specified, this distribution is assumed uniform.
#'
#' @param initial_node (optional, only when \code{type == 'collective'})
#' A character string giving the label of a node in the first
#' part of the input \code{path}. This node is then considered to have probability one, thus
#' being equivalent to specifying an \code{initial_distribution} with only zeros except for
#' one node. If not specified, no such node is defined and the initial distribution is then
#' assumed uniform.
#'
#' @param order A vector of positive floats (possibly including \code{Inf}) giving the orders
#' of the diversity measures to be computed. If neither \code{order} nor \code{measure} is
#' specified, a predefined list of 8 diversity measures is computed.
#' 
#' @param measure A vector of strings giving the names of the diversity measures to compute.
#' Possible values are \code{richness}, \code{entropy}, \code{herfindahl}, and
#' \code{bergerparker}.
#' 
#' @return A matrix (or a vector) of positive floats giving the individual diversity values 
#' of the random walks following the input \code{path} (or their mean, or the collective
#' diversity).
#'
#' @examples
#' data (tripartite_example)
#' graph <- get_multipartite (data=tripartite_example)
#' path <- c(1,2,3)
#' 
#' 
#' get_diversity_from_path (graph, path, 'individual', measure=c('entropy','herfindahl'))
#' get_diversity_from_path (graph, path, 'individual', order=c(0,1,Inf))
#'
#' # Mean of individual diversities
#' get_diversity_from_path (graph, path, 'mean')
#' get_diversity_from_path (graph, path, 'mean', mean_distribution=c(1/3,2/3))
#'
#' # Collective diversities
#' get_diversity_from_path (graph, path, 'collective')
#' get_diversity_from_path (graph, path, 'collective', initial_distribution=c(0.75,0.25))
#' 
#' @export
get_diversity_from_path <- function (graph,
                                     path,
                                     type='individual',
                                     mean_distribution=NULL,
                                     initial_distribution=NULL,
                                     initial_node=NULL,
                                     order=NULL,
                                     measure=NULL)
{
    path <- as.character (path)
    if (!all (path %in% graph$parts)) { stop ("'path' contains unknown part") }

    if (! type %in% c('individual', 'mean', 'collective')) { stop ("'type' unknown") }
    if (type == 'individual' && ! is.null (initial_distribution)) { stop ("'initial distribution' cannot be used when computing individual diversities") }
    if (type == 'individual' && ! is.null (mean_distribution)) { stop ("'mean distribution' cannot be used when computing individual diversities") }
    if (type == 'mean' && ! is.null (initial_distribution)) { stop ("'initial distribution' cannot be used when computing mean of individual diversities") }
    if (type == 'mean' && ! is.null (initial_node)) { stop ("'initial node' cannot be used when computing mean of individual diversities") }
    if (type == 'collective' && ! is.null (initial_node)) { stop ("'initial node' cannot be used when computing collective diversities") }
    if (type == 'collective' && ! is.null (mean_distribution)) { stop ("'mean distribution' cannot be used when computing collective diversities") }

    if (type == 'individual' && !is.null (initial_node)) { type <- 'collective' }
    
    diversity <- NULL
    if (type == 'individual')
    {
        compute_transition_from_path (graph, path)    
        transition <- graph$edges$Climb(path)$matrix
        if (is.null (transition)) { return (NULL) }
        
        diversity <- get_diversity_from_transition (transition, type='individual', order=order, measure=measure)
    }
    
    else if (type == 'mean')
    {
        compute_transition_from_path (graph, path)    
        transition <- graph$edges$Climb(path)$matrix
        if (is.null (transition)) { return (NULL) }

        diversity <- get_diversity_from_transition (transition, type='mean', mean_distribution=mean_distribution, order=order, measure=measure)
    }

    else if (type == 'collective') {
        distribution <- get_distribution_from_path (graph, path, initial_distribution=initial_distribution, initial_node=initial_node)
        if (is.null (distribution)) { return (NULL) }
        
        diversity <- get_diversity_from_distribution (distribution, order=order, measure=measure)
    }

    return (diversity)
}



#' @title Compute all diversity flows associated to different random walks between the parts
#' of a multipartite graph.
#' 
#' @description
#' \code{get_all_diversities} builds on \code{\link{get_diversity_from_path}} to compute
#' the diversity values associated to all possible paths of a given \code{length} between
#' the different parts of the input multipartite \code{graph}.
#'
#' @param graph A multipartite graph obtained by calling the \code{\link{get_multipartite}}
#' function.
#' 
#' @param length A positive integer giving the maximal length of the computed paths.
#'
#' @param cycles A boolean value indicating if the diversity associated to paths
#' containing cycles should be computed.
#'
#' @param type A vector specifying the types of diversity to compute, among \code{'mean'} and
#' \code{'collective'}. See \code{\link{get_diversity_from_path}}.
#'
#' @param order A vector of positive floats (possibly including \code{Inf}) giving the orders
#' of the diversity measures to be computed. If neither \code{order} nor \code{measure} is
#' specified, a predefined list of 8 diversity measures is computed.
#' 
#' @param measure A vector of strings giving the names of the diversity measures to compute.
#' Possible values are \code{richness}, \code{entropy}, \code{herfindahl}, and
#' \code{bergerparker}.
#' 
#' @return A dataframe summarising all the computed diversities.
#'
#' @examples
#' data (tripartite_example)
#' graph <- get_multipartite (data=tripartite_example)
#' get_all_diversities (graph)
#' 
#' @export
get_all_diversities <- function (graph, length=2, cycles=TRUE, type=c('mean','collective'), order=NULL, measure=NULL)
{
    total_time <- Sys.time()

    diversity <- NULL
    path_labels <- c()
    type_labels <- c()
    
    for (l in 0:length)
    {
        if (triversity.env$verbose) { message ('Computing diversity of paths of length ', l) }
        time <- Sys.time()

        paths <- expand.grid (rep (list (graph$parts), l+1), stringsAsFactors=FALSE)
        for (row in 1:nrow(paths))
        {
            path <- rev (unlist (paths[row,]))
            if (cycles || (length(unique(path)) == length(path)))
            {
                path_label <- paste (path, collapse='-')

                for (t in type)
                {
                    if (l > 0 || t == 'collective') {
                        div <- suppressMessages (get_diversity_from_path (graph, path=path, type=t, order=order, measure=measure))

                        if (!is.null(div) && !is.na (div[1]))
                        {
                            path_labels <- append (path_labels, path_label)
                            type_labels <- append (type_labels, t)                    
                            if (is.null (diversity)) { diversity <- div } else { diversity <- rbind (diversity, div) }
                        }
                    }
                }
            }
        }
        
        if (triversity.env$verbose) { message ('-> done in ', round (as.numeric (Sys.time() - time, units="secs"), 3), ' sec') }
    }

    rownames (diversity) <- 1:nrow(diversity)
    diversity <- data.frame (path=path_labels, type=type_labels, diversity, stringsAsFactors=FALSE)

    if (triversity.env$verbose) { message ('-> total time ', round (as.numeric (Sys.time() - total_time, units="secs"), 3), ' sec') }

    return (diversity)
}



#' @title An example of dataframe encoding a small tripartite graph.
#'
#' @description
#' \code{tripartite_example} is a \code{data.frame} containing raw data that encodes a small
#' tripartite graph. It has the proper format to be loaded by \code{\link{get_multipartite}}.
#' It contains four columns. Each row gives a link between two nodes belonging to two
#' different parts of the tripartite graph. The first column gives the label if the part of
#' the first node and the second column gives the label of this node (any character string).
#' Similarly, the third and fourth columns give the labels of the second part and of the
#' second node. A fifth column could eventually be added to give the weights of the links
#' (any positive integer or float value).
#'
#' @docType data
#'
#' @usage data(tripartite_example)
#'
#' @examples
#' data (tripartite_example)
#' head (tripartite_example)
#' 
#' graph <- get_multipartite (data=tripartite_example)
#' graph
"tripartite_example"
