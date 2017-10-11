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


#' @title Compute diversity measures on tripartite graphs.
#'
#' @description
#' \code{triversity} is an R package for the computation of diversity measures
#' on tripartite graphs. First, it implements a parametrized family of such diversity
#' measures which apply on probability distributions. Sometimes called "True Diversity",
#' this family contains famous measures such as the richness, the Shannon entropy, the
#' Herfindahl-Hirschman index, and the Berger-Parker index. Second, the package allows
#' to apply these measures on probability distributions resulting from random walks between
#' the levels of tripartite graphs. By defining an initial distribution at a given level of
#' the graph and a path to follow between the three levels, the probability of the walker's
#' position within the final level is then computed, thus providing a particular instance of
#' diversity to measure.
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
#' \code{triversity} is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see \url{http://www.gnu.org/licenses/}.
#' 
#' 
#' @docType package
#' 
#' @name triversity
#' 
#' @import Matrix data.tree
NULL


distribution_precision <- 1e-5

#' @title Build a properly-structured tripartite graph from raw data.
#' 
#' @description
#' \code{get_tripartite} builds a properly-structured tripartite graph
#' from a file or from a \code{data.frame} containing raw data. This object can then
#' be used by the other functions of this package. The structure
#' of the input data and of the resulting structure is detailed below.
#'
#' @param filename A string giving the path to the file containing the raw data.
#'
#' The input file should contain
#' at least four columns, separated by spaces. Each row gives a link between
#' two nodes belonging to two different levels of the tripartite graph.
#' The first column gives the level of the first node (any integer in
#' \{\code{1}, \code{2}, \code{3}\}) and the
#' second column gives its name (any character string). Similarly, the third
#' and fourth columns give the level and the name of the second node.
#' Note that no link is allowed between level \code{1} and level \code{3}.
#' The file can eventually contain a fifth column giving the weights of the links
#' (any positive integer or float value).
#' 
#' @param data A \code{data.frame} containing the raw data.
#'
#' This \code{data.frame} should have the same structure than the one described above
#' for the case of an input file: four columns indicating the
#' levels and the names of the two nodes constituting the link, and an optional
#' fifth column for its weight. At least \code{filename} or \code{data} should be
#' specified, but not both at the same time.
#' 
#' @return A properly-structured tripartite graph that can be used by the other functions
#' of the \code{triversity} package.
#' 
#' The resulting object encodes the names of the nodes at the three levels of the tripartite
#' graph, as well as the transition matrices of random walks following different paths
#' between levels (encoded as sparse matrices of floats in [\code{0},\code{1}], with all
#' rows summing to \code{1}). These transition matrices
#' are then used by functions such as \code{\link{get_distribution_from_path}} to
#' compute the probability distributions of such random walks, or such as
#' \code{\link{get_diversity_from_path}} to compute the diversity of these distributions.
#'
#' Assuming the object returned by \code{get_tripartite} is stored in the \code{tripartite}
#' variable, then:
#' \itemize{
#' \item \code{tripartite$nodes} is a list of string vectors giving the names of the nodes
#' constituting the three levels of the tripartite graph
#' (resp. \code{tripartite$nodes$level1}, \code{tripartite$nodes$level2}, and
#' \code{tripartite$nodes$level3}).
#' \item \code{tripartite$transitions} is a \code{data.tree} whose nodes each contains the
#' transition matrix of the corresponding random walk. For example,
#' \code{tripartite$transitions$level1$level2$mat} gives the
#' transition matrix from level 1 to level 2.
#' }
#'
#' @examples
#' data (tripartite_example)
#' tripartite <- get_tripartite (data=tripartite_example)
#'
#' tripartite$nodes$level1
#' tripartite$nodes$level2
#' tripartite$level1$level2$mat
#' tripartite
#' 
#' @export
get_tripartite <- function (filename=NULL, data=NULL)
{
    if (is.null (filename) && is.null (data)) { stop ("'filename' or 'data' has to be specified") }
    if (!is.null (filename) && !is.null (data)) { stop ("'filename' and 'data' cannot be both specified at the same time") }

    ## LOAD DATA
    if (!is.null (filename)) { data <- utils::read.table (filename, stringsAsFactors=FALSE) }

    if (length(colnames(data)) < 4) { stop ("input data should have at least 4 columns") }
    
    colnames(data)[1:4] <- c ('level1', 'id1', 'level2', 'id2')
    if (length(colnames(data)) > 4) { colnames(data)[5] <- 'weight' }

    ## EXTRACT NODES
    tripartite <- list()
    tripartite$nodes <- lapply (c (level1=1, level2=2, level3=3), function (level) { return (unique (append (as.character (data[data$level1 == level, 'id1']), as.character (data[data$level2 == level, 'id2'])))) })
    tripartite$nodes <- lapply (tripartite$nodes, function (labels) { return (labels [order (labels)]) })

    ## BUILD ONE-STEP MATRICES
    tripartite$transitions <- Node$new ('tree')
    level1 <- tripartite$transitions$AddChild ('level1')
    level2 <- tripartite$transitions$AddChild ('level2')
    level3 <- tripartite$transitions$AddChild ('level3')

    level12 <- level1$AddChild ('level2')
    level21 <- level2$AddChild ('level1')
    level23 <- level2$AddChild ('level3')
    level32 <- level3$AddChild ('level2')

    data12 <- data[data$level1 == 1 & data$level2 == 2,]
    data21 <- data[data$level1 == 2 & data$level2 == 1,]
    data21[,c(1,2,3,4)] <- data21[,c(3,4,1,2)]
    data12 <- rbind (data12, data21)

    i <- match (data12[,'id1'], tripartite$nodes$level1)
    j <- match (data12[,'id2'], tripartite$nodes$level2)
    level12$mat <- sparseMatrix (i, j, x=1, dims=c (length (tripartite$nodes$level1), length (tripartite$nodes$level2)), dimnames=list (tripartite$nodes$level1, tripartite$nodes$level2))

    data23 <- data[data$level1 == 2 & data$level2 == 3,]
    data32 <- data[data$level1 == 3 & data$level2 == 2,]
    data32[,c(1,2,3,4)] <- data32[,c(3,4,1,2)]
    data23 <- rbind (data23, data32)

    i <- match (data23[,'id1'], tripartite$nodes$level2)
    j <- match (data23[,'id2'], tripartite$nodes$level3)
    level23$mat <- sparseMatrix (i, j, x=1, dims=c (length (tripartite$nodes$level2), length (tripartite$nodes$level3)), dimnames=list (tripartite$nodes$level2, tripartite$nodes$level3))

    ## COMPLETE AND NORMALISE ONE-STEP MATRICES
    level21$mat <- Matrix::t(level12$mat)
    level32$mat <- Matrix::t(level23$mat)

    level12$mat <- level12$mat / Matrix::rowSums(level12$mat)
    level21$mat <- level21$mat / Matrix::rowSums(level21$mat)
    level23$mat <- level23$mat / Matrix::rowSums(level23$mat)
    level32$mat <- level32$mat / Matrix::rowSums(level32$mat)

    return (tripartite)
}



#' @title Compute the transition matrix of a random walk following a path between the levels of a
#' tripartite graph.
#'
#' @description
#' \code{get_transition_from_path} computes the transition matrix of a random walk following a
#' \code{path} between the different levels of the input \code{tripartite} graph.
#'
#' @param tripartite A tripartite graph obtained by calling the \code{\link{get_tripartite}} function.
#'
#' @param path A vector of integers in \{\code{1}, \code{2}, \code{3}\} giving the path of the random
#' walk between the different levels of the input \code{tripartite} graph. This path can be as long
#' as wanted. Two successive levels should always be adjacent, that is the random walk cannot go from
#' level \code{1} to level \code{3} (or conversely) without first going through level \code{2}.
#'
#' @return A matrix of floats in [\code{0},\code{1}], with all lines summing to \code{1},
#' giving the transition matrix of the random walk following the input \code{path}.
#'
#' @details
#' Note that the tripartite graph structure implemented in this package
#' stores in memory any computed transition matrix to avoid redundant computation
#' in the future. Hence, the first execution of \code{get_transition_from_path}, or of
#' any other function that builds on it, can be much slower than latter calls.
#' The transition matrices are stored within a \code{data.tree} in the input \code{tripartite}
#' variable (see \code{tripartite$transitions}).
#'
#' @examples
#' data (tripartite_example)
#' tripartite <- get_tripartite (data=tripartite_example)
#' 
#' get_transition_from_path (tripartite, c(2,1,2,3))
#' 
#' @export
get_transition_from_path <- function (tripartite, path)
{
    strpath <- paste ('level', path, sep='')
    
    current_step <- 1
    current_node <- tripartite$transitions$Climb(strpath[current_step])

    while (current_step < length (path))
    {
        past_step <- current_step
        past_node <- current_node

        current_step <- current_step + 1
        current_node <- current_node$Climb(strpath[current_step])
        
        if (is.null(current_node))
        {
            current_node <- past_node$AddChild(strpath[current_step])
            current_node$mat <- past_node$mat %*% tripartite$transitions$Climb(strpath[c(past_step,current_step)])$mat
        }
    }
    
    #current_node$mat <- current_node$mat / rowSums(current_node$mat)

    return (current_node$mat)
}


#' @title Compute the probability distribution associated to a random walk following a path
#' between the levels of a tripartite graph.
#'
#' @description
#' \code{get_distribution_from_path} computes the probability distribution of a random walk within
#' the different levels of the input \code{tripartite} graph. It starts at a given level with an
#' initial probability distribution, then randomly follows the links of the graph between the
#' different levels according to the input \code{path}, then stops at the last specified level.
#' 
#' @param tripartite A tripartite graph obtained by calling the \code{\link{get_tripartite}} function.
#'
#' @param path A vector of integers in \{\code{1}, \code{2}, \code{3}\} giving the path of the random
#' walk between the different levels of the input \code{tripartite} graph. This path can be as long
#' as wanted. Two successive levels should always be adjacent, that is the random walk cannot go from
#' level \code{1} to level \code{3} (or conversely) without first going through level \code{2}.
#' 
#' @param initial_distribution A vector of floats in [\code{0},\code{1}] and summing to
#' \code{1} giving the probability
#' distribution to start with at the first level of the input \code{path}. It should hence
#' contain as many values as there are nodes in the corresponding level. If not specified, this
#' distribution is assumed uniform.
#'
#' @param initial_node A string giving the name of a node in the first level of the input
#' \code{path}. This node is then considered to have probability one, thus being equivalent to
#' specifying an \code{initial_distribution} with only zeros except for one node. If not specified,
#' no such node is defined and the initial distribution is assumed uniform.
#'
#' @return A vector of floats in [\code{0},\code{1}] and summing to \code{1} giving the probability distribution of the
#' random walk when arriving at the last level, after having followed the input \code{path} within
#' the different levels of the graph.
#'
#' @examples
#' data (tripartite_example)
#' tripartite <- get_tripartite (data=tripartite_example)
#'
#' get_distribution_from_path (tripartite, path=c(2,1,2,3))
#' get_distribution_from_path (tripartite, path=c(2,1,2,3), initial_distribution=c(0.25,0,0.25,0.5))
#' get_distribution_from_path (tripartite, path=c(2,1,2,3), initial_node='i2')
#' 
#' @export
get_distribution_from_path <- function (tripartite, path, initial_distribution=NULL, initial_node=NULL)
{
    strpath <- paste ('level', path, sep='')

    ## INITIALISE DISTRIBUTION
    if (!is.null (initial_node) && !is.null (initial_distribution)) { stop ("'initial node' and 'initial distribution' cannot be both specified at the same time") }

    current_step <- 1
    current_length <- length (tripartite$nodes[[strpath[current_step]]])
        
    if (is.null (initial_node) && is.null (initial_distribution))
    {
        distribution <- rep (1 / current_length, current_length)
        names (distribution) <- tripartite$nodes[[strpath[current_step]]]
    }

    if (!is.null (initial_node))
    {
        if (!(initial_node %in% tripartite$nodes[[strpath[current_step]]])) { stop ("'initial node' was not found") }
        
        distribution <- rep (0, current_length)
        names (distribution) <- tripartite$nodes[[strpath[current_step]]]
        distribution[initial_node] <- 1
    }
    
    if (!is.null (initial_distribution))
    {
        if (length (initial_distribution) != current_length) { stop ("'initial distribution' has not the proper length") }
        if (abs (sum (initial_distribution) - 1) > distribution_precision) { stop ("'initial distribution' does not sum to 1") }

        distribution <- initial_distribution
        names (distribution) <- tripartite$nodes[[strpath[current_step]]]
    }

    ## RETURN RESULTING DISTRIBUTION
    if (length (path) > 1)
    {
        distribution <- distribution %*% get_transition_from_path (tripartite, path)
        distribution <- distribution[1,]
    }
    
    #distribution <- distribution / sum (distribution)
    
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
#' @param distribution A vector of floats in [\code{0},\code{1}] and summing to \code{1} giving
#' the probability distribution whose diversity is measured.
#' 
#' @param order A vector of positive floats (possibly including \code{Inf}) giving the
#' orders of the diversity measures to be computed. If neither \code{order} nor \code{measure} is
#' specified, a predefined list of 8 diversity measures is computed.
#'
#' @param measure A vector of strings giving the names of the diversity measures to compute.
#' Possible values are \code{richness}, \code{entropy}, \code{herfindahl},
#' and \code{bergerparker}.
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
    if (is.null (order) && is.null (measure))
    {
        order <- c (0, 1, 2, Inf)
        measure <- c ('richness', 'entropy','herfindahl','bergerparker')
    }
    
    if (!is.null (order) && any (order < 0)) { stop ("'order' should be positive") }
    if (!is.null (measure) && any (! measure %in% c ('richness', 'entropy','herfindahl','bergerparker'))) { stop ("'measure' unknown (possible measures are 'richness', 'entropy', 'herfindahl', and 'bergerparker')") }

    if (abs (sum (distribution) - 1) > distribution_precision) { stop ("'distribution' does not sum to 1") }

    distribution <- distribution [distribution != 0] # suppress null values in the distribution

    d <- 1
    diversity <- c()

    for (o in order)
    {
        if (o == 1) { diversity[d] <- 1 / prod (distribution ^ distribution) } # order tends towards 1
        else if (o == Inf) { diversity[d] <- 1 / max (distribution) } # order tends towards Inf
        else { diversity[d] <- sum (distribution ^ o) ^ (1/(1-o)) }
        
        names(diversity)[d] <- paste ('order', o)
        d <- d+1
    }
    
    for (m in measure)
    {
        if (m == 'richness') { diversity[d] <- get_diversity_from_distribution (distribution, 0) }
        else if (m == 'entropy') { diversity[d] <- log2 (get_diversity_from_distribution (distribution, 1)) }
        else if (m == 'herfindahl') { diversity[d] <- 1 / get_diversity_from_distribution (distribution, 2) }
        else if (m == 'bergerparker') { diversity[d] <- 1 / get_diversity_from_distribution (distribution, Inf) }
        
        names(diversity)[d] <- m
        d <- d+1
    }

    return (diversity)
}


#' @title Compute the conditional diversity of a transition matrix.
#'
#' @description
#' \code{get_conditional_diversity_from_transition} computes the geometric means of diversity
#' values associated to the lines of the input \code{transition} matrix, while weighting these values
#' according to an optional input \code{distribution}. This hence allows to compute conditional
#' diversity values associated to the matrix.
#' 
#' @param transition A matrix of floats in [\code{0},\code{1}], with all lines summing to \code{1},
#' giving the transition matrix from which the conditional diversity values are computed.
#' 
#' @param distribution A vector of floats in [\code{0},\code{1}] and summing to \code{1} giving
#' the probability
#' distribution that is used to weight the diversity values when computing their geometric means.
#' It should hence contain as many values as there are rows in the input \code{transition}. If not
#' specified, this distribution is assumed uniform.
#' 
#' @param order A vector of positive floats (possibly including \code{Inf}) giving the
#' orders of the diversity measures to be computed. If neither \code{order} nor \code{measure} is
#' specified, a predefined list of 8 diversity measures is computed.
#' 
#' @param measure A vector of strings giving the names of the diversity measures to compute.
#' Possible values are \code{richness}, \code{entropy}, \code{herfindahl},
#' and \code{bergerparker}.
#' 
#' @return A vector of positive floats giving the conditional diversity values of the input
#' \code{transition} matrix, that is the geometric means of the diversity values associated
#' to its rows.
#'
#' @examples
#' transition <- matrix (c (1/3, 1/3, 1/3, 0.9, 0.1, 0), nrow=2, ncol=3, byrow=TRUE)
#' get_conditional_diversity_from_transition (transition, order=c(0,Inf), measure='entropy')
#' get_conditional_diversity_from_transition (transition, distribution=c(0.75,0.25))
#' 
#' @export
get_conditional_diversity_from_transition <- function (transition, distribution=NULL, order=NULL, measure=NULL)
{
    if (is.null (order) && is.null (measure))
    {
        order <- c (0, 1, 2, Inf)
        measure <- c ('richness', 'entropy','herfindahl','bergerparker')
    }
    
    if (!is.null (order) && any (order < 0)) { stop ("'order' should be positive") }
    if (!is.null (measure) && any (! measure %in% c ('richness', 'entropy','herfindahl','bergerparker'))) { stop ("'measure' unknown (possible measures are 'richness', 'entropy', 'herfindahl', and 'bergerparker')") }

    if (is.null (dimnames (transition)))
    {
        dimnames(transition)[[1]] <- seq (1, dim(transition)[1])
        dimnames(transition)[[2]] <- seq (1, dim(transition)[2])
    }    
    
    if (is.null (distribution)) { distribution <- rep (1 / dim(transition)[1], dim(transition)[1]) }

    if (is.null (names (distribution))) { names (distribution) <- dimnames (transition)[[1]] }
    
    if (abs (sum (distribution) - 1) > distribution_precision) { stop ("'distribution' does not sum to 1") }

    d <- 1
    conditional_diversity <- c()
    
    for (o in order)
    {
        conditional_diversity[d] <- 1
        
        for (i in names (distribution)) { conditional_diversity[d] <- conditional_diversity[d] * get_diversity_from_distribution (transition[i,], order=o) ^ distribution[i] }
        
        names(conditional_diversity)[d] <- paste ('order', o)
        d <- d+1
    }
    
    for (m in measure)
    {
        if (m == 'richness') { o <- 0 }
        else if (m == 'entropy') { o <- 1 }
        else if (m == 'herfindahl') { o <- 2 }
        else if (m == 'bergerparker') { o <- Inf }

        conditional_diversity[d] <- 1
    
        for (i in names (distribution)) { conditional_diversity[d] <- conditional_diversity[d] * get_diversity_from_distribution (transition[i,], order=o) ^ distribution[i] }

        if (m == 'entropy') { conditional_diversity[d] <- log2 (conditional_diversity[d]) }
        else if (m == 'herfindahl') { conditional_diversity[d] <- 1 / conditional_diversity[d] }
        else if (m == 'bergerparker') { conditional_diversity[d] <- 1 / conditional_diversity[d] }

        names(conditional_diversity)[d] <- m
        d <- d+1
    }
    
    return (conditional_diversity)
}


#' @title Compute the diversity associated to a random walk following a path between the levels
#' of a tripartite graph.
#' 
#' @description
#' \code{get_diversity_from_path} computes diversity values of the probability distribution of a
#' random walk following a \code{path} between the different levels of the input \code{tripartite}
#' graph. It starts at a given level with an initial probability distribution, then randomly follows
#' the links of the graph between the different levels according to the input \code{path}, then
#' stops at the last specified level. The implemented diversity measures all belong to the
#' parametrized family of "True Diversity" measures. They can either be specified by their diversity
#' \code{order} in [\code{0},\code{Inf}[ or by their \code{measure} name when it corresponds to
#' classical instances such as the richness, the Shannon entropy, the Herfindahl-Hirschman index,
#' or the Berger-Parker index.
#'
#' @param tripartite A tripartite graph obtained by calling the \code{\link{get_tripartite}} function.
#' 
#' @param path A vector of integers in \{\code{1}, \code{2}, \code{3}\} giving the path of the random
#' walk between the different levels of the input \code{tripartite} graph. This path can be as long
#' as wanted. Two successive levels should always be adjacent, that is the random walk cannot go from
#' level \code{1} to level \code{3} (or conversely) without first going through level \code{2}.
#'
#' @param conditional_path A vector of integers in \{\code{1}, \code{2}, \code{3}\} eventually
#' giving another path to compute conditional diversity values instead of regular diversity
#' values. When specified, this conditional
#' path is first used to initiate the random walk. The resulting probability distribution is then
#' used to weight the individual diversity values obtained on the input \code{path} when computing
#' their geometric means (see
#' \code{\link{get_conditional_diversity_from_transition}}). This path can be as long
#' as wanted. The last level of the conditional path should be the same as the first level of the
#' input \code{path}. Moreover, two successive levels should always be adjacent, that is the random
#' walk cannot go from
#' level \code{1} to level \code{3} (or conversely) without first going through level \code{2}.
#'
#' @param initial_distribution A vector of floats in [\code{0},\code{1}] and summing to
#' \code{1} giving the probability
#' distribution to start with at the first level of the input \code{path}, or at the first level
#' of the input \code{conditional_path} if specified. It should hence
#' contain as many values as there are nodes in the corresponding level. If not specified, this
#' distribution is assumed uniform.
#'
#' @param initial_node A string giving the name of a node in the first level of the input
#' \code{path}, or at the first level of the input \code{conditional_path} if specified.
#' This node is then considered to have probability one, thus being equivalent to
#' specifying an \code{initial_distribution} with only zeros except for one node. If not specified,
#' no such node is defined and the initial distribution is assumed uniform.
#'
#' @param order A vector of positive floats (possibly including \code{Inf}) giving the
#' orders of the diversity measures to be computed. If neither \code{order} nor \code{measure} is
#' specified, a predefined list of 8 diversity measures is computed.
#'
#' @param measure A vector of strings giving the names of the diversity measures to compute.
#' Possible values are \code{richness}, \code{entropy}, \code{herfindahl},
#' and \code{bergerparker}.
#' 
#' @return A vector of positive floats giving the diversity values (or conditional diversity values)
#' of the random walk following the input \code{path}. 
#'
#' @examples
#' data (tripartite_example)
#' tripartite <- get_tripartite (data=tripartite_example)
#'
#' 
#' # COMPUTING DIFFERENT DIVERSITY VALUES FOR A GIVEN PATH
#' 
#' # Herfindahl-Hirschman index of nodes in level 3 wrt nodes in level 1
#' get_diversity_from_path (tripartite, path=c(1,2,3), measure='herfindahl')
#' 1 / get_diversity_from_path (tripartite, path=c(1,2,3), order=2)
#'
#' # Shannon entropy of nodes in level 3 wrt nodes in level 1
#' get_diversity_from_path (tripartite, path=c(1,2,3), measure='entropy')
#' log2 (get_diversity_from_path (tripartite, path=c(1,2,3), order=1))
#' 
#' # Some other diversity values of nodes in level 3 wrt nodes in level 1
#' get_diversity_from_path (tripartite, path=c(1,2,3), order=c(1,2,Inf),
#'                          measure=c('richness','bergerparker'))
#' 
#' # Eight of the main diversity values of nodes in level 3 wrt nodes in level 1
#' get_diversity_from_path (tripartite, path=c(1,2,3))
#'
#' 
#' # SPECIFYING THE INITIAL DISTRIBUTION
#' 
#' # Diversity of nodes in level 3 wrt nodes in level 1 (with non-uniform weights)
#' get_diversity_from_path (tripartite, path=c(1,2,3), initial_distribution=c(0.75,0.25))
#' 
#' # Individual diversity of nodes in level 3 wrt node 'u1' in level 1
#' get_diversity_from_path (tripartite, path=c(1,2,3), initial_node='u1')
#' get_diversity_from_path (tripartite, path=c(1,2,3), initial_distribution=c(1,0))
#'
#'
#' # COMPUTING THE MEAN OF INDIVIDUAL DIVERSITES
#' 
#' # Mean of individual diversities of nodes in level 3 wrt nodes in level 2 (with
#' # uniform weights)
#' get_diversity_from_path (tripartite, path=c(2,3), conditional_path=c(2))
#'
#' # Mean of individual diversities of nodes in level 3 wrt nodes in level 2 (weighted
#' # according to the path from level 1 to level 2, with a uniform distribution in level 1)
#' get_diversity_from_path (tripartite, path=c(2,3), conditional_path=c(1,2))
#'
#' # Mean of individual diversities of nodes in level 3 wrt nodes in level 2 (weighted
#' # according to the path from level 1 to level 2, with only node 'u1' in level 1)
#' get_diversity_from_path (tripartite, path=c(2,3), conditional_path=c(1,2),
#'                          initial_node='u1')
#' 
#' @export
get_diversity_from_path <- function (tripartite,
                                     path,
                                     conditional_path=NULL,
                                     initial_distribution=NULL,
                                     initial_node=NULL,
                                     order=NULL,
                                     measure=NULL)
{
    if (!is.null (conditional_path))
    {
        distribution <- get_distribution_from_path (tripartite, conditional_path, initial_distribution, initial_node)
        transition <- get_transition_from_path (tripartite, path)
        diversity <- get_conditional_diversity_from_transition (transition, distribution, order, measure)        
    } else {
        distribution <- get_distribution_from_path (tripartite, path, initial_distribution, initial_node)
        diversity <- get_diversity_from_distribution (distribution, order, measure)
    }
    
    return (diversity)
}



#' @title An example of dataframe encoding a small tripartite graph.
#'
#' @description
#' \code{tripartite_example} is a \code{data.frame} containing raw data
#' that encodes a small tripartite graph. It has the proper
#' format to be loaded by \code{\link{get_tripartite}}. It contains four columns.
#' Each row gives a link between
#' two nodes belonging to two different levels of the tripartite graph.
#' The first column gives the level of the first node (any integer in
#' \{\code{1}, \code{2}, \code{3}\}) and the
#' second column gives its name (any character string). Similarly, the third
#' and fourth columns give the level and the name of the second node.
#' A fifth column could eventually be added to give the weights of the links
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
#' # Load the data.frame into a proper data structure
#' tripartite <- get_tripartite (data=tripartite_example)
#'
#' # Get the names of the nodes in the second level of the tripartite graph
#' tripartite$nodes$level2
#'
#' # Get the transition matrix of a random walk going from the level 2 to level 1
#' tripartite$transitions$level2$level1$mat
"tripartite_example"
