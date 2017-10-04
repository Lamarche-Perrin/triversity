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


#' Compute diversity measures on tripartite graphs.
#'
#' \code{triversity} is an R package for the computation of diversity measures
#' on tripartite graphs. It first implements a parametrised family of such
#' diversity measures applying on probability distributions. Sometimes called
#' "True Diversity", this family contains famous measures such as the Richness,
#' the Shannon entropy, the Herfindahl-Hirschman index, and the Berger-Parker
#' index. This package then allows to apply such measures on probability
#' distributions resulting from random walks on tripartite graphs. By defining
#' an initial distribution at a given level in the graph, a path within the three
#' levels, and eventually a conditional step, the probability of the walker's
#' position within the final level is then computed, thus providing a particular
#' instance of diversity index.
#'
#' This package has been developed by researchers of the Complex Networks team,
#' (http://www.complexnetworks.fr/) within the Computer Science Laboratory of Paris 6
#' (https://www.lip6.fr/), for the AlgoDiv project (http://algodiv.huma-num.fr/),
#' founded by the French National Agency of Research
#' (http://www.agence-nationale-recherche.fr/) under grant ANR-15-CE38-0001.
#'
#' Contact: Robin Lamarche-Perrin <Robin.Lamarche-Perrin@lip6.fr>
#' 
#' See also my webpage: https://www-complexnetworks.lip6.fr/~lamarche/
#' 
#' List of main collaborators:
#' \itemize{
#' \item Lionel Tabourier
#' \item Fabien Tarissan
#' \item Rapha\"el Fournier S'niehotta
#' \item R\'emy Cazabet
#' }
#' 
#' Copyright \copyright 2017 Robin Lamarche-Perrin
#' 
#' \code{triversity} is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
#' 
#' 
#' @docType package
#' @name triversity
NULL

distribution_precision <- 1e-5

#' @title Get a properly-structured tripartite graph from raw data.
#' 
#' @description
#' \code{get_tripartite} returns a properly-structured tripartite graph
#' from a file or from a dataframe. The structure of the input data and
#' of the resultung data structure is detailed below.
#'
#' @param filename The path to the file containing raw data to build the tripartite graph.
#'
#' This input file should have
#' at least four columns, separated by spaces. Each row describe a link between
#' two nodes belonging to two different levels of the tripartite graph.
#' The first column indicates the level of the first node (any integer among
#' \code{1}, \code{2}, or \code{3}) and the
#' second column indicates its name (any character string). Similarly, the third
#' and fourth columns indicate the level and the name of the second node.
#' A fifth column can be added to indicate the eventual weights of the links
#' (any integer or float value).
#' 
#' @param data A \code{data.frame} containing the raw data to build the tripartite graph.
#'
#' This \code{data.frame} should have the same structure than the one described above
#' when using an input file: four columns indicating (in the same order) the
#' levels and the names of the two nodes constituting the link, and an optional
#' fifth column fot its weight.
#' 
#' @return A properly-structured tripartite graph that can be used by the other functions
#' of the \code{triversity} package.
#' 
#' The resulting data structure describe the different levels of the tripartite
#' graph, as well as its transition probabilities (encoded as sparse float matrices)
#' each following a given paths from one level to another. These transition matrices
#' are then used by functions such as \code{\link{get_distribution_from_path}} to
#' compute a distribution from a given path traveling through the different levels
#' of the graph.
#'
#' Once the \code{value} returned by \code{get_tripartite}:
#' \itemize{
#' \item \code{value$nodes} is a list of vectors contraining the names of the nodes
#' constituting the three levels of the tripartite graph
#' (resp. \code{value$nodes$level1}, \code{value$nodes$level2}, and \code{value$nodes$level3}).
#' \item \code{value$transitions} is a \code{data.tree} which nodes each contains a
#' transition matrix. For example, \code{value$transitions$level1$level2$mat} is the
#' transition matrix from level 1 to level 2.
#' }
#' 
#' @export
get_tripartite <- function (filename=NULL, data=NULL)
{
    if (is.null (filename) && is.null (data)) { stop ("'filename' or 'data' has to be specified") }
    if (!is.null (filename) && !is.null (data)) { stop ("'filename' and 'data' cannot be both specified at the same time") }

    tripartite <- list()

    ## LOAD DATA
    if (!is.null (filename)) {
        data <- read.table (filename, stringsAsFactors=FALSE)
    }

    if (length(colnames(data)) == 4) { colnames (data) <- c ('level1', 'id1', 'level2', 'id2') }
    if (length(colnames(data)) == 5) { colnames (data) <- c ('level1', 'id1', 'level2', 'id2', 'weight') }

    ## EXTRACT NODES
    tripartite$nodes <- lapply (c (level1=1, level2=2, level3=3), function (level) { return (unique (append (data[data$level1 == level, 'id1'], data[data$level2 == level, 'id2']))) })
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
    level21$mat <- t(level12$mat)
    level32$mat <- t(level23$mat)

    level12$mat <- level12$mat / rowSums(level12$mat)
    level21$mat <- level21$mat / rowSums(level21$mat)
    level23$mat <- level23$mat / rowSums(level23$mat)
    level32$mat <- level32$mat / rowSums(level32$mat)

    return (tripartite)
}



#' @title Get the probability transition matrix corresponding to a path
#' in a tripartite graph.
#'
#' @param tripartite A tripartite graph obtained by \code{\link{get_tripartite}}.
#' @param path The path to follow in the graph to compute the transition matrix.
#' @return The resulting transition matrix, that is a matrix of floats which rows sum to one.
#'
#' @details Note that the tripartite graph structure implemented in this package
#' store in memory any computed transition matrix to avoid redundant computation
#' in the future. Hence, the first execution of \code{get_transition_from_path} (or
#' any other function that is built on it), can be much slower than latter calls.
#' The transition matrices are stored in a \code{data.tree} (see \code{tripartite$transitions}).
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


#' @title Get the probability distribution corresponding to a path in the tripartite graph
#'
#' @param tripartite A tripartite graph obtained by the \code{\link{get_tripartite}}.
#' @param path The path to follow in the graph to compute the distribution.
#' @param initial_distribution A vector of floats describing the initial distribution to
#' start with at the first level of the path. It should sum to 1 and contains as many values
#' as there are nodes in the graph's level specified by the first step of the path.
#' @param initial_node The name of a node in the first level of the path.
#' @return The resulting probability distribution, that is a vector of floats which sum to 1.
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


#' @title Get the diversity value associated to a probability distribution.
#'
#' @param distribution A vector of floats describing the probability distribution to measure.
#' @param order A vector of positive floats (possibly including \code{Inf}) describing the
#' orders of the diversity measures to be computed.
#' @param measure A vector of strings giving the names of the diversity measures to compute.
#' The following possible values are \code{'richness'}, \code{'entropy'}, \code{'herfindahl'},
#' and \code{'bergerparker'}.
#' @return A vector containing the diversity values of the input distribution.
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


#' @title Get the conditional diversity values associated to a probability distribution and a transition matrix.
#'
#' @param distribution The probability distribution to measure.
#' @param transition The probability transition matrix to use.
#' @param order A vector of positive floats (possibly including \code{Inf}) describing the
#' orders of the diversity measures to be computed.
#' @param measure A vector of strings giving the names of the diversity measures to compute.
#' The following possible values are \code{'richness'}, \code{'entropy'}, \code{'herfindahl'},
#' and \code{'bergerparker'}.
#' @return The diversity of the input distribution.
#'
#' @export
get_conditional_diversity_from_distribution <- function (distribution, transition, order=NULL, measure=NULL)
{
    if (is.null (order) && is.null (measure))
    {
        order <- c (0, 1, 2, Inf)
        measure <- c ('richness', 'entropy','herfindahl','bergerparker')
    }
    
    if (!is.null (order) && any (order < 0)) { stop ("'order' should be positive") }
    if (!is.null (measure) && any (! measure %in% c ('richness', 'entropy','herfindahl','bergerparker'))) { stop ("'measure' unknown (possible measures are 'richness', 'entropy', 'herfindahl', and 'bergerparker')") }
    
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


#' @title Get diversity values associated to the distribution obtained through a path in a tripartite graph.
#'
#' @param tripartite A tripartite graph obtained by the \code{\link{get_tripartite}}.
#' @param path The path to follow to build the probability distribution
#' @param conditional_path A vector of integers (among \code{1}, \code{2}, \code{3}) giving
#' an eventual path to follow before computing and aggregating the diversity measures.
#' @param initial_distribution A vector of floats describing the initial distribution to
#' start with at the first level of the path. It should sum to 1 and contains as many values
#' as there are nodes in the graph's level specified by the first step of the path.
#' @param initial_node The name of a node in the first level of the path.
#' @param order A vector of positive floats (possibly including \code{Inf}) describing the
#' orders of the diversity measures to be computed.
#' @param measure A vector of strings giving the names of the diversity measures to compute.
#' The following possible values are \code{'richness'}, \code{'entropy'}, \code{'herfindahl'},
#' and \code{'bergerparker'}.
#' @return The diversity of the resulting probability distribution
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
        distribution <- get_distribution_from_path (tripartite, conditional_path, initial_node, initial_distribution)
        transition <- get_transition_from_path (tripartite, path)
        diversity <- get_conditional_diversity_from_distribution (distribution, transition, order, measure)        
    } else {
        distribution <- get_distribution_from_path (tripartite, path, initial_node, initial_distribution)
        diversity <- get_diversity_from_distribution (distribution, order, measure)
    }
    
    return (diversity)
}
