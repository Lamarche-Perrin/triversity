
distribution_precision <- 1e-5

#' Load a tripartite graph
#'
#' @param filename The file containing the tripartite graph
#' @param data A data.frame containing the tripartite graph
#' @return A tripartite graph
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



#' Get a probability transition matrix from a path in the tripartite graph
#'
#' @param path The path to follow to compute the transition matrix
#' @return The resulting transition matrix
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


#' Get a probability distribution from a path in the tripartite graph
#'
#' @param path The path to follow to compute the distribution
#' @return The resulting probability distribution
#'
#' @export
get_distribution_from_path <- function (tripartite, path, initial_node=NULL, initial_distribution=NULL)
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


#' Get the diversity associated to a probability distribution
#'
#' @param distribution The probability distribution to measure
#' @return The diversity of the input distribution
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


#' Get the conditional diversity associated to a probability distribution and a transition matrix
#'
#' @param distribution The probability distribution to measure
#' @param distribution The probability transition matrix to use
#' @return The diversity of the input distribution
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


#' Get the diversity associated to the distribution obtained through a path in the loaded tripartite graph
#'
#' @param path The path to follow to build the probability distribution
#' @return The diversity of the resulting probability distribution
#'
#' @export
get_diversity_from_path <- function (tripartite,
                                     path,
                                     conditional_path=NULL,
                                     initial_node=NULL,
                                     initial_distribution=NULL,
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
