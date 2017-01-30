##' Visualise distributions of generation times
##'
##' @export
##'
##' @param store A gensig output dataframe
##' @param param A gensig param list
##'
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_line
##' @importFrom ggplot2 labs
##' @importFrom ggplot2 aes_string

vis.w <- function(store, param = create.param()) {

    names <- names(param)

    df <- data.frame(days = seq_along(param[[names[1]]]$w))

    for(i in names) df[[i]] <- param[[i]]$w

    mlt <- reshape2::melt(df, id="days")
    mlt$variable <- factor(mlt$variable, levels = sort.gensig(store))

    p <- ggplot(mlt, aes_string('days', 'value', colour = 'variable')) + geom_line(size = 1.5) +
        labs(title = "Generation Times",
             subtitle = "Generation times of various pathogens",
             x = "Days",
             y = "Density")

    return(p)

}




##' Visualise distributions of genetic signatures
##'
##' @export
##'
##' @aliases vis.dist
##'
##' @param store A gensig output dataframe
##' @param config A gensig config list
##'
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 labs
##' @importFrom ggplot2 aes_string
##' @importFrom ggplot2 geom_bar
##' @importFrom ggplot2 facet_wrap
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 as_labeller

vis.dist <- function(store, config = create.config()) {

    ## Calculate proportions manually
    df <- by(store$gensig, store$disease, function(sig) prop.table(table(sig)))

    ## Collapse the list of dataframes into a single dataframe
    df <- plyr::ldply(df, data.frame)

    ## Sort the pathogen names by mean gensig
    df$.id <- factor(df$.id, levels = sort.gensig(store))

    p <- ggplot(df, aes_string(x = 'sig', y = 'Freq', fill = '.id', colour = '.id')) +
        geom_bar(stat = "identity") +
        facet_wrap( ~ .id,labeller = as_labeller(config$label)) +
        theme(legend.position = "none") +
        labs(title = "Genetic signatures", x = "Number of mutations", y = "Proportion")
    p

}




##' Compare mean genetic signatures between pathogens
##'
##' @export
##'
##' @param store A gensig output dataframe
##'
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes_string
##' @importFrom ggplot2 geom_bar

vis.mean <- function(store) {

    df <- plyr::ldply(by(store$gensig, store$disease, mean), data.frame)
    names(df) <- c("Pathogen", "mean")
    df$Pathogen <- factor(df$Pathogen, levels = sort.gensig(store))

    p <- ggplot(df, aes_string('Pathogen', 'mean', fill = 'Pathogen')) +
        geom_bar(stat='identity')
    p

}
