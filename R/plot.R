##' Visualise distributions of genetic signatures
##'
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_line
##' @importFrom ggplot2 labs
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 geom_bar
##' @importFrom ggplot2 facet_wrap
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 as_labeller

plot.w <- function(param) {

    names <- names(param)

    df <- data.frame(days=seq_along(param[[names[1]]]$w))

    for(i in names) df[[i]] <- param[[i]]$w

    mlt <- reshape2::melt(df,id="days")
    mlt$variable <- factor(mlt$variable, levels = sort.gensig(store))

    p <- ggplot(mlt,aes(days,value,colour=variable)) + geom_line(size=1.5) +
        labs(title = "Generation Times",
             subtitle = "Generation times of various pathogens",
             x = "Days",
             y = "Density")

    return(p)

}

## Create a bar chart of gensig distributions
plot.dist <- function(store, config = NULL) {

    ## Create config to access labeller
    if(is.null(config)) config <- create.config()

    ## Calculate proportions manually
    df <- by(store$gensig, store$disease, function(sig) prop.table(table(sig)))

    ## Collapse the list of dataframes into a single dataframe
    df <- plyr::ldply(df, data.frame)

    ## Sort the pathogen names by mean gensig
    df$.id <- factor(df$.id, levels = sort.gensig(store))

    p <- ggplot(df, aes(x = sig, y = Freq, fill = .id, colour = .id)) +
        geom_bar(stat = "identity") +
        facet_wrap( ~ .id,labeller = as_labeller(config$label)) +
        theme(legend.position = "none") +
        labs(title = "Genetic signatures", x = "Number of mutations", y = "Proportion")
    p

}

## Plot a bar chart of mean values
plot.mean <- function(store) {

    df <- plyr::ldply(by(store$gensig,store$disease,mean),data.frame)
    names(df) <- c("Pathogen","mean")
    df$Pathogen <- factor(df$Pathogen,levels=sort.gensig(store))

    p <- ggplot(df,aes(Pathogen,mean,fill=Pathogen)) +
        geom_bar(stat='identity')
    p

}
