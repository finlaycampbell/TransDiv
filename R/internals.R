###### CALCULATING GENETIC SIGNATURES OF TRANSMISSION #######

## Describe pathogen parameters
create.param <- function(param = NULL) {

    defaults <- list(
        ebola = list(R0 = 5.0 , mut = 3.40e-6*(2/3) , seql = 18058   , w.mean = 15.3 , w.sd = 9.3,  dist = "gamma" ),
        sars  = list(R0 = 3.5 , mut = 1.14e-5*(2/3) , seql = 29750   , w.mean = 8.4  , w.sd = 3.8,  dist = "weib"  ),
        mers  = list(R0 = 3.5 , mut = 0.25e-5*(2/3) , seql = 30115   , w.mean = 10.7 , w.sd = 6.0,  dist = "gamma" ),
        ifz   = list(R0 = 1.3 , mut = 1.19e-5*(2/3) , seql = 13155   , w.mean = 3.0  , w.sd = 1.5,  dist = "gamma" ),
        mrsa  = list(R0 = 1.3 , mut = 3.58e-9*(2/3) , seql = 2842618 , w.mean = 15.6 , w.sd = 10.0, dist = "gamma" ),
        klebs = list(R0 = 2.0 , mut = 6.40e-9*(2/3) , seql = 5305677 , w.mean = 62.7 , w.sd = 24.0, dist = "gamma" ),
        strep = list(R0 = 2.0 , mut = 5.44e-9*(2/3) , seql = 2126652 , w.mean = 6.6  , w.sd = 1.8 , dist = "gamma" ))

    for(disease in names(defaults)) defaults[[disease]]$w <- 0

    if(!is.null(param)) {
        for(i in 1:length(param)) {
            if(length(param[[i]]) != 7) {
                stop("Incorrect number of parameter values provided")
            }
            if(any(!names(param[[i]]) %in% names(defaults[[1]]))) {
                stop("Incorrect parameter name provided")
            }
        }
    } else {
        param <- defaults
    }

    ## Calculate discretised distributions from the mean and standard deviations
    for(pathogen in names(param)) {
        par <- param[[pathogen]]
        if (par$dist == "weib") param[[pathogen]]$w <- discr.weib(par$w.mean, par$w.sd)
        if (par$dist == "gamma") param[[pathogen]]$w <- discr.gamma(par$w.mean, par$w.sd)
    }

    return(param)

}

## Describe parameter values for other runs
create.config <- function(...) {

    config <- list(...)
    if(length(config) == 1L && is.list(config[[1]])) {
        config <- config[[1]]
    }

    ## If no user labeller is provided, simply use the given names
    label <- function(char) return(char)

    defaults <- list(n.hosts = 200,
                     dur = 500,
                     imp = 0.05,
                     min.n = 50,
                     label = label)

    config <- modify.defaults(defaults, config)

    return(config)

}

## Returns a discretized gamma distribution
discr.gamma <- function(mean, sd) {
    w <- sapply(1:100, EpiEstim::DiscrSI, mean, sd)
    return(w)
}

## Returns a discretised weibull distribution
discr.weib <- function(mean, sd) {
    shape <- (sd/mean)^-1.086
    scale <- mean/gamma(1+1/shape)
    w <- dweibull(1:100,shape=shape,scale=scale)
    return(w)
}

## Modifies default values of a function
modify.defaults <- function(default, modified) {

    not.found <- ! names(modified) %in% names(default)
    if(any(not.found)) stop(paste(paste(names(modified)[not.found], collapse = ", "),
                                  "is not a valid descriptor"))

    for (i in names(modified)) default[[i]] <- modified[[i]]
    return(default)

}

## Run the simulations
run.sim <- function(disease, param, config) {

    tmp <- param[[disease]]

    n <- 0
    while(n < config$min.n) {
        sim <- outbreaker::simOutbreak(R0 = tmp$R0,
                                       infec.curve = tmp$w,
                                       seq.length = tmp$seql,
                                       mu.transi = tmp$mut,
                                       n.hosts = config$n.hosts,
                                       duration = config$dur,
                                       rate.import.case = config$imp)
        n <- sim$n
    }

    return(sim)

}

## Returns a function calculating the genetic signature between an id and its
## ancestor. The simulation 'sim' is enlosed as it remains unchanged, avoiding
## unnecessary passing of the sim argument
get.gensig <- function(sim) {

    ids <- sim$id[!is.na(sim$ances)]
    gensig <- sapply(ids, function(id) {
        dna <- sim$dna[c(id, sim$ances[id]),]
        gensig <- as.numeric(ape::dist.dna(dna, model="N"))
    })

    return(gensig)

}

## Returns a vector of pathogen names ordered by mean genetic signature (high to low)
sort.gensig <- function(store) {

    means <- by(store$gensig,store$disease,mean)
    return(names(sort(-means)))

}


#===== Plotting functions =====#

## Plot generation time distributions
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
