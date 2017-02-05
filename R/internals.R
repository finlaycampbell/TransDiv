###### CALCULATING GENETIC SIGNATURES OF TRANSMISSION #######

## Describe pathogen parameters
create.param <- function(param = NULL) {

    defaults <- list(
        ebola = list(R0 = 1.7 , mut = 3.10e-6*(2/3) , seql = 18958   , w.mean = 14.4 , w.sd = 8.9,  dist = "gamma" ),
        sars  = list(R0 = 2.7 , mut = 1.14e-5*(2/3) , seql = 29714   , w.mean = 8.7  , w.sd = 3.6,  dist = "weib"  ),
        mers  = list(R0 = 3.5 , mut = 0.25e-5*(2/3) , seql = 30115   , w.mean = 10.7 , w.sd = 6.0,  dist = "gamma" ),
        ifz   = list(R0 = 1.3 , mut = 1.19e-5*(2/3) , seql = 13155   , w.mean = 3.0  , w.sd = 1.5,  dist = "gamma" ),
        mrsa  = list(R0 = 1.3 , mut = 3.58e-9*(2/3) , seql = 2842618 , w.mean = 15.6 , w.sd = 10.0, dist = "gamma" ),
        klebs = list(R0 = 2.0 , mut = 6.30e-9*(2/3) , seql = 5305677 , w.mean = 62.7 , w.sd = 24.0, dist = "gamma" ),
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
        tmp <- param[[pathogen]]
        if (tmp$dist == "weib") param[[pathogen]]$w <- discr.weib(tmp$w.mean, tmp$w.sd)
        if (tmp$dist == "gamma") param[[pathogen]]$w <- discr.gamma(tmp$w.mean, tmp$w.sd)
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
    w <- stats::dweibull(1:100,shape=shape,scale=scale)
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
