## Loading dependencies and data
library(ggplot2)

## Describe pathogen parameters
create.param <- function() {

    ebola.w.mean <- 15.3
    ebola.w.variance <- 9.3^2
    ebola.w <- dgamma(1:50, shape = ebola.w.mean^2/ebola.w.variance,
                      rate = ebola.w.mean/ebola.w.variance)

    sars.w.mean <- 8.4
    sars.w.sd <- 3.8
    sars.weib.shape <- (sars.w.sd/sars.w.mean)^-1.086
    sars.weib.scale <- sars.w.mean/gamma(1+1/sars.weib.shape)
    sars.w <- dweibull(1:30,shape=sars.weib.shape,scale=sars.weib.scale)

    out <- list(
        ebola = list(R0 = 5   , mut = 1.24e-3/365*(2/3) , seql = 18058 , w = ebola.w ),
        sars  = list(R0 = 3.5 , mut = 1.14e-5*(2/3)     , seql = 29750 , w = sars.w  ),
        const = list(n.hosts = 50 , dur = 50 , imp = 0.01)
    )

}

## Run the simulations
run.sim <- function(disease, param) {

    temp.p <- param[[disease]]

    sim <- outbreaker::simOutbreak(R0 = temp.p$R0,
                                   infec.curve = temp.p$w,
                                   seq.length = temp.p$seql,
                                   mu.transi = temp.p$mut,
                                   n.hosts = param$const$n.hosts,
                                   duration = param$const$dur,
                                   rate.import.case = param$const$imp)

    return(sim)

}

## Create a frequency table of R values
calc.R <- function(sim) {

    R <- as.data.frame(table(table(sim$ances)),stringsAsFactors = F)
    names(R) <- c("R","freq")
    zeroes <- c(0,sim$n - sum(R$freq))
    R <- rbind(zeroes,R)

}

## Returns a function calculating the genetic signature between an id and its
## ancestor. The simulation 'sim' is enlosed as it remains unchanged, avoiding
## unnecessary passing of the sim argument
create.calc.gensig <- function(sim) {

    out <- function(id) {
        dna <- sim$dna[c(id,sim$ances[id]),]
        gensig <- as.numeric(ape::dist.dna(dna,model="N"))
        return(gensig)
    }

}

## Create a histogram of gensig values
plot.bar <- function(gensig.vec) {

    p <- ggplot(data.frame(vec=gensig.vec),aes(vec)) +
        geom_bar(colour="white") +
        xlab("Genetic signature") + ylab("Frequency")
    p

}

## Run a simulation, calculate and plot the genetic signature
gensig <- function(disease) {

    ## Create list of epidemiological parameters
    param <- create.param()

    ## Simulate outbreak and sequence evolution
    sim <- run.sim(disease,param)

    ## Create calc.gensig with simulation data enclosed
    calc.gensig <- create.calc.gensig(sim)

    ## Remove imports with unknown ancestor for gensig calculation
    ids <- sim$id[!is.na(sim$id)]

    ## Return a vector with individual gensig values
    gensig.vec <- sapply(ids, calc.gensig)

    ## Create the bar plot of gensig values
    bar <- plot.bar(gensig.vec)

    ## Create a frequency table of gensig values
    gensig.freq <- as.data.frame(table(gensig.vec))
    names(gensig.freq) <- c("gensig","freq")

    out <- list(freq = gensig.freq, bar = bar)

    return(out)

}
