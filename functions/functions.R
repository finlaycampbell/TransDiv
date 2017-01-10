## Describe pathogen parameters
create.param <- function() {

    ebola.w.mean <- 15.3
    ebola.w.sd <- 9.3
    ebola.w <- sapply(1:100,EpiEstim::DiscrSI,ebola.w.mean,ebola.w.sd)

    sars.w.mean <- 8.4
    sars.w.sd <- 3.8
    sars.weib.shape <- (sars.w.sd/sars.w.mean)^-1.086
    sars.weib.scale <- sars.w.mean/gamma(1+1/sars.weib.shape)
    sars.w <- dweibull(1:100,shape=sars.weib.shape,scale=sars.weib.scale)

    mers.w.mean <- 10.7
    mers.w.sd <- 6.0
    mers.w <- sapply(1:100,EpiEstim::DiscrSI,mers.w.mean,mers.w.sd)

    ifz.w.mean <- 3.0
    ifz.w.sd <- 1.5
    ifz.w <- sapply(1:100,EpiEstim::DiscrSI,ifz.w.mean,ifz.w.sd)

    mrsa.w.mean <- 15.6
    mrsa.w.sd <- 10.0
    mrsa.w <- sapply(1:100,EpiEstim::DiscrSI,mrsa.w.mean,mrsa.w.sd)

    klebs.w.mean <- 62.7
    klebs.w.sd <- 24.0
    klebs.w <- sapply(1:100,EpiEstim::DiscrSI,klebs.w.mean,klebs.w.sd)

    strep.w.mean <- 6.6
    strep.w.sd <- 1.8
    strep.w <- sapply(1:100,EpiEstim::DiscrSI,strep.w.mean,strep.w.sd)

    out <- list(
        ebola = list(R0 = 5,0 , mut = 3.40e-6*(2/3) , seql = 18058   , w = ebola.w ),
        sars  = list(R0 = 3.5 , mut = 1.14e-5*(2/3) , seql = 29750   , w = sars.w  ),
        mers  = list(R0 = 3.5 , mut = 0.25e-5*(2/3) , seql = 30115   , w = mers.w  ),
        ifz   = list(R0 = 1.3 , mut = 1.19e-5*(2/3) , seql = 13155   , w = ifz.w   ),
        mrsa  = list(R0 = 1.3 , mut = 3.58e-9*(2/3) , seql = 2842618 , w = mrsa.w  ),
        klebs = list(R0 = 2.0 , mut = 6.40e-9*(2/3) , seql = 5305677 , w = klebs.w ),
        strep = list(R0 = 2.0 , mut = 5.44e-9*(2/3) , seql = 2126652 , w = strep.w ),
        const = list(n.hosts = 200 , dur = 500 , imp = 0.05 , min.n = 50)
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

## Plot generation time distributions
plot.w <- function(param) {

    names <- names(param)[-which(names(param)=="const")]

    df <- data.frame(days=seq_along(param[[names[1]]]$w))

    for(i in names) df[[i]] <- param[[i]]$w

    mlt <- reshape2::melt(df,id="days")

    p <- ggplot(mlt,aes(days,value,colour=variable)) + geom_line(size=1.5) +
        labs(title = "Generation Times",
             subtitle = "Generation times of various pathogens",
             x = "Days",
             y = "Density")

    return(p)

}

## Return a vector of realised gensig values
gensig <- function(disease) {

    ## Create list of epidemiological parameters
    param <- create.param()

    ## Simulate outbreak and sequence evolution
    repeat {
        sim <- run.sim(disease,param)
        print(sim$n)
        if(sim$n > param$const$min.n) break
    }

    ## Create calc.gensig with simulation data enclosed
    calc.gensig <- create.calc.gensig(sim)

    ## Remove imports with unknown ancestor for gensig calculation
    ids <- sim$id[!is.na(sim$id)]

    ## Return a vector with individual gensig values
    gensig.vec <- sapply(ids, calc.gensig)

    print(disease)

    return(gensig.vec)

}

## Run gensig over the pathogens provided in create.param
run.gensig <- function() {

    param <- create.param()

    names <- names(param)[-which(names(param)=="const")]

    df <- reshape2::melt(sapply(names,gensig))
    names(df) <- c("gensig","disease")

    dis.names <- c("ebola" = "Ebola",
                   "sars" = "SARS",
                   "mers" = "MERS",
                   "ifz" = "Influenza",
                   "mrsa" = "MRSA",
                   "klebs" = "Klebsiella",
                   "strep" = "Streptococcus")

    p <- ggplot(df,aes(gensig, fill = disease, colour = disease)) +
        geom_bar() + facet_wrap( ~ disease,labeller = as_labeller(dis.names)) +
        theme(legend.position = "none") +
        labs(title = "Genetic signatures",x = "Number of mutations", y = "Count")
    p

    return(list(p=p,df=df))

}
