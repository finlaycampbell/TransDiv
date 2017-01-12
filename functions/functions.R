###### CALCULATING GENETIC SIGNATURES OF TRANSMISSION #######

#===== Setting up the run ======#

## Describe pathogen parameters
create.param <- function() {

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

    label <- c("ebola" = "Ebola",
               "sars" = "SARS",
               "mers" = "MERS",
               "ifz" = "Influenza",
               "mrsa" = "MRSA",
               "klebs" = "Klebsiella",
               "strep" = "Streptococcus")

    param <- list(
        ebola = list(R0 = 5,0 , mut = 3.40e-6*(2/3) , seql = 18058   , w.mean = 15.3 , w.sd = 9.3 ),
        sars  = list(R0 = 3.5 , mut = 1.14e-5*(2/3) , seql = 29750   , w.mean = 8.4  , w.sd = 3.8 ),
        mers  = list(R0 = 3.5 , mut = 0.25e-5*(2/3) , seql = 30115   , w.mean = 10.7 , w.sd = 6.0 ),
        ifz   = list(R0 = 1.3 , mut = 1.19e-5*(2/3) , seql = 13155   , w.mean = 3.0  , w.sd = 1.5 ),
        mrsa  = list(R0 = 1.3 , mut = 3.58e-9*(2/3) , seql = 2842618 , w.mean = 15.6 , w.sd = 10.0),
        klebs = list(R0 = 2.0 , mut = 6.40e-9*(2/3) , seql = 5305677 , w.mean = 62.7 , w.sd = 24.0),
        strep = list(R0 = 2.0 , mut = 5.44e-9*(2/3) , seql = 2126652 , w.mean = 6.6  , w.sd = 1.8 ),
        const = list(n.hosts = 200 , dur = 500 , imp = 0.05 , min.n = 50),
        label = label
    )

    ## Calculate discretised distributions from the mean and standard deviations
    for(pathogen in names(label)) {
        par <- param[[pathogen]]
        if(pathogen=="sars") param[[pathogen]]$w <- discr.weib(par$w.mean, par$w.sd)
        else param[[pathogen]]$w <- discr.gamma(par$w.mean, par$w.sd)
    }

    return(param)

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

## Return a vector of realised gensig values for a given disease
run.gensig <- function(disease) {

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
    ids <- sim$id[!is.na(sim$ances)]

    ## Return a vector with individual gensig values
    gensig.vec <- sapply(ids, calc.gensig)

    print(disease)

    return(gensig.vec)

}

## run.gensig over the various pathogens
path.gensig <- function(param) {

    param <- create.param()

    names <- names(param)[-which(names(param)=="const")]

    df <- reshape2::melt(sapply(names,run.gensig))

    names(df) <- c("gensig","disease")

    return(df)

}


#====== Collect cluster results =====#

## Create a storage vector for gensig values
create.store <- function(obj,bundle.name,dir,load=FALSE,dl=TRUE) {

    if(load) load(paste0(dir,"store.RData"))
    else store <- NULL

    if(dl) {
        task_bundle <- obj$task_bundle_get(bundle.name)
        ids <- task_bundle$ids
        pb <- txtProgressBar(min=1,max=length(ids),style=3)
        for(i in seq_along(ids)){
            setTxtProgressBar(pb,i)
            task <- obj$task_get(ids[i])
            r <- task$result()
            store <- rbind(store,r)
        }
    }

    return(store)
}

## Create a table of parameter values for .Rmd
create.param.table <- function(store) {

    ## A function for quick accessing of parameter values
    access <- function(pathogen,factor) param[[pathogen]][[factor]]

    param <- create.param()

    ## Use consistent ordering of pathogen names throughout
    paths <- sort.gensig(store)

    df <- data.frame(matrix(nrow = length(paths), ncol = 4))
    names(df) <- c("path","R0","mut","seql")

    ## Access parameter values from param and fill them into df
    df$path <- param$label[paths]

    ## Fill in R0, mut and seql
    for(factor in names(df)[2:4])
        df[[factor]] <- sapply(paths, access, factor)

    ## Insert the generation time information in format: mean (sd)
    df$w <- sapply(paths, function(path)
        paste0(access(path,"w.mean"), " (",access(path,"w.sd"),")"))

    ## Calculate the expected genetic signature
    #df$prod <- sapply(paths, function(path)
    #    round(access(path, "mut") * access(path, "seql") * access(path, "w.mean"), 2))

    df$mut <- format(df$mut, digits = 3, scientific = TRUE)
    df$seql <- format(df$seql, digits = 3, scientific = TRUE)

    names(df) <- c("Pathogen",
                   "R0<br>",
                   "Mutation rate<br>(base<sup>-1</sup> day<sup>-1</sup>)",
                   "Genome length<br> (bases)",
                   "Generation time (SD)<br>(days) ")

    return(df)

}


#===== Analyse results ======#

## Returns a vector of pathogen names ordered by mean genetic signature (high to low)
sort.gensig <- function(store) {

    means <- by(store$gensig,store$disease,mean)
    return(names(sort(-means)))

}


#===== Plot results =====#

## Plot generation time distributions
plot.w <- function(param) {

    names <- names(param)[-which(names(param) %in% c("const","label"))]

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
plot.dist <- function(store) {

    ## Create param to access labeller
    param <- create.param()

    ## Calculate proportions manually
    df <- by(store$gensig, store$disease, function(sig) prop.table(table(sig)))

    ## Collapse the list of dataframes into a single dataframe
    df <- plyr::ldply(df, data.frame)

    ## Sort the pathogen names by mean gensig
    df$.id <- factor(df$.id, levels = sort.gensig(store))

    p <- ggplot(df,aes(x = sig, y = Freq, fill = .id, colour = .id)) +
        geom_bar(stat = "identity") +
        facet_wrap( ~ .id,labeller = as_labeller(param$label)) +
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
