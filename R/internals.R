#====== Running the simulations =====#

## Set up parameters for cluster run
cluster.param <- function() {
    diseases <- names(create.param())
    param <- data.frame(disease = rep(diseases, each = 5))
    param$min.n <- 30
    param$n.hosts <- 100
    param$runs <- 20
    param$imp <- 0
    param$dur <- 500
    param$disease <- as.character(param$disease)
    return(param)
}

## Run gensig, then run outbreaker2 on the simulation
run.cluster <- function(par) {

    config <- create.config(min.n = par$min.n,
                            n.hosts = par$n.hosts,
                            imp = par$imp,
                            dur = par$dur)

    disease <- par$disease

    param <- create.param()
    param <- list(tmp = param[[disease]])
    names(param) <- disease
    param[[disease]]$w[param[[1]]$w == 0] <- 1e-50

    store <- list()

    for(i in seq_len(par$runs)) {

        gensig.result <- gensig(param = param, config = config)

        store$param <- gensig.result$param
        store$config <- gensig.result$config
        store$gensig[[i]] <- gensig.result$store[[disease]]$gensig

        sim <- gensig.result$store[[disease]]$sim
        store$sim[[i]] <- sim
        store$sim[[i]]$dna <- NULL

        dna.outb.data <- list(w.dens = param[[disease]]$w,
                              dna = sim$dna,
                              dates = sim$onset)

        nodna.outb.data <- dna.outb.data
        nodna.outb.data$dna <- NULL

        outb.config <- list(n.iter = 1e5, sample.every = 200, find.import=FALSE,
                            move.kappa = FALSE, move.pi = FALSE, init.kappa = 1,
                            init.pi = 1, max.kappa = 1)

        dna.outb.result <- outbreaker2::outbreaker(dna.outb.data, outb.config)
        nodna.outb.result <- outbreaker2::outbreaker(nodna.outb.data, outb.config)

        dna.acc <- get.acc(dna.outb.result, sim)
        nodna.acc <- get.acc(nodna.outb.result, sim)

        store$dna.result[[i]] <- dna.outb.result
        store$nodna.result[[i]] <- nodna.outb.result

        store$dna.acc[[i]] <- dna.acc
        store$nodna.acc[[i]] <- nodna.acc

    }

    return(store)

}

## Returns the accuracy of transmission tree inference
get.acc <- function(result, outbreak) {

    id <- seq_len(outbreak$n)
    adder <- which(names(result)=="alpha.1") - 1
    samples <- length(result$step)

    #Determine the modal transmission network
    network <- data.frame(from=do.call(rbind, lapply(id,  function(i) ({
        modal.ances <- as.integer(names(which.max(table(result[[i+adder]]))))
        if(length(modal.ances)==0) return(NA) else return(modal.ances)
    }))),  to=id)

    import <- which(is.na(network$from))

    transmission.id <- id[!sapply(result[id+adder], function(i) any(is.na(i)))]

    #Determine the proportion of correctly inferred ancestries
    num.correct <-  sum(outbreak$ances==network$from, na.rm=TRUE)
    num.correct <- num.correct + sum(is.na(outbreak$ances[is.na(network$from)]))
    acc <- round(num.correct/nrow(network), 2)

    return(acc)

}

## Get the mean generation time
get.w <- function(sim) {

    out <- sim$onset - sim$onset[sim$ances]
    out <- out[!is.na(out)]
    return(out)

}

## Main run per pathogen
gensig <- function(param = NULL, config = NULL) {

    param <- create.param(param)

    config <- create.config(config)

    store <- list()

    for(disease in names(param)) {

        sim <- run.sim(disease, param, config)

        store[[disease]][["sim"]] <- sim

        store[[disease]][["gensig"]] <- get.gensig(sim)

    }

    return(list(param = param, config = config, store = store))

}

## Epidemiological and evolutionary parameters of pathogens
create.param <- function(param = NULL) {

    defaults <- list(
        ebola = list(R0 = 1.8 , mut = 3.10e-6 , seql = 18958   , w.mean = 14.4 , w.sd = 8.9,  dist = "gamma" ),
        sars  = list(R0 = 2.7 , mut = 1.14e-5 , seql = 29714   , w.mean = 8.7  , w.sd = 3.6,  dist = "gamma" ),
        mers  = list(R0 = 1.2 , mut = 0.25e-5 , seql = 30115   , w.mean = 10.7 , w.sd = 6.0,  dist = "gamma" ),
        ifz   = list(R0 = 1.5 , mut = 1.19e-5 , seql = 13155   , w.mean = 3.0  , w.sd = 1.5,  dist = "gamma" ),
       #ifz.h = list(R0 = 1.5 , mut = 1.09e-5 , seql = 1701    , w.mean = 3.0  , w.sd = 1.5,  dist = "gamma" ),
        mrsa  = list(R0 = 1.3 , mut = 5.21e-9 , seql = 2842618 , w.mean = 15.6 , w.sd = 10.0, dist = "gamma" ),
        klebs = list(R0 = 2.0 , mut = 6.30e-9 , seql = 5305677 , w.mean = 62.7 , w.sd = 24.0, dist = "gamma" ),
        strep = list(R0 = 1.4 , mut = 5.44e-9 , seql = 2126652 , w.mean = 6.6  , w.sd = 1.8 , dist = "gamma" ),
        shig  = list(R0 = 1.1 , mut = 1.64e-9 , seql = 4825265 , w.mean = 7.0  , w.sd = 2   , dist = "gamma" ),
        tb    = list(R0 = 1.8 , mut = 2.36e-10, seql = 4411621 , w.mean = 324  , w.sd = 385 , dist = "gamma" ),
        cdif  = list(R0 = 1.5 , mut = 8.76e-10, seql = 4290252 , w.mean = 27.7 , w.sd = 14.9, dist = 'gamma' ))

    defaults$tb$w.sd <- defaults$tb$w.sd/2

    defaults$klebs$mut <- defaults$klebs$mut*7
    for(i in c("w.mean", "w.sd")) defaults$klebs[[i]] <- defaults$klebs[[i]]/7


    defaults$tb$mut <- defaults$tb$mut*7
    for(i in c("w.mean", "w.sd")) defaults$tb[[i]] <- defaults$tb[[i]]/7

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
    w <- sapply(0:300, EpiEstim::DiscrSI, mean, sd)
    return(w)
}

## Returns a discretised weibull distribution
discr.weib <- function(mean, sd) {
    shape <- (sd/mean)^-1.086
    scale <- mean/gamma(1+1/shape)
    w <- stats::dweibull(0:200, shape = shape, scale = scale)
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

## Calculate the proportion of cases with gensig > 0
get.prop <- function(i) {
    1 - table(i)[1]/sum(table(i))
}

## Run the simulations
run.sim <- function(disease, param, config) {

    tmp <- param[[disease]]

    n <- 0
    while(n < config$min.n) {
        sim <- outbreaker::simOutbreak(R0 = tmp$R0,
                                       infec.curve = tmp$w,
                                       seq.length = tmp$seql,
                                       mu.transi = tmp$mut*(2/3),
                                       mu.transv = tmp$mut*(1/3),
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


#===== Analysing results =====#

## Create a storage vector for gensig values
create.store <- function(dir) {

    store <- list(gensig = data.frame(disease = character(),
                                      gensig = numeric(),
                                      dat = character()),
                  acc = data.frame(disease = character(),
                                   acc = numeric(),
                                   dat = character(),
                                   gensig = character()))

    files <- list.files(dir)
    n.files <- length(files)

    pb <- txtProgressBar(min = 1, max = length(files), style = 3)
    for(i in seq_along(files)) {
        setTxtProgressBar(pb, i)
        load(paste0(dir, files[i]))
        store <- add_r(store, r)
    }

    for(i in names(store)) {

        store[[i]]$disease <- factor(store[[i]]$disease,
                                     levels = rev(sort.gensig(store)))

        store[[i]] <- filter(store[[i]], disease != 'ifz.h')
    }

    return(store)

}

## Sort gensig from highest to lowest
sort.gensig <- function(store) {

    store$gensig %>%
        group_by(disease) %>%
        summarise(gensig = mean(gensig)) %$%
        disease[order(gensig)]

}

## Create a labeller for linking variable names with disease labels
create.lab <- function(input) {

    c(ebola = "EBOV",
      sars = "SARS-CoV",
      mers = "MERS-CoV",
      ifz = "Influenza A",
      mrsa = "MRSA",
      klebs = "K. pneumoniae",
      strep = "S. pneumoniae",
      shig = "S. sonnei",
      tb = "M. tuberculosis",
      cdif = "C. difficile")

}

## Create a labeller for x and y axes from variable names
create.axlab <- function(input) {

    reference <- c(gensig = 'Average transmission divergence of outbreak',
                   improv = 'Change in accuracy of outbreak reconstruction',
                   dna = 'Accuracy of outbreak reconstruction',
                   nodna = 'Accuracy of outbreak reconstruction',
                   prop = 'Proportion of genetically distinct transmission pairs')

    return(reference[[input]])

}

## Adds r loaded from cluster or file to store
add_r <- function(store, r) {

    gensig <- data.frame(disease = names(r$param),
                         gensig = unlist(r$gensig))

    acc <- data.frame(disease = names(r$param),
                      dna = r$dna.acc,
                      nodna = r$nodna.acc)

    sim <- data.frame(disease = names(r$param),
                      n = sapply(seq_along(r$sim), function(i) r$sim[[i]]$n))

    acc$improv <- with(acc, dna - nodna)
    acc$gensig <- ldply(r$gensig, mean)$V1
    acc$prop <- unlist(ldply(r$gensig, get.prop))

    store$gensig <- rbind(store$gensig, gensig)
    store$acc <- rbind(store$acc, acc)
    store$sim <- rbind(store$sim, sim)

    return(store)

}

## Calculate the proportion of cases with gensig > 0
get.prop <- function(i) {
    1 - table(i)[1]/sum(table(i))
}

## Returns a vector of pathogen names ordered by mean genetic signature (high to low)
sort.gensig <- function(store) {

    means <- by(store$gensig$gensig, store$gensig$disease, mean)
    return(names(sort(-means)))

}

## Proto package functions that weren't downloading properly
as.proto.list <- function(x, envir, parent, all.names = FALSE, ...,
                          funEnvir = envir, SELECT = function(x) TRUE) {
  if (missing(envir)) {
    if (missing(parent))
      parent <- parent.frame()
    envir <- if (is.proto(parent))
      parent$proto(...)
    else
      proto(parent, ...)
  }
  for (s in names(x))
    if (SELECT(x[[s]])) {
      assign(s, x[[s]], envir = envir)
      if (is.function(x[[s]]) && !identical(funEnvir, FALSE))
        environment(envir[[s]]) <- funEnvir
    }
  if (!missing(parent))
    parent.env(envir) <- parent
  as.proto.environment(envir)  # force refresh of .that and .super
}

as.proto.environment <- function(x, ...) {
  assign(".that", x, envir = x)
  assign(".super", parent.env(x), envir = x)
  structure(x, class = c("proto", "environment"))
}

as.lm <- function(object, ...) UseMethod("as.lm")

as.lm.nls <- function(object, ...) {
    if (!inherits(object, "nls")) {
		w <- paste("expected object of class nls but got object of class:",
			paste(class(object), collapse = " "))
		warning(w)
	}

	gradient <- object$m$gradient()
	if (is.null(colnames(gradient))) {
		colnames(gradient) <- names(object$m$getPars())
	}

	response.name <- if (length(formula(object)) == 2) "0" else
		as.character(formula(object)[[2]])

	lhs <- object$m$lhs()
	L <- data.frame(lhs, gradient)
	names(L)[1] <- response.name

	fo <- sprintf("%s ~ %s - 1", response.name,
		paste(colnames(gradient), collapse = "+"))
	fo <- as.formula(fo, env = as.proto.list(L))

	do.call("lm", list(fo, offset = substitute(fitted(object))))

}
