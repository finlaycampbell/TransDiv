##====== Prepare workspace ======##

## Libraries
load.libs <- function() {

  libs <- c('ggplot2', 'reshape2', 'scales', 'plyr', 'dplyr', 'magrittr',
            'tidyr', 'ggrepel', 'lazyeval', 'nls2', 'ape', 'phybreak',
            'gridExtra', 'rlang', 'xtable')

  for(i in libs) {
    if(!require(i, character.only = TRUE)) {
      install.packages(i)
    }
  }

  if(!require('outbreaker2')) {
    devtools::install_github('reconhub/outbreaker2')
  }
  
}

## Load results
load.store <- function() {

  dirs <- c('~/OneDrive/output/gensig/rev/phyb/', 
            '~/OneDrive/output/gensig/rev/outb/',
            '~/OneDrive/output/gensig/rev/o.p/')

  for(i in dirs) load(paste0(i, 'store.RData'), envir = .GlobalEnv)

}



##===== Define parameters, set up run =====##

## Set up parameters for cluster run
## Describe epidemiological parameters of diseases
pathogens.param <- function(param = NULL) {
  
  param <- list(
    ebola = list(R0 = 1.8 , mut = 3.10e-6 , seql = 18958   , w.mean = 14.4 , w.sd = 8.9,  dist = "gamma" ),
    sars  = list(R0 = 2.7 , mut = 1.14e-5 , seql = 29714   , w.mean = 8.7  , w.sd = 3.6,  dist = "gamma" ),
    mers  = list(R0 = 1.2 , mut = 0.25e-5 , seql = 30115   , w.mean = 10.7 , w.sd = 6.0,  dist = "gamma" ),
    ifz   = list(R0 = 1.5 , mut = 1.19e-5 , seql = 13155   , w.mean = 3.0  , w.sd = 1.5,  dist = "gamma" ),
    mrsa  = list(R0 = 1.3 , mut = 5.21e-9 , seql = 2842618 , w.mean = 15.6 , w.sd = 10.0, dist = "gamma" ),
    klebs = list(R0 = 2.0 , mut = 6.30e-9 , seql = 5305677 , w.mean = 62.7 , w.sd = 24.0, dist = "gamma" ),
    strep = list(R0 = 1.4 , mut = 5.44e-9 , seql = 2126652 , w.mean = 6.6  , w.sd = 1.8 , dist = "gamma" ),
    shig  = list(R0 = 1.1 , mut = 1.64e-9 , seql = 4825265 , w.mean = 7.0  , w.sd = 2   , dist = "gamma" ),
    tb    = list(R0 = 1.8 , mut = 2.36e-10, seql = 4411621 , w.mean = 324  , w.sd = 385 , dist = "gamma" ),
    cdif  = list(R0 = 1.5 , mut = 8.76e-10, seql = 4290252 , w.mean = 27.7 , w.sd = 14.9, dist = 'gamma' ))

  param$klebs$mut <- param$klebs$mut*7
  for(i in c("w.mean", "w.sd")) param$klebs[[i]] <- param$klebs[[i]]/7

  param$tb$mut <- param$tb$mut*7
  for(i in c("w.mean", "w.sd")) param$tb[[i]] <- param$tb[[i]]/7

  for(disease in names(param)) param[[disease]]$w <- 0

  ## Calculate discretised distributions from the mean and standard deviations
  for(pathogen in names(param)) {
    tmp <- param[[pathogen]]
    if (tmp$dist == "weib") param[[pathogen]]$w <- discr.weib(tmp$w.mean, tmp$w.sd)
    if (tmp$dist == "gamma") param[[pathogen]]$w <- discr.gamma(tmp$w.mean, tmp$w.sd)
    
    ## Find first w that hits 0 (from decreasing slope), set the rest to 0
    ## Avoids confusion when w has lots of low but non-zero numbers
    ind <- which(c(FALSE, diff(param[[pathogen]]$w)) < 0 & param[[pathogen]]$w < 1e-10)
    if(length(ind) == 0) {
      tend <- length(param[[pathogen]]$w)
    } else {
      tend <- min(ind)
    }
    param[[pathogen]]$w[tend:length(param[[pathogen]]$w)] <- 1e-50

    param[[pathogen]]$sd.mean <- tmp$w.sd/tmp$w.mean
    
  }

  return(param)

}

## Describe parameter values for other runs
create.config <- function() {

  return(list(n.hosts = 20,
              dur = 500,
              imp = 0,
              min.n = 5))

}




##====== Analysis =====##
## Most functions are written once for outbreaker, and once for phybreak

## Main run: simulate outbreaks, and reconstruct them
run.analysis <- function(disease, runs) {

  config <- create.config()
  param <- pathogens.param()[[disease]]
  
  param$w[param$w == 0] <- 1e-50
  param$disease <- disease

  store <- list()

  for(i in seq_len(runs)) {

    sim <- run.sim(param, config)

    dna.outb.data <- list(w_dens = param$w,
                          dna = sim$dna,
                          dates = sim$onset)

    nodna.outb.data <- dna.outb.data
    nodna.outb.data$dna <- NULL

    outb.config <- list(n_iter = 1e5, sample_every = 200, find_import=FALSE,
                        move_kappa = FALSE, move_pi = FALSE, init_kappa = 1,
                        init_pi = 1, max_kappa = 1)

    dna.outb.result <- outbreaker2::outbreaker(dna.outb.data, outb.config)
    nodna.outb.result <- outbreaker2::outbreaker(nodna.outb.data, outb.config)

    store$param <- param
    store$config <- config
    store$gensig[[i]] <- get.gensig(sim)
    store$sim[[i]] <- sim
    store$sim[[i]]$dna <- NULL
    
    store$dna.result[[i]] <- dna.outb.result
    store$nodna.result[[i]] <- nodna.outb.result

    store$dna.acc[[i]] <- get.acc(dna.outb.result, sim)
    store$nodna.acc[[i]] <- get.acc(nodna.outb.result, sim)

    store$uniq[[i]] <- get.uniq(sim)
    
  }

  return(store)
  
}
run.phyb.analysis <- function(disease, runs) {

  config <- create.config()
  param <- pathogens.param()[[disease]]
  
  param$w[param$w == 0] <- 1e-50
  param$disease <- disease

  store <- list()

  for(i in seq_len(runs)) {
    
    param$shape.gen <- param$w.mean^2/param$w.sd^2
    param$mean.sample <- param$w.mean
    param$shape.sample <- param$w.mean^2/param$w.sd^2

    if(param$seql > 1e6) {
      param$seql <- round(param$seql/100, 0)
      param$mut <- param$mut*100
    }
    
    sim <- run.phyb.sim(param, config)

    dna.res <- phybreak(sim,
                        gen.mean = param$w.mean,
                        gen.shape = param$shape.gen,
                        sample.mean = param$mean.sample,
                        sample.shape = param$shape.sample,
                        est.gen.mean = FALSE,
                        est.sample.mean = FALSE,
                        est.wh.slope = FALSE) %>%
      burnin.phybreak(ncycles = 1e3) %>% 
      sample.phybreak(nsample = 500, thin = 20) 

    ## for the no.dna run make all sequences the same
    nodna.sim <- sim
    for(j in seq_along(nodna.sim$sequences)) {
      nodna.sim$sequences[[j]] <- nodna.sim$sequences[[1]]
    }
    
    nodna.res <- phybreak(nodna.sim,
                          gen.mean = param$w.mean,
                          gen.shape = param$shape.gen,
                          sample.mean = param$mean.sample,
                          sample.shape = param$shape.sample,
                          est.gen.mean = FALSE,
                          est.sample.mean = FALSE,
                          est.wh.slope = FALSE) %>%
      burnin.phybreak(ncycles = 1e3) %>% 
      sample.phybreak(nsample = 500, thin = 20) 
    
    store$param <- param
    store$config <- config
    
    store$phyb.sim[[i]] <- sim
    store$phyb.gensig[[i]] <- get.phyb.gensig(sim)
    
    store$phyb.dna.res[[i]] <- dna.res
    store$phyb.nodna.res[[i]] <- nodna.res
    
    store$phyb.dna.acc[[i]] <- get.phyb.acc(dna.res, sim)
    store$phyb.nodna.acc[[i]] <- get.phyb.acc(nodna.res, sim)
    
    store$phyb.uniq[[i]] <- get.phyb.uniq(sim)

  }

  return(store)
}

## Returns the accuracy of transmission tree inference
get.acc <- function(res, sim, burnin = 1000) {

  inferred <- summary(res, burnin = burnin)$tree$from
  true <- sim$ances

  inferred[is.na(inferred)] <- 0
  true[is.na(true)] <- 0

  acc <- mean(inferred == true)

  return(acc)

}
get.phyb.acc <- function(res, sim) {

  mean(transtree(res, "edmonds")$infector == sim$sim.infectors)

}

## Calculate the mean entropy of an outbreak run
get.ent <- function(res) {

  i <- grep('alpha', names(res))
  
  ent <- sapply(res[i], calc.ent)

  return(round(mean(ent), 2))
  
}
get.phyb.ent <- function(res) {

  nsamples <- res$d$nsamples
  obsize <- res$p$obs
  samplerange <- seq_along(res$s$mu)

  ## Extract posterior distribution of ancestries (col = samples, row = id)
  mat <- res$s$nodehosts[nsamples:(nsamples + obsize - 1), samplerange]

  out <- round(mean(apply(mat, 1, calc.ent)), 2)

  return(out)
  
}

## Run the simulation
run.sim <- function(param, config) {
 
  n <- 0
  while(n < config$min.n) {
    sim <- outbreaker::simOutbreak(R0 = param$R0,
                                   infec.curve = param$w,
                                   seq.length = param$seql,
                                   mu.transi = param$mut*(2/3),
                                   mu.transv = param$mut*(1/3),
                                   n.hosts = config$n.hosts,
                                   duration = config$dur,
                                   rate.import.case = config$imp)

    sim$inf <- sim$onset
    
    sim$onset <- sim$onset + sample(0:300,
                                    length(sim$onset),
                                    replace = TRUE,
                                    prob = param$w)
    
    n <- sim$n
  }

  return(sim)

}
run.phyb.sim <- function(param, config) {

  n <- 0

  while(n < config$min.n) {
      
    sim <- sim.phybreak(obsize = NA,
                        popsize = config$n.hosts,
                        R0 = param$R0,
                        mean.gen = param$w.mean,
                        shape.gen = param$shape.gen,
                        mean.sample = param$mean.sample,
                        shape.sample = param$shape.sample,
                        sequence.length = param$seql,
                        mu = param$mut,
                        wh.model = 3,
                        wh.slope = 1)

    ## Check for no-transmission
    if(inherits(sim, 'character')) next
    
    n <- length(sim$sample.hosts)

  }

  return(sim)

}

## Calculate the genetic signature of a simulation
get.gensig <- function(sim) {

  ids <- sim$id[!is.na(sim$ances)]
  gensig <- sapply(ids, function(id) {
    dna <- sim$dna[c(id, sim$ances[id]),]
    gensig <- as.numeric(ape::dist.dna(dna, model="N"))
  })

  return(gensig)

}
get.phyb.gensig <- function(sim) {

  ids <- which(sim$sim.infectors != "index")
  all.dna <- as.DNAbin(sim$sequences)
  gensig <- sapply(ids, function(id) {
    dna <- all.dna[c(id, which(sim$sample.hosts == sim$sim.infectors[id])),]
    gensig <- as.numeric(ape::dist.dna(dna, model="N"))
  })

  return(as.vector(gensig))

}

## Calculate the proportion of unique sequences
get.uniq <- function(sim) {

  sim$dna <- phangorn::as.phyDat(sim$dna)
  length(unique(sim$dna))/sim$n
  
}
get.phyb.uniq <- function(sim) {

  length(unique(sim$sequences))/length(sim$sample.hosts)

}

## Get the mean serial interval of a simulation
get.w <- function(sim) {

  out <- sim$onset - sim$onset[sim$ances]
  out <- out[!is.na(out)]
  return(out)

}
get.phyb.w <- function(sim) {

  ances.inftime <- sim$sample.times[match(sim$sim.infectors, sim$sample.hosts)]
  out <- sim$sample.times - ances.inftime
  out <- as.vector(out[!is.na(out)])
  return(out)
  
}


##===== Collate useful results into one object =====##
## Create a storage vector for gensig values
create.store <- function(dir, mod = 'ob') {

  cur.wd <- getwd()
  on.exit(setwd(cur.wd))

  store <- list(gensig = data.frame(disease = character(),
                                    gensig = numeric(),
                                    dat = character()),
                acc = data.frame(disease = character(),
                                 acc = numeric(),
                                 dat = character(),
                                 gensig = character()))

  files <- list.files(dir)
  rem <- grep("store", files)
  if(length(rem) > 0) files <- files[-grep("store", files)]
  n.files <- length(files)

  if(mod == 'ob') {
    adder <- add_r
  } else if(mod == 'phyb') {
    adder <- add_phyb_r
  }
  
  for(i in seq_along(files)) {
    load(paste0(dir, files[i]))
    store <- adder(store, r)
  }
  
  for(i in names(store)) {

    store[[i]]$disease <- factor(store[[i]]$disease,
                                 levels = sort.gensig(store))

    store[[i]] <- filter(store[[i]], disease != 'ifz.h')

  }

  return(store)

}

## Extract useful information from raw results
add_r <- function(store, r) {

  gensig <- data.frame(disease = r$param$disease,
                       gensig = unlist(r$gensig))

  acc <- data.frame(disease = r$param$disease,
                    dna = r$dna.acc,
                    nodna = r$nodna.acc)

  acc$improv <- with(acc, dna - nodna)
  acc$gensig <- ldply(r$gensig, mean)$V1
  acc$prop <- unlist(ldply(r$gensig, get.prop))
  acc$uniq <- r$uniq

  acc$dna.ent <- sapply(seq_along(r$dna.res),
                        function(i) get.ent(r$dna.res[[i]]))
  
  acc$nodna.ent <- sapply(seq_along(r$nodna.res),
                          function(i) get.ent(r$nodna.res[[i]]))
  
  acc$w.neg <- sapply(seq_along(r$sim),
                      function(i) mean(get.w(r$sim[[i]]) < 0))
  
  acc$n <- sapply(seq_along(r$sim),
                  function(i) length(r$sim[[i]]$n))

  acc$improv.ent <- with(acc, dna.ent - nodna.ent)
  
  store$gensig <- rbind(store$gensig, gensig)
  store$acc <- rbind(store$acc, acc)

  return(store)

}
add_phyb_r <- function(store, r) {

  gensig <- data.frame(disease = r$param$disease,
                       gensig = unlist(r$phyb.gensig))

  acc <- data.frame(disease = r$param$disease,
                    dna = r$phyb.dna.acc,
                    nodna = r$phyb.nodna.acc)

  acc$improv <- with(acc, dna - nodna)
  acc$gensig <- ldply(r$phyb.gensig, mean)$V1
  acc$prop <- unlist(ldply(r$phyb.gensig, get.prop))
  acc$uniq <- r$phyb.uniq

  acc$dna.ent <- sapply(seq_along(r$phyb.dna.res),
                        function(i) get.phyb.ent(r$phyb.dna.res[[i]]))
  
  acc$nodna.ent <- sapply(seq_along(r$phyb.nodna.res),
                          function(i) get.phyb.ent(r$phyb.nodna.res[[i]]))
  
  acc$w.neg <- sapply(seq_along(r$phyb.sim),
                      function(i) mean(get.phyb.w(r$phyb.sim[[i]]) < 0))
  
  acc$n <- sapply(seq_along(r$phyb.sim),
                  function(i) length(r$phyb.sim[[i]]$sample.hosts))

  acc$improv.ent <- with(acc, dna.ent - nodna.ent)
  
  store$gensig <- rbind(store$gensig, gensig)
  store$acc <- rbind(store$acc, acc)

  return(store)

}



##===== Helper functions =====##
## Calculate the entropy of a vector
calc.ent <- function(x) {
  x <- as.character(x)
  fk <- table(x)/sum(table(x))
  return(-sum(log(fk)*fk))
}

## Calculate the proportion of cases with gensig > 0
get.prop <- function(i) {
  mean(i > 0)
}

## Returns a vector of pathogen names ordered by mean genetic signature
sort.gensig <- function(store) {

  means <- by(store$gensig$gensig, store$gensig$disease, mean)
  return(names(sort(-means)))

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
  w <- stats::dweibull(0:300, shape = shape, scale = scale)
  return(w)
}
