#===== TransDiv =====#

#===== Loading libraries and functions =====#
libs <- c('ggplot2', 'reshape2', 'scales', 'gensig', 'plyr', 'dplyr', 'magrittr',
          'tidyr', 'ggrepel', 'lazyeval', 'nls2', 'ape', 'vegan')
lapply(libs, require, character.only = TRUE)

source("R/internals.R")
source("R/plots.R")

#===== Use the data from the manuscript ======#
store <- create.store(dir = 'data/')

#===== OR =====#

#===== Run the simulations yourself =====#
param <- cluster.param()
for(i in nrow(param)) {
    r <- run.cluster(param[i,])
    save(r, file = paste0('output/store_', i, '.RData'))
}

store <- create.store(dir = 'output/')

#===== Create figures =====#
## Fig. 1A
vis.gensig(store)

## Fig. 1B
vis.prop(store)

## Fig. 2
vis.rel(store)

## Fig. 3
vis.atleast(store)
