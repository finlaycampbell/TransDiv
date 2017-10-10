##===== TransDiv =====##

source("R/internals.R")
source("R/plots.R")
load.libs()

##===== Use the data from the manuscript ======##
load.store()


##===== Or run the simulations yourself =====##
## This will take a long time - reduce the number of runs
pathogens <- names(pathogens.param())
for(i in seq_along(pathogens)) {

  # Run the analysis using the outbreak model
  r <- run.analysis(pathogens[i], runs = 100)
  save(r, file = paste0('output/outbreaker/run_', i, '.RData'))

  # Run the analysis using the phybreak model
  r <- run.phyb.analysis(pathogens[i], runs = 100)
  save(r, file = paste0('output/phybreak/run_', i, '.RData'))
  
}

## Extract results from simulation data
o.store <- create.store(dir = 'output/outbreaker/', mod = 'ob')
p.store <- create.store(dir = 'output/phybreak/', mod = 'phyb')

##===== Create figures =====##
## These will be saved to the /figs folder
## See R/plots.R for code
create.figs()
