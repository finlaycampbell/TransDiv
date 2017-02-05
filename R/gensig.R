
##' Simulating genetic signatures of disease transmission
##'
##' @author Finlay Campbell <f.campbell15@@imperial.ac.uk>
##'
##' @return Returns a dataframe of two columns
##'
##' @export
##'
##' @param param a list of lists describing the parameter values for each pathogen of interest
##' @param config an optional list of parameter values for the run

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
