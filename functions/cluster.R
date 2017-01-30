#===== Libraries =====#

library(ggplot2)
library(reshape2)
library(scales)

#===== Setting directory and loading functions =====#

setwd("Q:\\gensig")
source("functions.R")

#===== Setting up the cluster context and object =====#

didewin::didewin_config_global(credentials="fc1915",cluster="fi--dideclusthn")
didewin::web_login()

our.pkgs <- c("ggplot2","reshape2","outbreaker","plyr","ape","scales",'EpiEstim')

our.sources <- "functions.R"

ctx <- context::context_save("contexts",packages=our.pkgs,sources=our.sources)

obj <- didewin::queue_didewin(ctx)


#===== Sending off the tasks =====#

param <- data.frame(run = 1:20)

grp <- queuer::enqueue_bulk(obj, param, path.gensig, do.call=FALSE, timeout=0)


#====== Checking the status ======#
bundle.name <- 'brimstone_blueandgoldmackaw'
bundle <- obj$task_bundle_get(bundle.name)
bundle$status()

#===== Accessing runs on the cluster  =====#
dir <- "C:\\Users\\Finlay\\OneDrive - Imperial College London\\R\\gensig\\output\\"
store <- create.store(obj,bundle.name,dir,load=TRUE,dl=TRUE)
