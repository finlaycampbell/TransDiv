#' Distribution of genetic signatures of common pathogens
#'
#' Genetic signatures were simulated using the gensig function for the following
#' pathogens: Ebola, SARS-CoV, MERS-CoV, Influenza A (H1N1), MRSA, Klebsiella
#' pneumoniae and Streptococcus pneumoniae. Simulations were run on the DIDE
#' computer cluster, and the empirical genetic signature of >2000 transmission
#' pairs per pathogen determined.
#'
#' @docType data
#'
#' @format {
#' A data frame with 32898 rows and 2 columns
#' \describe{
#'   \item{disease}{Disease}
#'   \item{gensig}{Genetic signature of a single transmission event}
#' }
#' }
#'
#' @rdname store
#'
#' @examples
#' ## visualise distribution
#' plot.dist(store)
#'
"store"
