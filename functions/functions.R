#===== Loading dependencies and data =====#
library(ggplot2)

#===== Describe pathogen parameters =====#

make.param <- function() {

    out <- data.frame(matrix(c(

    # | R0  | Mutation rate (/day) | Seq length | Disease  |
        5   , 1.24e-3/365*(2/3)    , 18058      , # ebola
        3.5 , 1.14e-5*(2/3)        , 29750      ) # sars

        ,nrow=2,byrow=TRUE,dimnames=list(NULL,c("r0","mut","seql"))))

    out$disease <- c("ebola","sars")
    out <- out[,c(4,1,2,3)]
}
