% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HiC2Tree.R
\name{normalizeCount}
\alias{normalizeCount}
\title{Normalize count matrix using cyclic loess}
\usage{
normalizeCount(count_matrice)
}
\arguments{
\item{count_matrice}{The count matrix to normalize.}
}
\value{
\describe{
  \item{count_matrice}{ A data.frame of the normalized count matrix.}
 }
}
\description{
This function normalizes the count matrix using loess
regression.
}
\examples{
nb_row <- 120
chromosome <- rep(1, nb_row)
index1 <- sample(1:100, nb_row, replace = TRUE)
index2 <- sample(1:100, nb_row, replace = TRUE)

m <- data.frame("mat_1" =  sample(1:500, nb_row, replace = TRUE),
                "mat_2" =  sample(1:500, nb_row, replace = TRUE),
                "mat_3" =  sample(1:500, nb_row, replace = TRUE),
                "mat_4" =  sample(1:500, nb_row, replace = TRUE))

mat <- cbind(chromosome, index1, index2, m)

normalizeCount(mat)

}
