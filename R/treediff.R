#' @title Treediff
#'
#' @description Perform the treediff test to compare two sets of trees.
#'
#' @aliases summary.treeTest
#' @aliases print.treeTest
#'
#' @details This function compares two sets of trees using a p-value aggregation
#' method. The p-values are obtained by the treediff method, as described in
#' (Neuvial \emph{et al.}, 2023).
#'
#' @param trees1 A list of trees corresponding to the first condition (set).
#' Trees are structured into groups (or clusters) with the same number of
#' replicates in each group. Trees are ordered by groups and then by replicates:
#' \{group1+rep1, group1+rep2, ...\}. One test is performed for each group.
#' @param trees2 A list of trees corresponding to the second condition. Trees
#' are also structured in groups (or clusters) that are exactly the same than
#' for the first condition. The number of replicates in each group can be
#' different from that of \code{trees1}.
#' @param replicates A numeric vector of length 2 with the number of replicates
#' for each condition.
#' @param scale Logical. If \code{TRUE}, the trees are all rescaled to have a
#' minimum height equal to 0 and a maximum height equal to 1. Default to
#' \code{FALSE}.
#' @param order_labels Logical. If \code{TRUE}, align leaves ordering in all
#' trees (required if your trees don't have their leaves ordered identically).
#' Default to \code{FALSE}.
#'
#' @return An object of class \code{treeTest} with the following entries:
#' \itemize{
#'   \item{p.value}{ the p-value for the treediff test.}
#'   \item{statistic}{ the value of the Student's statistic of each leaf pair of
#'   the tree test.}
#'   \item{p.value.indiv}{ the p-value of the Student's test for each leaf pair
#'   of the tree test.}
#'   \item{method}{ a character string indicating what type of test was
#'   performed.}
#'   \item{data.name}{ a character string giving the names of the tree
#'   conditions.}
#' }
#'
#' @author Gwendaëlle Cardenas\cr
#' Marie Chavent \email{marie.chavent@u-bordeaux.fr}\cr
#' Sylvain Foissac \email{sylvain.foissac@inrae.fr}\cr
#' Pierre Neuvial \email{pierre.neuvial@math.univ-toulouse.fr}\cr
#' Nathanaël Randriamihamison\cr
#' Nathalie Vialaneix \email{nathalie.vialaneix@inrae.fr}
#'
#' @references Neuvial Pierre, Randriamihamison Nathanaël, Chavent Marie,
#' Foissac Sylvain and Vialaneix Nathalie (2023) Testing differences in
#' structure between families of trees. \emph{Preprint submitted for
#' publication}.
#'
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom limma squeezeVar
#' @importFrom reshape2 colsplit
#' @importFrom stats cophenetic
#' @importFrom stats pt
#' @importFrom rlang .data
#' @import testthat
#'
#' @examples
#' leaves <- c(100, 120, 50, 80)
#'
#' trees <- sapply(leaves, FUN = function(leaf) {
#'   base_data <- matrix(rnorm(2000), nrow = leaf, ncol = 200)
#'
#' ## generates two sets of trees with 4 clusters with 100, 120, 50 and 80
#' ## leaves respectively
#' ## 4 replicates in the first condition and 6 in the second condition
#'
#'   set1 <- replicate(4, sample(1:100, 50, replace = FALSE))
#'   set2 <- replicate(6, sample(101:200, 50, replace = FALSE))
#'
#'   trees1 <- apply(set1, 2, function(asample) {
#'     samples <- base_data[, asample]
#'     out <- hclust(dist(samples), method = "ward.D2")
#'     return(out)
#'   })
#'
#'   trees2 <- apply(set2, 2, function(asample) {
#'     samples <- base_data[, asample]
#'     out <- hclust(dist(samples), method = "ward.D2")
#'     return(out)
#'   })
#'   return(list("trees1" = trees1, "trees2" = trees2))
#' })
#'
#' trees1 <- unlist(trees[1, ], recursive = FALSE)
#' trees2 <- unlist(trees[2, ], recursive = FALSE)
#' replicates <- c(4, 6)
#'
#' tree_pvals <- treediff(trees1, trees2, replicates)
#' ## 4 p-values, one for each cluster
#' tree_pvals$p.value

treediff <- function(trees1, trees2, replicates, scale = FALSE,
                     order_labels = FALSE) {

  # Check if `replicates` is numeric vector
  if (!is.numeric(replicates)) {
    stop("`replicates` is not a numeric vector")
  }

  # Check if the length of replicates is 2
  if (length(replicates) != 2) {
    stop("`replicates` must be a vector of length 2.")
  }

  # Check if the number of clusters is equal for both conditions
  if (length(trees1) / replicates[1] != length(trees2) / replicates[2]) {
    stop("The number of clusters is different between conditions (or ",
         "`replicates` is not correct).")
  }

  # Check the number of leaves is the same for each cluster
  tree_order1 <- lapply(trees1, "[[", "order")
  leaves1 <- sapply(tree_order1, length)

  tree_order2 <- lapply(trees2, "[[", "order")
  leaves2 <- sapply(tree_order2, length)

  if (!identical(unique(leaves1), unique(leaves2))) {
    stop("the number of leaves in one or more clusters is different between ",
    "the two sets of trees.")
  }
  
  # Check if 'scale' and 'order_labels' are logical
  if (!is.logical(scale)) stop("'scale' must be logical")
  if (!is.logical(order_labels)) stop("'order_labels' must be logical")

  # Merge trees from both conditions
  trees <- c(trees1, trees2)

  ## if order_labels is TRUE, take a reference tree and keep track of leave order
  if (order_labels) {
    ref_order <- trees1[[1]]$labels
    labels_perm <- lapply(trees, function(atree) match(ref_order, atree$labels))
  }

  # Compute cophenetic distance
  coph_dist <- sapply(trees, cophenetic, simplify = FALSE)

  # Normalize
  if (scale) {
    coph_dist <- normalize_trees(coph_dist)
  }

  # Convert cophenetic distances to vector
  if (order_labels) {
    coph_vect <- mapply(function(adist, anorder) {
      adist <- as.matrix(adist)
      adist <- adist[anorder, anorder]
      return(adist[upper.tri(adist)])
    }, coph_dist, labels_perm)

    # Convert the result to a list of vectors
    coph_vect <- lapply(1:ncol(coph_vect), function(acol) {
      return(as.vector(coph_vect[, acol]))
    })
  } else {
    coph_vect <- lapply(coph_dist, function(adist) {
      adist <- as.matrix(adist)
      return(adist[upper.tri(adist)])
    })
  }

  # Compute squeeze factor
  outs <- compute_squeeze(coph_vect, replicates)

  # Compute p-values
  outp <- compute_pvalue(outs$average_coph, outs$squeezed_var, replicates)

  # Aggregate p-values
  out_aggr <- suppressWarnings(outp %>%
    group_by(.data$cluster) %>%
    summarise("p.value" = min(sort(.data$p.value) / (1:.data$p)) * .data$p,
              .groups = "keep") %>%
    unique())

  # Store results in a list
  data_name <- paste(substitute(trees1), "and", substitute(trees2))
  out <- list("method" = "Tree test based on aggregated t-tests",
              "data.name" = data_name,
              "p.value" = out_aggr$p.value,
              "statistic" = outp$statistics,
              "p.value.indiv" = outp$p.value)

  # Assign class to the list
  class(out) <- "treeTest"

  # Return result
  return(out)
}

#' @export
#' @param x a \code{treeTest} object to print
#' @param ... not used
#' @rdname treediff

print.treeTest <- function(x, ...) {

  # print method
  cat("\n")
  cat(strwrap(x$method, prefix = "\t"), sep = "\n")
  cat("\n")

  # Print the name of the data used in the test
  cat("data: ", x$data.name, "\n", sep = "")

  # print alternative hypothesis
  cat("alternative hypothesis: the two sets of trees have different structure.")
  cat ("\n\n")

  # print the first 5 p-values
  print(colsplit(x$p.value, " ", "p-values:"), max = 5)
}

#' @method summary treeTest
#' @param object a \code{treeTest} object to print
#' @export
#' @rdname treediff

summary.treeTest <- function(object, ...) {

  # Print a header for the summary
  cat("\nSummary\n")

  # Print the class of the object
  cat("\tClass: ", class(object))
  cat("\n")

  # Print the test output
  print(object)
  cat("\n")

  # Print a summary of the p.value
  cat(" p-value summary:\n")
  print(summary(object$p.value))
  cat("\n")

  # Print a summary of the `statistic` and `p.value.indiv`
  summary(data.frame("indiv. statistics" = object$statistic,
                     "indiv. p-values" = object$p.value.indiv,
                     check.names = FALSE))
}

compute_squeeze <- function(dist_coph, replicates) {

  # Calculate number of clusters
  nb_cluster <- length(dist_coph) / sum(replicates)

  clusters1 <- rep(1:nb_cluster, each = replicates[1])
  clusters2 <- rep(1:nb_cluster, each = replicates[2])

  # Store replicates values for each group
  n1 <- replicates[1]
  n2 <- replicates[2]

  # Indices for each group
  set1 <- 1:length(clusters1)
  set2 <- 1:length(clusters2) + length(set1)

  # Average per groups and conditions
  average_coph_trees1 <- lapply(unique(clusters1), function(acluster) {
    where_clust <- which(clusters1 == acluster)
    colMeans(Reduce(rbind, dist_coph[where_clust]))
  })
  average_coph_trees2 <- lapply(unique(clusters1), function(acluster) {
    where_clust <- which(clusters2 == acluster) + length(clusters1)
    colMeans(Reduce(rbind, dist_coph[where_clust]))
  })

  # Merge average values and cluster vector
  cluster_length <- sapply(average_coph_trees1, length)
  average_coph <- data.frame("set1" = Reduce(c, average_coph_trees1),
                             "set2" = Reduce(c, average_coph_trees2),
                             "cluster" = rep(unique(clusters1), cluster_length))

  # Calculate variance
  sq_average_coph <- sweep(average_coph[-3]^2, 2, replicates, "*")

  # Sum of squared values for each group
  sum_sq_coph_trees1 <- lapply(unique(clusters1), function(acluster) {
    where_clust <- which(clusters1 == acluster)
    colSums(Reduce(rbind, dist_coph[where_clust])^2)
  })
  sum_sq_coph_trees2 <- lapply(unique(clusters1), function(acluster) {
    where_clust <- which(clusters2 == acluster) + length(clusters1)
    colSums(Reduce(rbind, dist_coph[where_clust])^2)
  })

  # `data.frame` with the sum of squared values for both groups and the cluster
  # number
  sum_sq_coph <- data.frame("set1" = Reduce(c, sum_sq_coph_trees1),
                            "set2" = Reduce(c, sum_sq_coph_trees2),
                            "cluster" = rep(unique(clusters1), cluster_length))

  # Compute the variance using the sum of squared values and the average values
  variances <- sum_sq_coph[, 1:2] - sq_average_coph
  variances <- rowSums(variances)
  variances <- as.vector(variances / (n1 + n2 - 2))

  # Calculate squeezed variance
  squeezed_var <- suppressWarnings(squeezeVar(variances, df = n1 + n2 - 2,
                                              robust = FALSE))

  # Return list of average_coph and squeezed_var
  return(list("average_coph" = average_coph, "squeezed_var" = squeezed_var))
}

compute_pvalue <- function(average_coph, squeezed_var, replicates) {

  # Extract the cluster information from the average_coph data frame
  cluster <- average_coph$cluster

  # Calculate the numerator for the t-statistic
  numerator <- Reduce("-", average_coph[, -3])

  # Store replicates values for each group
  n1 <- replicates[1]
  n2 <- replicates[2]

  # Calculate the denominator for the t-statistic
  denominator <- squeezed_var$var.post * (n1 + n2) / (n1 * n2)

  # Calculate the t-statistic
  statistics <- numerator / sqrt(denominator)

  # Extract the degrees of freedom
  d0 <- squeezed_var$df.prior

  # Calculate the p-value
  p.value <- 2*(1 - pt(abs(statistics), d0 + n1 + n2 - 2))

  # Create a vector of the number of pair for each cluster
  p = rep(table(cluster), table(cluster))

  # Create output data frame
  out <- data.frame ("statistics" = statistics,
                     "p.value" = p.value,
                     "cluster" = factor(cluster),
                     "p" = p)

  # Return the output data frame
  return(out)
}

normalize_trees <- function(coph_dist) {

  # Find the minimum value for each tree and store in all_minima
  all_minima <- sapply(coph_dist, min)

  # Subtract all_minima from each tree in coph_dist and store in out
  out <- mapply("-", coph_dist, all_minima, SIMPLIFY = FALSE)

  # Find the maximum value for each tree and store in all_maxima
  all_maxima <- sapply(out, max)

  # Divide each distance by all_maxima for each tree and store in out
  out <- mapply("/", out, all_maxima, SIMPLIFY = FALSE)

  # Return coph_dist normalized
  return(out)
}
