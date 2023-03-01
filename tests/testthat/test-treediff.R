# 4 clusters
leaves <- c(100,120,50,80)

# Generate 4 trees from trees1 and  trees from trees2 for each number of leaves
trees <- sapply(leaves, FUN = function(leaf){
  base_data <- matrix(rnorm(2000), nrow = leaf, ncol = 200)

  set1 <- replicate(4, sample(1:100, 50, replace = FALSE))
  set2 <- replicate(6, sample(101:200, 50, replace = FALSE))

  # Compute hierarchical clustering for each tree in set1
  trees1 <- apply(set1, 2, function(asample) {
    samples <- base_data[, asample]
    out <- hclust(dist(samples), method = "ward.D2")
    return(out)
  })

  # Compute hierarchical clustering for each tree in set2
  trees2 <- apply(set2, 2, function(asample) {
    samples <- base_data[, asample]
    out <- hclust(dist(samples), method = "ward.D2")
    return(out)
  })
  return(list("trees1" = trees1, "trees2" = trees2))
})

# Unlist the trees
trees1 <- unlist(trees[1,], recursive = FALSE)
trees2 <- unlist(trees[2,], recursive = FALSE)

# Set the number of replicates
replicates <- c(4, 6)

test_that("'treediff' works for simple cases", {

  # Call the function to be tested
  res <- treediff(trees1, trees2, replicates)

  # Check the output object has the expected names
  expect_named(res, c("method", "data.name", "p.value",
                      "statistic", "p.value.indiv"))

  # Check the output object has the expected class
  expect_s3_class(res, "treeTest")

  # Check that `p.value` is a numeric value
  expect_is(res$p.value, "numeric")

  # Check that `statistic` is a numeric value
  expect_is(res$statistic, "numeric")

  # Check that `p.value.indiv` is a numeric value
  expect_is(res$p.value.indiv, "numeric")

  # Check that `method` is a character string
  expect_is(res$method, "character")

  # Check that `data.name` is a character string
  expect_is(res$data.name, "character")

  # Check the length of `statistic` matches `p.value.indiv`
  expect_length(res$statistic, length(res$p.value.indiv))

  # Check the length of `p-value` matches the number of clusters
  expect_length(res$p.value, 4)

  # Check the length `statistic` matches the expected number of pairwise
  # comparisons between the clusters
  expect_length(res$statistic, sum(choose(leaves, 2)))
})

## Test errors
test_that("Test errors", {
  # Test if replicates is numeric vector
  expect_error(treediff(trees1 = trees1, trees2 = trees2, replicates = "abc"),
               "`replicates` is not a numeric vector")

  # Test if replicates is a vector of length 2
  expect_error(treediff(trees1 = trees1, trees2 = trees2, replicates = 1:3),
               "`replicates` must be a vector of length 2.")

  # Test if output is a list with class treeTest
  out <- treediff(trees1, trees2, c(4, 6))
  expect_true(inherits(out, "treeTest"))

  # Test if the number of leaves in each cluster is the same between the two
  # sets of trees
  trees2[[5]]$order <- c(trees2[[5]]$order, 101)
  expect_error(treediff(trees1, trees2, c(4, 6)),
               "the number of leaves in one or more clusters is different between the two sets of trees.")
})


# compute_squeeze
test_that("Test for compute_squeeze", {

  # Create a test data matrix
  test_matrix <- matrix(rnorm(24500), nrow = 1225, ncol = 20)
  replicates <- c(10,10)

  # Call the function
  res <- compute_squeeze(test_matrix, replicates)

  # Compare the result with the expected value
  expect_named(res, c("average_coph", "squeezed_var"))
  expect_is(res, "list")
  expect_is(res$average_coph, "data.frame")
  expect_is(res$squeezed_var, "list")
  expect_is(res$squeezed_var$var.post, "numeric")
})

# compute_pvalue
test_that("Test for compute_pvalue", {

  # Set up example inputs
  test_matrix <- matrix(rnorm(24500), nrow = 1225, ncol = 20)
  replicates <- c(10,10)

  # Calculate the average and variance of the test matrix
  res <- compute_squeeze(test_matrix, replicates)

  # Call the function
  pvals <- compute_pvalue(res$average_coph, res$squeezed_var, replicates)

  expect_named(pvals, c("statistics", "p.value", "cluster", "p"))

  # Check that the output is a data frame
  expect_is(pvals, "data.frame")
  expect_length(pvals$p.value, 1225)
  expect_length(pvals, 4)
})

# Test for the scale argument
test_that("Test for the scale argument", {

  # Perform the treediff test with scaling
  res <- treediff(trees1, trees2, replicates, scale = TRUE)

  # Check the output object has the expected names
  expect_named(res, c("method", "data.name", "p.value", "statistic", 
                      "p.value.indiv"))

  # Perform the treediff test without scaling
  result1 <- treediff(trees1, trees2, replicates)
  result2 <- treediff(trees1, trees2, replicates, scale = FALSE)
  expect_error({
    result3 <- treediff(trees1, trees2, replicates, scale = 5)
  }, "'scale' must be logical")

  expect_equal(result1, result2)
})

# Test for the order_labels argument
test_that("Test for the order_labels argument", {

  # Perform the treediff test with ordering the labels without labels
  expect_error(treediff(trees1, trees2, replicates, order_labels = TRUE))

  # Create a set of trees with different leaf orders
  trees1 <- list(hclust(dist(mtcars[, 1:2]), method = "ward.D2"),
                 hclust(dist(mtcars[, 3:4]), method = "ward.D2"))
  trees2 <- list(hclust(dist(mtcars[, 5:6]), method = "ward.D2"),
                 hclust(dist(mtcars[, 7:8]), method = "ward.D2"))

  # Perform the treediff test with and without ordering the labels
  res1 <- treediff(trees1, trees2, c(2, 2), order_labels = TRUE)
  res2 <- treediff(trees1, trees2, c(2, 2), order_labels = FALSE)
  expect_error({
    res3 <- treediff(trees1, trees2, c(2, 2), order_labels = 5)
  }, "'order_labels' must be logical")

  # Test that the p-values are the same for both tests
  expect_equal(res1, res2)
})


