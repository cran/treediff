replicat <- 1:3
cond <- c("90", "110")

all_begins <- interaction(expand.grid(replicat, cond), sep = "-")
all_begins <- as.character(all_begins)

nb_chr <- 2
chromosomes <- 1:nb_chr

all_mat_chr <- lapply(chromosomes, function(chr){
  all_mat <- lapply(all_begins, function(ab) {
    mat_file <- paste0("Rep", ab, "-chr", chr, "_200000.bed")
    return(mat_file)
  })
  all_mat <- unlist(all_mat)
  return(all_mat)
})

files <- system.file("extdata", unlist(all_mat_chr), package = "treediff")
format <- rep(rep("HiC-Pro", 6), nb_chr)
binsize <- 200000
index <- system.file("extdata", "index.200000.longest18chr.abs.bed",
                     package = "treediff")
replicates <- c(3,3)

res <- HiC2Tree(files, format, binsize, index, chromosomes, replicates)

test_that("'HiC2Tree' works for simple cases", {

  # Check the output object has the expected names
  expect_named(res, c("trees", "metadata", "index", "testRes"))

  # Check that `trees` is a list
  expect_is(res$trees, "list")

  # Check that the first `trees` is hclust
  expect_s3_class(res$trees[[1]], "hclust")

  # Check that `metadata` is a data.frame
  expect_is(res$metadata, "data.frame")

  # Check that `index` is a data.frame
  expect_is(res$index, "data.frame")

  # Check that `testRes` is a treeTest
  expect_s3_class(res$testRes, "treeTest")

  # Check the length of `trees` matches `metadata`
  expect_length(res$trees, nrow(res$metadata))
})

## Test errors
test_that("Test errors", {

  files_test <- NULL
  # Test if `files` and `format` is correctly filled in
  expect_error(HiC2Tree(files_test, format, binsize, index, chromosomes,
  replicates),"`files`, `chromosomes` or `replicates` is incorrectly filled in.")

  binsize_test <- "test"
  # Test if `binzise` is a numeric vector
  expect_error(HiC2Tree(files, format, binsize_test, index, chromosomes,
                        replicates), "`binsize` is not a numeric vector.")

  format_test <- format
  format_test[4] <- "test"

  # Test if `format` is a correctly filled in
  expect_error(HiC2Tree(files, format_test, binsize, index, chromosomes,
                  replicates), "The `format` vector is incorrectly filled in.")

  index_test <- NULL
  # Test if `index` is not empty
  expect_error(HiC2Tree(files, format, binsize, index_test, chromosomes,
                        replicates), "`index` is empty.")

  binsize_test <- NULL
  # Test if `binsize` is not empty
  expect_error(HiC2Tree(files, format, binsize_test, index, chromosomes,
                        replicates), "`binsize` is empty.")
})

data <- HiCDOCDataSet(files, format, binsize, chromosomes, index)
norm_count <- normalizeCount(data$HiCDOCDataSet)
resTree <- clusterTree(norm_count)
tail(resTree$metadata)

trees1 <- resTree$trees[which(resTree$metadata$mat == "mat_1" |
                              resTree$metadata$mat == "mat_2" |
                              resTree$metadata$mat == "mat_3" |
                              resTree$metadata$mat == "mat_7" |
                              resTree$metadata$mat == "mat_8" |
                              resTree$metadata$mat ==  "mat_9")]
trees2 <- resTree$trees[which(resTree$metadata$mat == "mat_4" |
                              resTree$metadata$mat == "mat_5" |
                              resTree$metadata$mat == "mat_6" |
                              resTree$metadata$mat == "mat_10" |
                              resTree$metadata$mat == "mat_11" |
                              resTree$metadata$mat == "mat_12")]

resTest <- treediff(trees1, trees2, replicates)

test_that("Comparison HiC2Tree and sub-functions works for simple cases", {

  # Check the output object has the expected names
  expect_named(data, c("HiCDOCDataSet", "indexData", "index_mat_chr"))

  # Check that the column names of the `norm_count` data frame are as expected
  expect_equal(colnames(norm_count[[1]]), c("chromosome", "index1", "index2",
                  "mat_1", "mat_2", "mat_3", "mat_4", "mat_5", "mat_6"))

  # Check that the length of `norm_count` is as expected
  expect_length(norm_count, nb_chr)

  # Check that `resTree` is a list object
  expect_is(resTree, "list")

  # Check that the `mat_1` element of the first element of `resTree` is an
  # "hclust" object
  expect_is(resTree$trees$tree1, "hclust")

  # Check that `resTest` is of class "treeTest"
  expect_s3_class(resTest, "treeTest")

  # Check that `res$testRes` is identical to `resTest`
  expect_identical(res$testRes, resTest)

  # Check that the length of  `trees1` and `trees2` is equal to the number of
  # rows in `res$metadata`
  expect_length(c(trees1, trees2), nrow(res$metadata))
})

test_that("HiCDOCDataSet function", {

  # Check the output object has the expected names
  expect_named(data, c("HiCDOCDataSet", "indexData", "index_mat_chr"))

  # Check the output object has the expected class
  expect_s3_class(data$HiCDOCDataSet, "HiCDOCDataSetList")
  expect_s3_class(data$indexData, "data.frame")
  expect_s3_class(data$index_mat_chr, "data.frame")
})

test_that("normalizeCount function", {

  # Check the output object has the expected class
  expect_type(norm_count, "list")
  expect_length(norm_count, length(chromosomes))
})

test_that("clusterTree function", {

  # Check the output object has the expected class
  expect_type(resTree, "list")

  # Check the output object has the expected names
  expect_named(resTree, c("trees", "metadata"))
  expect_length(resTree$trees, nrow(resTree$metadata))
})

test_that("create_cluster function", {

  res <- create_cluster(norm_count[[1]])
  # Check the output object has the expected class
  expect_type(res, "list")
})

test_that("create_trees function", {

  cluster <- create_cluster(norm_count[[1]])

  # Get the indices of the rows and columns in the cluster
  selected <- rownames(cluster)[cluster == 1]

  # Select the sub-matrix
  red_mat <- ((norm_count[[1]]$index1 %in% selected) &
                (norm_count[[1]]$index2 %in% selected))
  red_mat <- norm_count[[1]][red_mat, ]

  res <- create_trees(red_mat, selected)
  # Check the output object has the expected class
  expect_type(res, "list")
  expect_s3_class(res$mat_1, "hclust")

  expect_length(selected, length(res$mat_1$order))
})
