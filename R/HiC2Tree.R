#' @title Convert Hi-C to trees
#'
#' @description This function converts Hi-C data into trees, using 
#' \link[adjclust]{adjClust}. It takes as input a file path vector, the format 
#' of the input data, the bin size of the Hi-C array, the chromosomes to be 
#' included in the analysis, and the number of replicates. It returns a list 
#' containing all trees, metadata, index and treediff results.
#'
#' @param files A character vector containing the file paths of the input data.
#' @param format A character vector indicating the format of the input data:
#' "tabular", "cooler", "juicer", or "HiC-Pro".
#' @param binsize An integer indicating the bin size of the Hi-C matrix.
#' @param index A character indicating the path of the index for the input data.
#' Required (and used) only with the "HiC-Pro" format.
#' @param chromosomes A vector containing the chromosomes to be included in the
#' analysis.
#' @param replicates An integer indicating the number of replicates to be used
#' in \code{treediff}.
#' 
#' @return A list containing:
#' \item{trees}{ A list of all trees.}
#' \item{metadata}{ A data frame containing the following columns: names (name
#' of each tree), chromosome, cluster, and file.}
#' \item{index}{ A data table containing the correspondence of each bin in the
#' genome.}
#' \item{testRes}{ A list of treediff results for each cluster.}
#' 
#' @references
#' Christophe Ambroise, Alia Dehman, Pierre Neuvial, Guillem Rigaill, and 
#' Nathalie Vialaneix (2019) Adjacency-constrained hierarchical clustering of a
#' band similarity matrix with application to genomics. \emph{Algorithms for 
#' Molecular Biology}, \strong{14}(22), 363–389.
#'
#' @examples
#' replicates <- 1:3
#' cond <- c("90", "110")
#' all_begins <- interaction(expand.grid(replicates, cond), sep = "-")
#' all_begins <- as.character(all_begins)
#'
#' # single chromosome
#' nb_chr <- 1
#' chromosomes <- 1:nb_chr
#' all_mat_chr <- lapply(chromosomes, function(chr) {
#'   all_mat <- lapply(all_begins, function(ab) {
#'     mat_file <- paste0("Rep", ab, "-chr", chr, "_200000.bed")
#'   })
#'   all_mat <- unlist(all_mat)
#' })
#' index <- system.file("extdata", "index.200000.longest18chr.abs.bed",
#'                      package = "treediff")
#' format <- rep("HiC-Pro", length(replicates) * length(cond) * nb_chr)
#' binsize <- 200000
#' files <- system.file("extdata", unlist(all_mat_chr), package = "treediff")
#' replicates <- c(3, 3)
#' HiC2Tree(files, format, binsize, index, chromosomes, replicates)
#' 
#' \dontrun{
#' # two chromosomes
#' nb_chr <- 2
#' chromosomes <- 1:nb_chr
#' all_mat_chr <- lapply(chromosomes, function(chr) {
#'   all_mat <- lapply(all_begins, function(ab) {
#'     mat_file <- paste0("Rep", ab, "-chr", chr, "_200000.bed")
#'   })
#'   all_mat <- unlist(all_mat)
#' })
#' files <- system.file("extdata", unlist(all_mat_chr), package = "treediff")
#' format <- rep("HiC-Pro", length(replicates) * length(cond) * nb_chr)
#' replicates <- c(3, 3)
#' HiC2Tree(files, format, binsize, index, chromosomes, replicates)
#' }
#'
#' @export

HiC2Tree <- function(files, format, binsize = NULL, index = NULL,
                     chromosomes, replicates) {

  if (length(files) != length(chromosomes) * sum(replicates)){
    stop("`files`, `chromosomes` or `replicates` is incorrectly filled in.")
  }

  # create an HiCDOCDataSet object
  HiCDOCDataSet <- HiCDOCDataSet(files, format, binsize, chromosomes, index)

  # Extract the HiCDOCDataSet and index_mat_chr variables
  res <- HiCDOCDataSet$HiCDOCDataSet
  index_mat_chr <- HiCDOCDataSet$index_mat_chr

  # Normalize the count data using cyclic loess
  norm_matrices <- normalizeCount(res)

  res_clus <- clusterTree(norm_matrices)

  nb_cluster <- sapply(chromosomes, function(chr) {
    max(res_clus$metadata$cluster[which(res_clus$metadata$chromosome == chr)])
  })

  clus_sum <- sum(nb_cluster)
  T1 <- rep(rep(c(TRUE, FALSE), replicates), clus_sum)
  T2 <- rep(rep(c(FALSE, TRUE), replicates), clus_sum)

  trees1 <- res_clus$trees[T1]
  trees2 <- res_clus$trees[T2]

  testRes <- treediff(trees1, trees2, replicates)

  # Return a list containing the trees, metadata, and indexData
  return(list("trees" = res_clus$trees, "metadata" = res_clus$metadata,
              "index" = HiCDOCDataSet$indexData, "testRes" = testRes))
}

#' @title Create a HiCDOCDataSet object from a set of files
#'
#' @description This function creates a count matrix from a set of files
#' in different formats, such as tabular, cooler, juicer or HiC-Pro. It returns
#' a list of interaction matrices.
#'
#' @param files A character vector of file paths.
#' @param format A character vector of file formats corresponding to the files
#' in \code{file}. Supported formats are "tabular", "cooler", "juicer", and
#'   "HiC-Pro".
#' @param binsize An integer representing the bin size to use for cooler and
#'   juicer formats. Ignored for tabular and HiC-Pro formats.
#' @param chromosomes A character vector specifying the chromosomes to include
#'   in the output.
#' @param index A character vector of file paths to the index files required for
#'   HiC-Pro format. Ignored for other formats.
#'
#' @return A list containing the following objects:
#' \describe{
#'   \item{HiCDOCDataSet}{A list of interaction matrices of the HiCDOCDataSet
#'   class of the HiCDOC package, one for each file}
#'   \item{indexData}{A data frame of index data for each interaction in the
#'   matrices.}
#'   \item{index_mat_chr}{A data frame containing the name of the matrices and
#'   the corresponding chromosome.}
#' }
#'
#' @importFrom BiocGenerics cbind
#' @importFrom data.table as.data.table
#' @importFrom data.table setnames
#' @importFrom HiCDOC HiCDOCDataSetFromTabular
#' @importFrom HiCDOC HiCDOCDataSetFromCool
#' @importFrom HiCDOC HiCDOCDataSetFromHiC
#' @importFrom HiCDOC HiCDOCDataSetFromHiCPro
#' @importFrom InteractionSet interactions
#' @importFrom methods new
#' @importFrom purrr reduce
#' @importFrom SummarizedExperiment assay
#'
#' @examples
#' \dontrun{
#' replicates <- 1:2
#' cond <- "90"
#' all_begins <- interaction(expand.grid(replicates, cond), sep = "-")
#' all_begins <- as.character(all_begins)
#' 
#' nb_chr <- 2
#' chromosomes <- 1:nb_chr
#' all_mat_chr <- lapply(chromosomes, function(chr) {
#'   all_mat <- lapply(all_begins, function(ab) {
#'     mat_file <- paste0("Rep", ab, "-chr", chr, "_200000.bed")
#'   })
#'   all_mat <- unlist(all_mat)
#' })
#' index <- system.file("extdata", "index.200000.longest18chr.abs.bed",
#'                      package = "treediff")
#' format <- rep("HiC-Pro", length(replicates) * length(cond) * nb_chr)
#' binsize <- 200000
#' files <- system.file("extdata", unlist(all_mat_chr), package = "treediff")
#' HiCDOCDataSet(files, format, binsize, chromosomes, index)
#' }
#'
#' @export

HiCDOCDataSet <- function(files, format, binsize = NULL, chromosomes,
                          index = NULL) {

  if (is.null(chromosomes))stop("`chromsomes` is empty.")

  if (length(files) != length(format)) {
    stop("`files` or `format` is incorrectly filled in.")
  }

  lapply(unique(format), function(var) {
    if (var != "tabular" & var != "cooler" & var != "juicer" & var !="HiC-Pro"){
      stop("The `format` vector is incorrectly filled in.")
    }
  })

  lapply(unique(format), function(var) {
    if (var == "HiC-Pro" & is.null(index)) stop("`index` is empty.")
  })

  lapply(unique(format), function(var) {
    if (var != "tabular" & is.null(binsize)) stop("`binsize` is empty.")
    else {
      if (!is.null(binsize)){
        if (!is.numeric(binsize)) stop("`binsize` is not a numeric vector.")
      }
    }
  })

  # Set the conditions to "mat"
  conditions <- "mat"

  # Create a list of HiCDOCDataSet objects
  HiCDOCDataSet <- sapply(seq_along(files), function(i) {
    if (format[i] == "tabular") {
      HiCDOCDataSet <- HiCDOCDataSetFromTabular(files[i])
    } else if (format[i] == "cooler") {
      HiCDOCDataSet <- HiCDOCDataSetFromCool(files[i], i, conditions,
                                             binsize)
    } else if (format[i] == "juicer") {
      HiCDOCDataSet <- HiCDOCDataSetFromHiC(file[i], i, conditions,
                                            binsize)

    } else if (format[i] == "HiC-Pro") {
      HiCDOCDataSet <- HiCDOCDataSetFromHiCPro(files[i], bedPaths = index, i,
                                               conditions)
    }
    return(HiCDOCDataSet)
  })

  if (length(HiCDOCDataSet) > 1) {
    class(HiCDOCDataSet) <- "HiCDOCDataSetList"
  } else {
    HiCDOCDataSet <- HiCDOCDataSet[[1]]
  }

  # Create a matrix of chromosomes and matrix
  index_mat_chr <- sapply(HiCDOCDataSet, function(data) {
    chromosome <- data@chromosomes
    variables <- paste(data@colData$condition, data@colData$replicate,
                       sep = "_")
    return(list("chromosome" = chromosome, "variables" = variables ))
  })

  index_mat_chr <- as.data.frame(t(index_mat_chr))
  index_mat_chr$chromosome <- unlist(index_mat_chr$chromosome)
  index_mat_chr$variables <- unlist(index_mat_chr$variables)

  indexData <- lapply(chromosomes, function(chr) {

    dataSet <- HiCDOCDataSet[index_mat_chr$chromosome == chr]

    indexData <- lapply(dataSet, function(data) {
      index <- as.data.frame(data@interactions@regions)
    })

    indexData <- Reduce(rbind, indexData)

    indexData <- unique(indexData)
    indexData <- indexData[order(indexData$index),]
  })

  indexData <- Reduce(rbind, indexData)

  return(list("HiCDOCDataSet" = HiCDOCDataSet, "indexData" = indexData,
              "index_mat_chr" = index_mat_chr))
}

#' @title Normalize count matrix using cyclic loess
#'
#' @description This function normalizes the count matrix using loess
#' regression.
#'
#' @param count_matrice The count matrix to normalize.
#'
#' @return \describe{
#'   \item{count_matrice}{ A data.frame of the normalized count matrix.}
#'  }
#'
#' @importFrom csaw normOffsets
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @export
#'
#' @examples
#' nb_row <- 120
#' chromosome <- rep(1, nb_row)
#' index1 <- sample(1:100, nb_row, replace = TRUE)
#' index2 <- sample(1:100, nb_row, replace = TRUE)
#'
#' m <- data.frame("mat_1" =  sample(1:500, nb_row, replace = TRUE),
#'                 "mat_2" =  sample(1:500, nb_row, replace = TRUE),
#'                 "mat_3" =  sample(1:500, nb_row, replace = TRUE),
#'                 "mat_4" =  sample(1:500, nb_row, replace = TRUE))
#'
#' mat <- cbind(chromosome, index1, index2, m)
#'
#' normalizeCount(mat)
#' 
normalizeCount <- function(count_matrice){
  UseMethod("normalizeCount")
}

#' @export
normalizeCount.default <- function(count_matrice) {

  if (!identical(names(count_matrice[ ,c(1:3)]), c("chromosome", "index1",
                                                   "index2"))) {
    message("Check that the first 3 columns contain in order `chromsome`,
            `index1` and `index2`")
  }

  # Extract the count values as a matrix
  counts <- as.matrix(count_matrice[ ,-c(1:3)])
  colnames(counts) <- NULL

  # Compute the total number of reads for each sample
  cur_dge <- SummarizedExperiment(list(counts = counts))
  cur_dge$totals <- colSums(count_matrice[ ,-c(1:3)])

  lib.sizes <- cur_dge$totals

  # Normalize the counts
  offsets <- normOffsets(cur_dge, se.out = FALSE)

  offsets <- offsets - mean(log(lib.sizes))

  counts <- counts / exp(offsets)
  count_matrice <- as.data.frame(count_matrice)
  count_matrice[ ,-c(1:3)] <- data.frame(counts)

  return(count_matrice)
}

#' @export
normalizeCount.list <- function(count_matrice) {
  count_matrice <- lapply(count_matrice, normalizeCount.default)
  return(count_matrice)
}

#' @export
normalizeCount.HiCDOCDataSet <- function(count_matrice){

  chromosome <- count_matrice@chromosomes

  interaction <- interactions(count_matrice)
  interaction <- as.data.frame(interaction)
  interaction_mat <- cbind(interaction[, c(6,12)], assay(count_matrice))

  names(interaction_mat) <- c("index1", "index2", "mat_1")

  interaction_mat[is.na(interaction_mat)] <- 0
  chr <- rep(chromosome, nrow(interaction_mat))
  count_matrice <- cbind("chromosome" = chr, interaction_mat)

  return(count_matrice)
}

#' @export
normalizeCount.HiCDOCDataSetList <- function(count_matrice){

  matChr <- sapply(count_matrice, function(data) data@chromosomes)
  chromosomes <- unique(matChr)

  mat_count <- lapply(chromosomes, function(chr) {
    index <- matChr == chr
    dataSet <- count_matrice[index]

    data <- lapply(dataSet, function(data) {
      interaction <- interactions(data)
      interaction <- as.data.frame(interaction)
      interaction <- cbind(interaction[, c(6,12)], assay(data))
    })

    interaction_mat <- suppressWarnings(reduce(data, merge,
                                    by = c("index1", "index2"), all = TRUE))
    names_mat <- paste("mat", which(index), sep = "_")
    names(interaction_mat) <- c("index1", "index2", names_mat)

    interaction_mat[is.na(interaction_mat)] <- 0
    chr <- rep(chr, nrow(interaction_mat))
    interaction_mat <- cbind("chromosome" = chr, interaction_mat)

    return(interaction_mat)
  })

  mat_count <- normalizeCount.list(mat_count)
  return(mat_count)
}

#' @title Create hierarchical clustering trees for each cluster in a given
#' matrix
#'
#' @description This function creates a hierarchical clustering tree for each
#' group in a given matrix. The function breaks the chromosome into clusters
#' using the "broken stick" method and then converts the clusters into trees.
#'
#' @param mat A matrix containing the data to cluster. It should have columns
#' named 'index1', 'index2', 'chromosome and one column for each matrices.
#'
#' @return A list containing the following objects:
#' \describe{
#'   \item{trees}{ A list of hierarchical clustering trees, one for each cluster
#'   in the matrix.}
#'   \item{metadata}{ A data frame containing the following columns: names (name
#'   of each tree), chromosome, cluster, and file. }
#' }
#'
#' @importFrom adjclust adjClust
#' @importFrom adjclust select
#' @importFrom stats hclust
#'
#' @export
#'
#' @examples
#' n <- 5
#'
#' matrice <- matrix(runif(n*n), nrow = n, ncol = n)
#' matrice[lower.tri(matrice)] <- t(matrice)[lower.tri(matrice)]
#' matrice <- as.data.frame(matrice)
#' names(matrice) <- c("mat_1", "mat_2", "mat_3", "mat_4", "mat_5")
#'
#' chromosome <- rep(1, n)
#' index1 <- sample(1:100, n, replace = TRUE)
#' index2 <- sample(1:100, n, replace = TRUE)
#'
#' mat <- cbind(chromosome, index1, index2, matrice)
#'
#' res <- clusterTree(mat)
clusterTree <- function(mat){
  UseMethod("clusterTree")
}

#' @export
clusterTree.default <- function(mat) {

  if (!inherits(mat, "data.frame")){
    stop("`mat` is not a `data.frame`")
  }

  if (!identical(names(mat[ ,c(1:3)]), c("chromosome", "index1",
                                                   "index2"))) {
    message("Check that the first 3 columns contain in order `chromsome`,
            `index1` and `index2`")
  }

  # Create cluster for each chromosome
  cur_clust <- create_cluster(mat)

  # List of all clusters
  all_clusters <- unique(cur_clust)

  # Create a tree for each cluster
  trees <- sapply(seq_along(all_clusters$merged_clust), function(ac) {

    # Get the indices of the rows and columns in the cluster
    selected <- rownames(cur_clust)[cur_clust == ac]

    # Select the sub-matrix
    red_mat <- ((mat$index1 %in% selected) & (mat$index2 %in% selected))
    red_mat <- mat[red_mat, ]

    # Create trees
    trees <- create_trees(red_mat, selected)

    return(trees)
  }, simplify = FALSE)

  chromosome <- unique(mat$chromosome)
  nb_cluster <- length(trees)

  nb_mat_chr <- unique(sapply(trees, length))

  chromosomes <- rep(chromosome, nb_cluster * nb_mat_chr)
  clusters <- rep(1:nb_cluster, each = nb_mat_chr)

  trees <- unlist(trees, recursive = FALSE)
  names_mat <- names(trees)

  names <- paste0("tree", seq_along(trees))
  names(trees) <- names

  metadata <- data.frame("names" = names, "chromosome" = chromosomes,
                         "cluster" = clusters, "mat" = names_mat)

  return(list("trees" = trees, "metadata" = metadata))

}

#' @export
clusterTree.list <- function(mat){
  res <- lapply(mat, clusterTree.default)

  metadata <- lapply(res, function(data) data$metadata)
  metadata <- Reduce(rbind, metadata)

  trees <- lapply(res, function(data) data$trees)
  trees <- unlist(trees, recursive = FALSE)

  names <- paste0("tree", seq_along(trees))
  names(trees) <- names

  metadata$names <- names

  return(list("trees" = trees, "metadata" = metadata))
}

create_cluster <- function(res) {

  # extract all unique bin
  all_bins <- unique(c(res$index1, res$index2))

  # create index mapping for the bins
  res$bindex1 <- match(res$index1, all_bins)
  res$bindex2 <- match(res$index2, all_bins)

  nb_col <- length(res)
  merged_mat <- rowSums(res[ ,-c(1:3, nb_col -1, nb_col)])
  cur_mat <- matrix(0, ncol = length(all_bins), nrow = length(all_bins))
  cur_mat[cbind(res$bindex1, res$bindex2)] <- merged_mat
  cur_mat[cbind(res$bindex2, res$bindex1)] <- merged_mat

  # log transformation
  cur_mat <- log(cur_mat + 1)

  # perform constrained hierarchical clustering
  merged_res <- adjClust(cur_mat, type = "similarity", h = length(all_bins) - 1)

  merged_res$labels <- as.character(all_bins)

  # select the number of clusters with broken stick method
  merged_clust <- adjclust::select(merged_res, type = "bstick")
  names(merged_clust) <- all_bins

  return(data.frame(merged_clust))
}

create_trees <- function(red_mat, selected){

  # match bin indices in the reduced similarity matrix to selected indices
  match1 <- match(red_mat$index1, selected)
  match2 <- match(red_mat$index2, selected)

  # perform hierarchical clustering for each column
  adjclust_out <- sapply(red_mat[, -c(1:3)], function(arep) {

    mat_crep <- matrix(0, ncol = length(selected), nrow = length(selected))
    mat_crep[cbind(match1, match2)] <- arep
    mat_crep[cbind(match2, match1)] <- arep

    # perform hierarchical clustering on the similarity matrix
    out <- adjClust(log(mat_crep + 1), type = "similarity")
    class(out) <- "hclust"

    return(out)
  }, simplify = FALSE)

  return(adjclust_out)
}
