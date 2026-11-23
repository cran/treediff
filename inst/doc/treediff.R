## ----setup, include = FALSE---------------------------------------------------
hicream_available <- requireNamespace("hicream", quietly = TRUE)
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", message = FALSE,
                      warning = FALSE, eval = hicream_available)

## ----nodata, echo = FALSE, eval = !hicream_available, results = "asis"--------
# cat("**The package `hicream` is not available on this machine: the code of",
#     "this vignette is displayed but has not been evaluated.** `hicream`",
#     "provides the Hi-C example files used below.")

## ----library------------------------------------------------------------------
library("treediff")

## ----files--------------------------------------------------------------------
replicates <- 1:3
conditions <- c("90", "110")
chromosomes <- 1:2

all_begins <- as.character(interaction(expand.grid(replicates, conditions),
                                       sep = "-"))
allMat <- lapply(chromosomes, function(chr) {
  paste0("Rep", all_begins, "-chr", chr, "_200000.bed")
})

files <- system.file("extdata", unlist(allMat), package = "hicream")
index <- system.file("extdata", "index.200000.longest18chr.abs.bed",
                     package = "hicream")
format <- rep("HiC-Pro", length(files))
binsize <- 200000
nb_replicates <- c(3, 3)
basename(files)

## ----hic2tree-----------------------------------------------------------------
res <- HiC2Tree(files, format, binsize, index, chromosomes, nb_replicates)
names(res)

## ----metadata-----------------------------------------------------------------
res$metadata

## ----testres------------------------------------------------------------------
res$testRes

## ----summary------------------------------------------------------------------
summary(res$testRes)

## ----hicdataset---------------------------------------------------------------
data <- HiCDataSet(files, format, binsize, chromosomes, index)
names(data)
length(data$HiCDataSet)
head(data$indexData)
data$index_mat_chr

## ----normalize----------------------------------------------------------------
norm_counts <- normalizeCount(data$HiCDataSet)
length(norm_counts)
head(norm_counts[[1]])

## ----clustertree--------------------------------------------------------------
clust <- clusterTree(norm_counts)
names(clust)
class(clust$trees[[1]])
clust$metadata

## ----treediff-----------------------------------------------------------------
nb_clusters <- nrow(unique(clust$metadata[, c("chromosome", "cluster")]))
in_cond1 <- rep(rep(c(TRUE, FALSE), nb_replicates), nb_clusters)
trees1 <- clust$trees[in_cond1]
trees2 <- clust$trees[!in_cond1]

test <- treediff(trees1, trees2, replicates = nb_replicates)
test$p.value

## ----identical----------------------------------------------------------------
all.equal(test$p.value, res$testRes$p.value)

## ----scale--------------------------------------------------------------------
treediff(trees1, trees2, replicates = nb_replicates, scale = TRUE)$p.value

## ----hist, fig.width = 6, fig.height = 4--------------------------------------
hist(test$p.value.indiv, breaks = 50, main = "Individual p-values",
     xlab = "p-value")

## ----coords-------------------------------------------------------------------
aggregate(cbind(start, end) ~ seqnames, data = res$index,
          FUN = function(x) c(min = min(x), max = max(x)))

## ----sessionInfo--------------------------------------------------------------
sessionInfo()


