% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/odin_knn.R
\name{odin_knn}
\alias{odin_knn}
\title{k-Nearest Neighbors (kNN) search using core cells as reference and unknown cells as quesry.}
\usage{
odin_knn(
  scData_sub,
  reference_obj,
  knn_method = "annoy",
  npcs = 50,
  dims = 1:50,
  normalization_method = "LogNormalize",
  verbose = TRUE,
  n_trees = 50,
  k_weight = 50,
  k_filter = NA,
  mapping_score_k = TRUE
)
}
\arguments{
\item{scData_sub}{A Seurat object subset containing the cells to analyze. For instance subsetted to CD4 T cells.}

\item{reference_obj}{A Seurat object containing the reference data (e.g., downsampled core cell types).}

\item{knn_method}{The kNN method to use, either "annoy" or "rann" (default is "annoy").}

\item{npcs}{Number of principal components to use (default is 50).}

\item{dims}{Dimensions to use for kNN search (default is \code{1:50}).}

\item{normalization_method}{The normalization method for the Seurat object (default is \code{"LogNormalize"}).}

\item{verbose}{Logical flag to indicate if progress should be printed (default is TRUE).}

\item{n_trees}{Number of trees to use for the kNN method (default is 50 for annoy).}

\item{k_weight}{The number of neighbors to use when predicting labels (default is 50).}

\item{k_filter}{Whether to filter anchors by distance (default is \code{NA}).}

\item{mapping_score_k}{Logical, whether to compute the mapping score (default is TRUE).}
}
\value{
A Seurat object with transferred labels and prediction scores stored in the "predictions" assay.
}
\description{
This function performs a k-Nearest Neighbors (kNN) search on a Seurat object
to classify undefined cells using a reference dataset of core cell types. It wraps around Seurat's
\code{FindTransferAnchors} and \code{TransferData} functions with additional parameters.
}
\examples{
\dontrun{
transferred_obj <- odin_knn(
  scData_sub = scData_sub, 
  reference_obj = reference_obj, 
  knn_method = "annoy", 
  npcs = 50, 
  dims = 1:50)
}

}
