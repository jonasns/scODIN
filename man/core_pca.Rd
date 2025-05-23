% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_pca.R
\name{core_pca}
\alias{core_pca}
\title{Wrapper for Seurat functions: FindVariableFeatures, ScaleData, and RunPCA}
\usage{
core_pca(
  seurat_obj,
  assay = "RNA",
  selection_method = "vst",
  nfeatures = 2000,
  verbose = TRUE,
  npcs = 50,
  reduction_name = "core_PCA"
)
}
\arguments{
\item{seurat_obj}{A Seurat object on which the functions will be applied.}

\item{assay}{Assay to use (default is \code{"RNA"}).}

\item{selection_method}{Method for selecting variable features in FindVariableFeatures (default is \code{"vst"}).}

\item{nfeatures}{Number of variable features to select in FindVariableFeatures (default is 2000).}

\item{verbose}{Logical indicating whether to print progress (default is TRUE).}

\item{npcs}{Number of principal components to compute in RunPCA (default is 50).}

\item{reduction_name}{Name of the dimensional reduction set to store (default is \code{"core_PCA"}).}
}
\value{
A Seurat object with variable features found, data scaled, and PCA run.
}
\description{
This function sequentially applies FindVariableFeatures, ScaleData, and RunPCA on a given Seurat object,
with user-defined input for parameters such as the number of variable features, number of principal components, etc.
}
\examples{
\dontrun{
seurat_obj_processed <- core_pca(seurat_obj, assay = "RNA", selection_method = "vst", nfeatures = 2000, npcs = 50, verbose = TRUE, reduction_name = "core_PCA")
}

}
