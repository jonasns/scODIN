% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/odin_merge.R
\name{odin_merge}
\alias{odin_merge}
\title{Merge core cell type Seurat object and kNN predicted Seurat object}
\usage{
odin_merge(core_object, knn_object)
}
\arguments{
\item{core_object}{A Seurat object containing core cell types with original labels.}

\item{knn_object}{A Seurat object containing kNN results with predicted labels.}
}
\value{
A merged Seurat object with updated final labels and a source column indicating
whether the cell originated from the core or the kNN object.
}
\description{
This function merges a core Seurat object with a kNN Seurat object,
adds a source column to indicate the origin of each cell, and corrects
the final labels by replacing dashes with underscores based on
the original labels.
}
\examples{
\dontrun{
merged_seurat <- odin_merge(core_object = pbmc_CD4_core, 
                                           knn_object = pbmc_CD4_knn)
}

}
