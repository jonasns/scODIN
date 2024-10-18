#' Wrapper for Seurat functions: FindVariableFeatures, ScaleData, and RunPCA
#'
#' This function sequentially applies FindVariableFeatures, ScaleData, and RunPCA on a given Seurat object, 
#' with user-defined input for parameters such as the number of variable features, number of principal components, etc.
#'
#' @param seurat_obj A Seurat object on which the functions will be applied.
#' @param assay Assay to use (default is `"RNA"`).
#' @param selection_method Method for selecting variable features in FindVariableFeatures (default is `"vst"`).
#' @param nfeatures Number of variable features to select in FindVariableFeatures (default is 2000).
#' @param verbose Logical indicating whether to print progress (default is TRUE).
#' @param npcs Number of principal components to compute in RunPCA (default is 50).
#' @param reduction_name Name of the dimensional reduction set to store (default is `"core_PCA"`).
#'
#' @return A Seurat object with variable features found, data scaled, and PCA run.
#' 
#' @examples
#' \dontrun{
#' seurat_obj_processed <- core_pca(seurat_obj, assay = "RNA", selection_method = "vst", nfeatures = 2000, npcs = 50, verbose = TRUE, reduction_name = "core_PCA")
#' }
#' 
#' @export
core_pca <- function(seurat_obj, 
                     assay = "RNA", 
                     selection_method = "vst", 
                     nfeatures = 2000, 
                     verbose = TRUE, 
                     npcs = 50, 
                     reduction_name = "core_PCA") {
  
  # Step 1: Find variable features
  seurat_obj <- FindVariableFeatures(seurat_obj, 
                                     selection.method = selection_method, 
                                     nfeatures = nfeatures, 
                                     verbose = verbose, 
                                     assay = assay)
  
  # Step 2: Scale the data
  seurat_obj <- ScaleData(seurat_obj, 
                          assay = assay, 
                          verbose = verbose)
  
  # Step 3: Run PCA
  seurat_obj <- RunPCA(seurat_obj, 
                       npcs = npcs, 
                       reduction.name = reduction_name, 
                       assay = assay, 
                       verbose = verbose)
  
  return(seurat_obj)
}