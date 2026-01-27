#' k-Nearest Neighbors (kNN) search using core cells as reference and unknown cells as quesry.
#'
#' This function performs a k-Nearest Neighbors (kNN) search on a Seurat object 
#' to classify undefined cells using a reference dataset of core cell types. It wraps around Seurat's 
#' `FindTransferAnchors` and `TransferData` functions with additional parameters.
#'
#' @param scData_sub A Seurat object subset containing the cells to analyze. For instance subsetted to CD4 T cells.
#' @param reference_obj A Seurat object containing the reference data (e.g., downsampled core cell types).
#' @param knn_method The kNN method to use, either "annoy" or "rann" (default is "annoy").
#' @param npcs Number of principal components to use (default is 50).
#' @param dims Dimensions to use for kNN search (default is `1:50`).
#' @param normalization_method The normalization method for the Seurat object (default is `"LogNormalize"`).
#' @param verbose Logical flag to indicate if progress should be printed (default is TRUE).
#' @param n_trees Number of trees to use for the kNN method (default is 50 for annoy).
#' @param k_weight The number of neighbors to use when predicting labels (default is 50).
#' @param k_filter Whether to filter anchors by distance (default is `NA`).
#' @param mapping_score_k Logical, whether to compute the mapping score (default is TRUE).
#' @param k.score Number of neighbors (k) used when scoring anchors during label transfer. Must be less than the number of query cells (default is 30).
#' 
#' @return A Seurat object with transferred labels and prediction scores stored in the "predictions" assay.
#' 
#' @examples
#' \dontrun{
#' transferred_obj <- odin_knn(
#'   scData_sub = scData_sub, 
#'   reference_obj = reference_obj, 
#'   knn_method = "annoy", 
#'   npcs = 50, 
#'   dims = 1:50)
#' }
#' 
#' @export
odin_knn <- function(scData_sub, 
                     reference_obj, 
                     knn_method = "annoy", 
                     npcs = 50, 
                     dims = 1:50, 
                     normalization_method = "LogNormalize", 
                     verbose = TRUE, 
                     n_trees = 50, 
                     k_weight = 50, 
                     k_filter = NA, 
                     mapping_score_k = TRUE,
                     k.score = 30) {
  
  # Step 1: Selecting undefined cells from the subset
  message("Selecting undefined cells")
  
  # Select "unknown" cells (the undefined cells)
  pbmc_undefined <- scData_sub[, scData_sub$final_labels %in% "unknown"]
  
  message("Commencing kNN search using ", knn_method)
  
  # Track the start time
  start_time <- Sys.time()
  
  # Step 2: Perform the FindTransferAnchors function
  anchors <- FindTransferAnchors(
    reference = reference_obj, 
    query = pbmc_undefined, 
    normalization.method = normalization_method, 
    reference.reduction = "core_PCA", 
    npcs = npcs, 
    dims = dims, 
    k.filter = k_filter,
    k.score = k.score,
    mapping.score.k = mapping_score_k, 
    nn.method = knn_method
  )
  
  # Step 3: Transfer the data using TransferData
  pbmc_undefined_transf <- TransferData(
    anchorset = anchors, 
    refdata = list(celltype.l1 = reference_obj$final_labels), 
    reference = reference_obj, 
    query = pbmc_undefined, 
    weight.reduction = "pcaproject", 
    l2.norm = FALSE, 
    dims = NULL, 
    k.weight = k_weight, 
    sd.weight = 1, 
    eps = 0, 
    n.trees = n_trees, 
    verbose = verbose, 
    slot = "data", 
    prediction.assay = TRUE, 
    only.weights = FALSE, 
    store.weights = TRUE
  )
  
  # Step 4: Store prediction scores in the Seurat object
  pbmc_undefined_transf[["predictions"]] <- pbmc_undefined_transf@assays$prediction.score.celltype.l1
  
  # Calculate the time elapsed
  end_time <- Sys.time()
  message("kNN calculation has finished in ", format(end_time - start_time))
  
  return(pbmc_undefined_transf)
}
