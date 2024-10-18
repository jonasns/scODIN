#' Merge core cell type Seurat object and kNN predicted Seurat object
#'
#' This function merges a core Seurat object with a kNN Seurat object,
#' adds a source column to indicate the origin of each cell, and corrects 
#' the final labels by replacing dashes with underscores based on 
#' the original labels.
#'
#' @param core_object A Seurat object containing core cell types with original labels.
#' @param knn_object A Seurat object containing kNN results with predicted labels.
#' 
#' @return A merged Seurat object with updated final labels and a source column indicating 
#'         whether the cell originated from the core or the kNN object.
#' 
#' @examples
#' \dontrun{
#' merged_seurat <- odin_merge(core_object = pbmc_CD4_core, 
#'                                            knn_object = pbmc_CD4_knn)
#' }
#' 
#' @import stringr
#' 
#' @export
odin_merge <- function(core_object, knn_object) {
  
  # below contains the original cell types in final_labels
  core_object = core_object[, !core_object$final_labels %in% c("unknown")]
  
  # Ensure the final labels from knn_object are assigned
  knn_object$final_labels <- knn_object$predicted_id_knn
  
  # Merge the Seurat objects
  combined_seurat <- merge(
    x = core_object,
    y = list(knn_object),
    add.cell.ids = c("core", "knn"),
    project = "CombinedProject"  # Set default project name directly here
  )
  
  # Create a source column to indicate the origin of cells
  combined_seurat$source <- ifelse(is.na(combined_seurat$predicted_id_knn), "core", "knn")
  
  # Create a vector of new labels
  new_labels <- combined_seurat$final_labels
  
  # Create a vector of original labels from the core object
  original_labels <- core_object$final_labels
  
  # Function to replace dashes with underscores if corresponding original label exists
  replace_dashes <- function(new_label, original_labels) {
    # Check if new label contains dashes
    if (str_detect(new_label, "-")) {
      # Replace dashes with underscores
      potential_label <- str_replace_all(new_label, "-", "_")
      # Check if the potential label exists in the original labels
      if (potential_label %in% original_labels) {
        return(potential_label)
      }
    }
    return(new_label)
  }
  
  # Apply the function to each new label
  revised_labels <- sapply(new_labels, replace_dashes, original_labels = original_labels)
  
  # Update the new labels with the revised ones
  combined_seurat$final_labels <- revised_labels
  
  return(combined_seurat)
}