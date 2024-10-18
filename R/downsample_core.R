#' Downsample the identified core cells to a specified cell count per final label
#'
#' This function downsamples a Seurat object to have a similar number of core cells per final label.
#' It randomly selects cells up to a user-defined target count or the median cell count if specified.
#'
#' @param scData A Seurat object containing the single-cell data and the core cell types as `final_labels`.
#' @param target_count An integer specifying the target count for downsampling. It is overruled by use_median = TRUE.
#' @param use_median A logical indicating whether to use the median cell count (TRUE) or the target count (FALSE).
#' @param labels_to_exclude A character vector specifying the labels to exclude (default is `c("unknown")`).
#' @param seed An integer value for setting the seed to ensure reproducibility (default is 1984).
#'
#' @return A downsampled Seurat object with the specified cell level and a similar number of cells per final label.
#' 
#' @examples
#' downsampled_seurat <- downsample_core(scData, target_count = 100, use_median = FALSE, 
#'                                       labels_to_exclude = c("unknown"), seed = 1984)
#'
#'@export
downsample_core <- function(scData, target_count = NULL, use_median = TRUE, labels_to_exclude = c("unknown"), seed = 1984) {
  
  # Further subset to remove unwanted labels (including "unknown" by default)
  scData_sub <- scData[, !scData$final_labels %in% labels_to_exclude]
  
  # Get the number of cells per final label
  cell_counts <- table(scData_sub$final_labels)
  
  # Determine the target count based on user input
  if (use_median) {
    target_count <- median(cell_counts)
  }
  
  # Find the maximum cell count across all final labels
  max_cell_count <- max(cell_counts)
  
  # Check if the target count is greater than the maximum cell count
  if (!is.null(target_count) && target_count > max_cell_count) {
    subset_name <- names(cell_counts)[which.max(cell_counts)]
    message(sprintf("The downsampling input number is larger than any subset. Your highest cell count is %d for %s. Downsampling will not take place.", 
                    max_cell_count, subset_name))
    return(scData_sub)  # Return the original subset if no downsampling occurs
  }
  
  # Initialize a vector to hold the downsampled cells
  cells_to_keep <- c()
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Loop through each label and downsample if needed
  for (label in names(cell_counts)) {
    # Get cells for the current label
    label_cells <- colnames(scData_sub)[scData_sub$final_labels == label]
    
    # Downsample if the number of cells is greater than the target count
    if (length(label_cells) > target_count) {
      downsampled_cells <- sample(label_cells, target_count)
    } else {
      downsampled_cells <- label_cells
    }
    
    # Combine cells to keep
    cells_to_keep <- c(cells_to_keep, downsampled_cells)
  }
  
  # Subset the Seurat object with the downsampled cells
  scData_sub_downsampled <- subset(scData_sub, cells = cells_to_keep)
  
  return(scData_sub_downsampled)
}