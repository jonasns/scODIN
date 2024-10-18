#' Simplify final double labels in scRNA-seq data
#'
#' This function simplifies the final double labels in a Seurat object
#' by replacing accepted double combinations of cell types with a simplified name provided by the user.
#'
#' @param scData A Seurat object that contains `final_labels` in its metadata.
#' @param accepted_doubles_table A tibble for accepted double labels. Usually coming from an Excel sheet.
#' @param cell_level Character string specifying the subset level (e.g., "CD4_T").
#' 
#' @return The Seurat object with updated `final_labels` in the metadata.
#'
#' @examples
#' \dontrun{
#' scData <- simplify_double_labels(scData, accepted_doubles_table, "CD4_T")
#' }
#' 
#' @export
simplify_double_labels <- function(scData, accepted_doubles_table, cell_level) {
  
  # Subset accepted_doubles_table to subset level
  accepted_doubles_table <- accepted_doubles_table[accepted_doubles_table$cell_level == cell_level, ]
  
  # Create a combined column in the accepted_doubles_table DataFrame
  accepted_doubles_table$combined <- paste(accepted_doubles_table$cell_type1, accepted_doubles_table$cell_type2, sep = "_")
  accepted_doubles_table$combined_rev <- paste(accepted_doubles_table$cell_type2, accepted_doubles_table$cell_type1, sep = "_")
  accepted_combinations <- c(accepted_doubles_table$combined, accepted_doubles_table$combined_rev)
  
  # Create a named vector for label replacements
  replacement_vector <- data.frame(accepted_combinations, accepted_doubles_table$simple_name)
  replacement_vector <- setNames(rep(accepted_doubles_table$simple_name, 2), accepted_combinations)
  
  # Replace labels in scData$final_labels
  scData$final_labels <- ifelse(scData$final_labels %in% names(replacement_vector),
                                replacement_vector[scData$final_labels],
                                scData$final_labels)
  
  return(scData)
}
