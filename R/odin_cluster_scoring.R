#' Combine cell number and scODIN score for labeling cells on the cluster level
#'
#' This function calculates a combined score based on cell number and scODIN score per cluster, 
#' and assigns the top scoring label to each cluster in the given Seurat object.
#' 
#' @param scData A Seurat object containing the clustering data and scODIN scores.
#' @param clustering_column The name of the column that contains the clustering information (default is "seurat_clusters").
#' 
#' @return The Seurat object with an updated `odin_classification` column, 
#'         containing the highest scoring label for each cluster.
#' 
#' @examples
#' scData <- odin_cluster_scoring(scData, "seurat_clusters")
#' 
#' @import dplyr
#' 
#' @export
odin_cluster_scoring <- function(scData, clustering_column = "seurat_clusters") {
  
  # Create a data frame with counts of cell types per cluster, excluding "unknown"
  cell_number_summary <- scData@meta.data %>%
    mutate(final_labels = scData$final_labels) %>%
    filter(final_labels != "unknown") %>%
    group_by(!!sym(clustering_column), final_labels) %>%
    summarize(ncells = n(), .groups = 'drop') %>%
    arrange(!!sym(clustering_column), desc(ncells))
  
  # Calculate the top 10 ODIN scores for each cluster
  odin_score_summary <- do.call(rbind, lapply(unique(scData@meta.data[[clustering_column]]), function(cluster_id) {
    # Identify the cells in the current cluster
    cluster_cells <- scData@meta.data[scData@meta.data[[clustering_column]] == cluster_id, ]
    
    # Calculate the scores for these cells
    scores_per_cell <- rowSums(t(odin_score_all)[, rownames(cluster_cells)])
    
    # Sort scores in decreasing order
    sorted_scores <- sort(scores_per_cell, decreasing = TRUE)
    
    # Create a data frame with the top scores
    top_scores <- head(data.frame(
      cluster = cluster_id,
      type = names(sorted_scores),
      scores = sorted_scores,
      ncells = nrow(cluster_cells)
    ), 10)
    
    return(top_scores)
  }))
  
  # Merge the datasets based on the corresponding columns
  merged_data <- cell_number_summary %>%
    inner_join(odin_score_summary, by = c(setNames("cluster", clustering_column), setNames("type", "final_labels")))
  
  # Multiply the scores by sqrt(ncells)
  merged_data <- merged_data %>%
    mutate(result = scores * sqrt(ncells.x))
  
  # Calculate enrichment ratios for each cluster
  result_data <- merged_data %>%
    group_by(!!sym(clustering_column)) %>%
    mutate(enrichment = {
      # Find the two highest values in the 'result' column
      sorted_results <- sort(result, decreasing = TRUE)
      highest <- sorted_results[1]
      second_highest <- sorted_results[2]
      # Compute the ratio
      ratio <- highest / second_highest
      # Repeat the ratio for each row in the group
      ratio
    }) %>%
    ungroup()
  
  # Select the top score for each cluster
  odin_cluster_scores <- result_data %>%
    group_by(!!sym(clustering_column)) %>%
    top_n(n = 1, wt = result)
  
  # Initialize the odin_classification column in the Seurat object
  scData@meta.data$odin_classification <- ""
  
  # Loop over each unique cluster and assign the highest scoring label
  for (j in unique(scData@meta.data[[clustering_column]])) {
    # Subset the cluster from odin_cluster_scores
    cluster_label_subset <- odin_cluster_scores[odin_cluster_scores[[clustering_column]] == j, ]
    # Assign the cluster name
    scData@meta.data$odin_classification[scData@meta.data[[clustering_column]] == j] <- cluster_label_subset$final_labels[1]
  }
  
  # Return the updated Seurat object
  return(scData)
}
