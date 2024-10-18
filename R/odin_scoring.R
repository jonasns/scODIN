#' scODIN scoring function
#'
#' This function calculates scODIN scores for cell types based on priority genes from a provided gene priority table and a Seurat object with gene expression values.
#'
#' @param gene_priority_table A tibble containing gene priorities for different cell types. Usually provided in an Excel sheet.
#' @param scData A Seurat object containing normalized and scaled scRNA-seq gene expression data.
#' @param core_cell_cutoff A numeric value indicating the cutoff for what is considered a core cell type. Default is 5, but should be adjusted lower if using less deeply sequenced dataset.
#' @param similarity_threshold A numeric value defining the similarity threshold for determining double labels. Default is 1.
#' @param accepted_doubles_table A tibble for accepted double labels. Usually provided in an Excel sheet.
#' @param cell_level A string telling which level the analysis should be done on (e.g., "Top", "CD4_T").
#'
#' @return A Seurat object with scODIN scores and metadata added.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' gene_priority_table <- read_excel("~/Desktop/241007_TableS1_gene_priority_table.xlsx")
#' scData <- your_seurat_object # Load your Seurat object
#' result <- scODIN_scoring(gene_priority_table, scData, core_cell_cutoff = 5, similarity_threshold = 1,  accepted_doubles_table)
#'}
#'
#' @export
odin_scoring <- function(gene_priority_table, scData, core_cell_cutoff = 5, similarity_threshold = 1, accepted_doubles_table, cell_level = "Top") {
  
  start_time <- Sys.time()
  message("Preparing gene priority table for ", cell_level, " cell level")
  
  # Subset to the specified cell level
  gene_priority_table <- gene_priority_table[gene_priority_table$cell_level == cell_level, ]
  
  # Check if there are any missing genes and print a message if there are
  missing_genes <- gene_priority_table$gene_id[!gene_priority_table$gene_id %in% rownames(scData)]
  if (length(missing_genes) > 0) {
    message("The following priority genes are not found in the scRNAseq dataset: ", paste(missing_genes, collapse = ", "), ". Please check if your genes are labelled correctly.")
  } else {
    message("All priority genes are found in the scRNAseq dataset (~˘▾˘)~.")
  }
  
  message("Preparing gene expression table")
  
  # Select genes in the gene_priority_table found in the data
  gene_priority_table <- gene_priority_table[gene_priority_table$gene_id %in% rownames(scData), ]
  
  # Downsample Seurat object to same genes
  scData_ds <- scData[gene_priority_table$gene_id, ]
  
  # Select just the data that we need
  scData_ds <- scData_ds[["RNA"]]$scale.data
  
  # Extract cell types from the priority table
  unique_cell_types <- unique(gene_priority_table$cell_type)
  
  # Generate an empty matrix, same length as the data
  odin_score_all <- matrix(nrow = 0, ncol = dim(scData)[2])
  colnames(odin_score_all) <- colnames(scData)
  
  message("Initiating scODIN scoring")
  
  for (j in unique_cell_types) {
    # Select only genes present per cell type
    gene_priority_table_1 <- subset(gene_priority_table, cell_type == j)
    
    # Check for duplicates in gene_priority_table
    if (anyDuplicated(gene_priority_table_1$gene_id) > 0) {
      stop(" There are duplicates in gene_priority_table for cell type ", j, ". Please correct your gene_priority table and try again.")
    } else {
      message(" Now processing ", j)
    }
    
    # Downsample Seurat object to same genes
    scData_sub <- scData_ds[gene_priority_table_1$gene_id, ]
    
    if (dim(gene_priority_table_1)[1] > 1) {
      # Sort the two data objects
      scData_sub <- scData_sub[order(rownames(scData_sub)), ]
      gene_priority_table_1 <- gene_priority_table_1[order(gene_priority_table_1$gene_id), ]
    }
    
    # Make a new matrix to store the result
    odin_score_matrix <- matrix(nrow = length(gene_priority_table_1$gene_id), ncol = dim(scData)[2])
    
    # Multiply expression values with gene_priority and divide by the square-root of the number of genes per cell type
    if (length(gene_priority_table_1$gene_id) > 1) {
      for (i in 1:length(gene_priority_table_1$gene_id)) {
        odin_score_matrix[i, ] <- scData_sub[i, ] * gene_priority_table_1$gene_priority[i] / sqrt(length(gene_priority_table_1$gene_id))
      }
    } else {
      odin_score_matrix[1, ] <- scData_sub * gene_priority_table_1$gene_priority / sqrt(1)
    }
    
    # Convert the matrix to a data frame and add to combined dataframe
    odin_score_all <- rbind(odin_score_all, colSums(odin_score_matrix, na.rm = TRUE))
  }
  
  # Label the rownames according to the cell types
  rownames(odin_score_all) <- unique_cell_types
  
  message("Applying core cell type cutoff")
  
  # Apply core_cell_cutoff to all values in the dataframe
  odin_score_all <- ifelse(odin_score_all < core_cell_cutoff, 0, odin_score_all)
  
  # Transpose the df for easier handling of cell-wise operations
  odin_score_all <- t(odin_score_all)
  
  message("Naming cells based on tier")
  
  # Initialize columns for single labels, double labels, and unknowns
  single_labels <- rep("unknown", nrow(odin_score_all))
  double_labels <- rep("not_double_label", nrow(odin_score_all))
  
  # Create a list of cell types for each tier
  tier_list <- split(gene_priority_table$cell_type, gene_priority_table$tier)
  
  # Function to assign labels based on tier
  assign_labels <- function(scores, tier_cell_types, similarity_threshold) {
    cell_types <- names(which(scores > 0 & names(scores) %in% tier_cell_types))
    num_labels <- length(cell_types)
    
    if (num_labels == 1) {
      return(list(single = cell_types[1], double = "not_double_label"))
    } else if (num_labels == 2) {
      top_two_scores <- sort(scores[scores > 0], decreasing = TRUE)
      if ((top_two_scores[1] - top_two_scores[2]) < similarity_threshold) {
        return(list(single = "unknown", double = paste(cell_types, collapse = "_")))
      } else {
        return(list(single = names(which.max(scores)), double = "not_double_label"))
      }
    }
    return(list(single = "unknown", double = "not_double_label"))
  }
  
  # Initialize a vector to track labeled cells
  labeled_cells <- rep(FALSE, nrow(odin_score_all))
  
  for (tier in sort(unique(gene_priority_table$tier))) {
    tier_cell_types <- tier_list[[as.character(tier)]]
    
    for (i in 1:nrow(odin_score_all)) {
      if (!labeled_cells[i]) {
        scores <- odin_score_all[i, ]
        labels <- assign_labels(scores, tier_cell_types, similarity_threshold)
        if (labels$single != "unknown") {
          single_labels[i] <- labels$single
          labeled_cells[i] <- TRUE
        } else if (labels$double != "not_double_label") {
          double_labels[i] <- labels$double
          labeled_cells[i] <- TRUE
        }
      }
    }
  }
  
  message("Processing accepted double labels")
  
  # Create a combined column in the accepted_doubles_table
  accepted_doubles_table$combined <- paste(accepted_doubles_table$cell_type1, accepted_doubles_table$cell_type2, sep = "_")
  accepted_doubles_table$combined_rev <- paste(accepted_doubles_table$cell_type2, accepted_doubles_table$cell_type1, sep = "_")
  accepted_combinations <- c(accepted_doubles_table$combined, accepted_doubles_table$combined_rev)
  
  # Create a new column for combined labels
  combined_labels <- single_labels
  
  # Check each double label and see if it is in the accepted list
  for (i in 1:nrow(odin_score_all)) {
    if (double_labels[i] != "not_double_label") {
      if (double_labels[i] %in% accepted_combinations) {
        combined_labels[i] <- double_labels[i]
      }
    }
  }
  
  message("Adding scODIN information to Seurat object")
  
  # Add the odin.scores and new metadata to the Seurat object
  for (i in 1:length(unique_cell_types)) {
    scData <- AddMetaData(scData, odin_score_all[, i], col.name = unique_cell_types[i])
  }
  
  scData <- AddMetaData(scData, single_labels, col.name = "single_labels")
  scData <- AddMetaData(scData, double_labels, col.name = "double_labels")
  scData <- AddMetaData(scData, combined_labels, col.name = "final_labels")
  
  end_time <- Sys.time()
  message("scODIN calculation has finished in ", format(end_time - start_time))
  
  # place the octim_score_all matrix in the global environment to use in cluster ID function 
  assign("odin_score_all", odin_score_all, envir = .GlobalEnv)
  
  
  return(scData)
}
