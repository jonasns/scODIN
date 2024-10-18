#' Apply a score filter to kNN predictions
#'
#' This function applies a score filter to the predictions in a specified assay of a Seurat object.
#' It assigns "Unassigned" to any cell where the maximum prediction score is below the specified 
#' threshold and returns the name of the cell type with the highest score for those above the threshold.
#'
#' @param scData A Seurat object containing the predictions.
#' @param assay A character string specifying the assay from which to retrieve the predictions. 
#'        Defaults to "predictions".
#' @param slot A character string specifying the slot from which to retrieve the data within the assay.
#'        Defaults to "data".
#' @param score.filter A numeric value indicating the threshold below which cells are labeled as "Unassigned". 
#'        Defaults to 0.5. This should be adjusted empirically, as few cell types in the core requires a higher cutoff
#'        than a core with many cell types.
#' 
#' @return A character vector of predicted cell types, with "Unassigned" for cells that do not meet 
#'         the score threshold.
#' 
#' @examples
#' predictions <- apply_score_filter(scData = scData, 
#'                                    assay = "predictions", 
#'                                    slot = "data", 
#'                                    score.filter = 0.5)
#' 
#' @export
apply_score_filter <- function(scData, assay = "predictions", slot = "data", score.filter = 0.5) {
  dat <- GetAssayData(scData[[assay]], layer = slot)
  predictions <- apply(
    X = dat,
    MARGIN = 2,
    FUN = function(x) {
      if (max(x) < score.filter) {
        return("Unassigned")
      } else {
        x <- x[which.max(x)]
        return(names(x)[which.max(x)])
      }
    }
  )
  return(predictions)
}