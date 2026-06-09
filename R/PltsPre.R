#' PltsPre is a function to predict platelet cells  in single-cell RNA sequencing (scRNA-seq) data.
#'
#' @param suerat_object  A Seurat object that metadata contains nFeature_RNA.
#' @param cutoff Numeric. Cutoff value for binary prediction. Default is 0.6911.
#'
#' @return  A Seurat object with predicted cell labels and scores.
#'
#' @examples
#' # Example usage of the function
#' data("pbmc_small")
#' pbmc_small <- PltsPre(pbmc_small)
#'
#' @import xgboost
#' @import Seurat
#' @importFrom stats predict
#' @importFrom utils data
#'
#' @export
#'

PltsPre <- function(suerat_object, cutoff = 0.6911){

  # Load models 
  pltmodels <- load_models_from_inst()


  # Feature list 
  feature_list <- c(
    'CBR3','NQO1','TINAGL1','BAG3','UCHL1','PGF',
    'SACS','NFE2','PADI4','BANK1','S100A8','PPBP',
    'CXCL5','CA2','Gene_number'
  )

  # Cell-level feature (nFeature_RNA rank)
  feature_rank <- rank(suerat_object$nFeature_RNA)
  feature_rank <- feature_rank / length(suerat_object$nFeature_RNA)
  suerat_object@meta.data$feature_rank <- feature_rank
  gene_numb <- t(as.matrix(feature_rank))

  # Expression matrix (Seurat v4/v5 safe)
  overlap_gene <- intersect(feature_list, rownames(suerat_object))
  exp_mat <- get_seurat_matrix(suerat_object, overlap_gene)

  # Handle missing genes
  if (length(overlap_gene) != length(feature_list)) {
    lack_of_gene <- setdiff(feature_list, overlap_gene)
    if (length(lack_of_gene) > 0) {
      lack_mat <- matrix(
        0,
        nrow = length(lack_of_gene),
        ncol = ncol(exp_mat)
      )
      rownames(lack_mat) <- lack_of_gene
      exp_mat <- rbind(exp_mat, lack_mat)
    }
  }

  # Add Gene_number
  exp_mat <- rbind(exp_mat, gene_numb)
  rownames(exp_mat)[nrow(exp_mat)] <- "Gene_number"

  # Feature ordering
  exp_mat <- exp_mat[feature_list, , drop = FALSE]
  exp_mat <- t(exp_mat)
  exp_mat <- as.matrix(exp_mat)

  # XGBoost prediction 
  dtest <- xgboost::xgb.DMatrix(exp_mat)
  pre_score <- matrix(
    NA,
    nrow = nrow(exp_mat),
    ncol = length(pltmodels)
  )
  colnames(pre_score) <- paste0("Model", seq_along(pltmodels))
  rownames(pre_score) <- rownames(exp_mat)
  for (i in seq_along(pltmodels)) {
    pre_score[, i] <- predict(pltmodels[[i]], dtest)
  }

  #  Aggregate score
  plts_score <- rowMeans(pre_score)
  #  Add to Seurat object
  suerat_object <- Seurat::AddMetaData(
    suerat_object,
    metadata = plts_score,
    col.name = "PltsScore"
  )
  suerat_object$PltsLabel <- ifelse(
    suerat_object$PltsScore >= cutoff,
    "Platelet",
    "Other_cell"
  )

  return(suerat_object)
}
