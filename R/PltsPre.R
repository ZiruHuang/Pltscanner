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
PltsPre<- function(suerat_object,cutoff = 0.6911){
  # Load the XGBoost models and Define the  features
  pltmodels<- Pltscanner::pltmodels
  feature_list <- c('CBR3','NQO1','TINAGL1','BAG3','UCHL1','PGF','SACS','NFE2','PADI4','BANK1','S100A8','PPBP','CXCL5','CA2')
  # Get the features from the Seurat object
  features <- pltmodels[[1]]$feature_names
  # Rank the features based on nFeature_RNA and add it to the metadata of the Seurat object
  feature_rank <- rank(suerat_object$nFeature_RNA)
  feature_rank <- feature_rank/length(suerat_object$nFeature_RNA)
  suerat_object@meta.data$feature_rank <- feature_rank
  # Get the number of genes expressed in each cell
  gene_numb<- t(as.matrix(suerat_object$feature_rank))
  # Extract the expression data for the overlapping genes
  overlap_gene<- intersect(feature_list,rownames(suerat_object))
  exp_mat <- as.matrix(Seurat::GetAssayData(suerat_object[overlap_gene,],slot = 'data'))
  # If there are missing genes, add them and fill with zeros
  if(length(overlap_gene)!=length(feature_list)){
    print("Warning message:Not enough features")
    print(length(overlap_gene))
    lack_of_gene <- setdiff(feature_list,overlap_gene)
    lack_of_gene_mat <- matrix(0,nrow=length(lack_of_gene),ncol=ncol(exp_mat))
    rownames(lack_of_gene_mat) <- as.matrix(lack_of_gene)
    exp_mat  <- rbind(exp_mat,lack_of_gene_mat)
  }
  # Add the gene count to the expression data and transpose the matrix
  exp_mat <- rbind(exp_mat,gene_numb)
  rownames(exp_mat)[nrow(exp_mat)] <- "Gene_number"
  exp_mat <- as.data.frame(t(exp_mat))
  # Reorder the expression data to match the features
  exp_mat <- exp_mat[,features]
  # Make predictions for each model
  pre_score <- matrix(NA,nrow =nrow(exp_mat),ncol =10 )
  rownames(pre_score) <- rownames(exp_mat)
  colnames(pre_score) <- paste0("Model",c(1:10))
  for ( i in (1:length(pltmodels))){
    exp_mat <- as.matrix(exp_mat)
    pre_score[,i] <- predict(pltmodels[[i]],newdata = exp_mat)
  }
  # Get the average prediction score across all models
  pre_label <- as.matrix(rowMeans(pre_score))
  # Add the predicted labels to the metadata of the Seurat object
  suerat_object <- Seurat::AddMetaData(suerat_object,pre_label,col.name = 'PltsScore')
  suerat_object$PltsLabel <- ifelse(suerat_object$PltsScore>= cutoff,"Platelet","Other_cell")
  # Return the updated Seurat object
  return(suerat_object)
}
