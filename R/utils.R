get_seurat_matrix <- function(obj, genes) {

  Seurat::DefaultAssay(obj) <- "RNA"

  mat <- tryCatch(
    {
      Seurat::GetAssayData(obj[genes, ], assay = "RNA", slot = "data")
    },
    error = function(e) {
      Seurat::GetAssayData(obj[genes, ], assay = "RNA")
    }
  )

  as.matrix(mat)
}
load_models_from_inst <- function() {

  model_files <- system.file(
    "extdata/models",
    paste0("model_", 1:10, ".json"),
    package = "Pltscanner"
  )

  lapply(model_files, xgboost::xgb.load)
}
