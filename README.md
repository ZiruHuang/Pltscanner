# Pltscanner

Pltscanner is an R package for identifying platelet cells in single-cell RNA-seq (scRNA-seq) data using a machine learning approach.

The method is based on an ensemble of XGBoost models and performs cell-level prediction directly on gene expression data, without requiring clustering or dimensionality reduction.

---

## Installation

```r
install.packages("devtools")
devtools::install_github("ZiruHuang/Pltscanner")
```

---

## Dependencies

```r
install.packages("Seurat")
install.packages("xgboost")
```

The package has been tested with XGBoost versions from 1.4.x to 3.2.x.

---

## Usage

```r
library(Pltscanner)
library(Seurat)

data("pbmc_small")

pbmc_small <- PltsPre(pbmc_small)

DimPlot(pbmc_small, group.by = "PltsLabel")
FeaturePlot(pbmc_small, features = "PltsScore")
```

To use a custom cutoff:

```r
pbmc_small <- PltsPre(pbmc_small, cutoff = 0.6911)
```

---

## Output

The function adds two columns to the Seurat object metadata:

* `PltsScore`: prediction score (0–1)
* `PltsLabel`: predicted label ("Platelet" or "Other_cell")

---

## Notes

* Input must be a Seurat object with an RNA assay.
* `nFeature_RNA` is required and used as an additional feature.
* Missing genes are filled with zeros during prediction.

---

## Changelog

### v0.1.1 (June 2026)

* Improved compatibility with XGBoost (1.4.x–3.2.x)
* Switched to JSON model format
* Updated for Seurat v5

### v0.1.0 (April 2023)

* Initial release

---

## Contributors

Pltscanner was developed by Ziru Huang. Please contact Ziru Huang: ziru_huang@163.com for any questions or suggestions.
