# Ptscanner

**Updated: April 2023**

Pltscanner is a machine learning-based R package designed for automated annotation of platelets/megakaryocytes in scRNA-seq data. 

# Install
```R
devtools::install_github('ZiruHuang/Ptscanner')
```

# Dependency
```R
library(Seurat)
library(xgboost)
```

# Usage

```R
library(Ptscanner)
# Pltscanner accepts Seurat objects as input data and returns a Seurat object with platelet cell labels and prediction scores for each cell.
data("pbmc_small")
pbmc_small <- PltsPre(pbmc_small)
# The prediction results will be stored in the Metadata.
DimPlot(pbmc_small,group.by = "PltsLabel")
FeaturePlot(pbmc_small,features = "PltsScore")
```
# Contributors

Ptscanner was developed by Ziru Huang. Please contact Ziru Huang: ziru_huang@163.com for any questions or suggestions.