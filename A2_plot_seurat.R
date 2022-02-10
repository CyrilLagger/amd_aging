####################################################
##
## Project: amd_aging
##
## D.Dhirachaikulpanich@liverpool.ac.uk
## cyril.lagger@liverpool.ac.uk
##
####################################################
##

## DimPlots ####

DimPlot(
  amd_seurat,
  reduction = "umap",
  label = TRUE
)
DimPlot(
  amd_seurat,
  reduction = "umap",
  group.by = "sex"
)
DimPlot(
  amd_seurat,
  reduction = "umap",
  group.by = "age"
)
DimPlot(
  amd_seurat,
  reduction = "umap",
  group.by = "condition"
)
DimPlot(
  amd_seurat,
  reduction = "umap",
  group.by = "location"
)