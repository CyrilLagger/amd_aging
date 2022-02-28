####################################################
##
## Project: amd_aging
##
## D.Dhirachaikulpanich@liverpool.ac.uk
## cyril.lagger@liverpool.ac.uk
##
####################################################
##

## Library ####

library(Seurat)
library(scDiffCom)
library(CellChat)
library(data.table)
library(ggplot2)
library(future)
library(ComplexUpset)
library(escape)
library(ggVennDiagram)
library(ggExtra)
library(limma)

## Options ####

plan(multisession, workers = 20)

## Directory paths ####

path_data <- "../data/"
path_results <- "../results/"

## Read Seurat object #####

amd_seurat <- readRDS(
  paste0(
    path_data,
    "/A1_amd_seurat.rds"
  )
)
amd_seurat

## Rescale data to make the object lighter ####

amd_seurat <- ScaleData(amd_seurat)

GetAssayData(amd_seurat)[1:5, 1:5]
GetAssayData(amd_seurat, slot = "counts")[1:5, 1:5]
head(rownames(amd_seurat))

## Understand the meta.data ####

# TODO ask why this strange convention is used

table(amd_seurat$AB) # cell type
table(amd_seurat$AD) # cell type
table(amd_seurat$H) # age
table(amd_seurat$F) # sex
table(amd_seurat$B) # species
table(amd_seurat$J) # left/right eye
table(amd_seurat$L) # death
table(amd_seurat$N) # condition
table(amd_seurat$V) # condition
table(amd_seurat$X) # region
table(amd_seurat$R) # region

table(amd_seurat$orig.ident) # nothing
table(amd_seurat$E) # nothing
table(amd_seurat$I) # nothing
table(amd_seurat$Q) # nothing
table(amd_seurat$U) # nothing
table(amd_seurat$AA) # nothing
table(amd_seurat$C) # url?
table(amd_seurat$G) # url?
table(amd_seurat$K) # url?
table(amd_seurat$M) # url?
table(amd_seurat$O) # url?
table(amd_seurat$S) # url?
table(amd_seurat$W) # url?
table(amd_seurat$Y) # url?
table(amd_seurat$AC) # url?
table(amd_seurat$AE) # url?
table(amd_seurat$D) # ??
table(amd_seurat$P) # ??
table(amd_seurat$T) # ??
table(amd_seurat$Z) # ??

## Rename informative meta.data ####

amd_seurat$cell_type <- Idents(amd_seurat)
amd_seurat$age <- amd_seurat$H
amd_seurat$sex <- amd_seurat$F
amd_seurat$condition <- amd_seurat$N
amd_seurat$location <- amd_seurat$X
amd_seurat$is_rpe <- ifelse(
  amd_seurat$cell_type == "RPE", TRUE, FALSE
)

## Extract informative meta.data in data.table ####

amd_md <- setDT(copy(amd_seurat[[]]))
amd_md <- amd_md[, c("cell_type", "age", "sex", "condition", "location")]

amd_md_summary <- amd_md[, .N, by = c("age", "sex", "condition", "location")]
amd_md_summary_ct <- amd_md[
  ,
  .N,
  by = c("cell_type", "age", "sex", "condition", "location")
]

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
