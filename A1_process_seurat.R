####################################################
##
## Project: amd_aging
##
## D.Dhirachaikulpanich@liverpool.ac.uk
## cyril.lagger@liverpool.ac.uk
##
####################################################
##

## Create Seurat object #####

# We read the Seurat object retrieved from the galaxy portal and
# extract read counts and metadata. As it is unclear how this object
# has been originally processed, we process it again following the
# standard Seurat workflow. Importantly, gene names need to be
# standardized.

# read object retrieved from human cell atlas galaxy portal

amd_seurat_unprocessed <- readRDS(
  paste0(
    path_data,
    "A1_amd_seurat_unprocessed.rds"
  )
)

# understand and extract meta.data

table(amd_seurat_unprocessed$AB) # cell type
table(amd_seurat_unprocessed$AD) # cell type
table(amd_seurat_unprocessed$H) # age
table(amd_seurat_unprocessed$F) # sex
table(amd_seurat_unprocessed$B) # species
table(amd_seurat_unprocessed$J) # left/right eye
table(amd_seurat_unprocessed$L) # death
table(amd_seurat_unprocessed$N) # condition
table(amd_seurat_unprocessed$V) # condition
table(amd_seurat_unprocessed$X) # region
table(amd_seurat_unprocessed$R) # region
table(amd_seurat_unprocessed$orig.ident) # nothing
table(amd_seurat_unprocessed$E) # nothing
table(amd_seurat_unprocessed$I) # nothing
table(amd_seurat_unprocessed$Q) # nothing
table(amd_seurat_unprocessed$U) # nothing
table(amd_seurat_unprocessed$AA) # nothing
table(amd_seurat_unprocessed$C) # url?
table(amd_seurat_unprocessed$G) # url?
table(amd_seurat_unprocessed$K) # url?
table(amd_seurat_unprocessed$M) # url?
table(amd_seurat_unprocessed$O) # url?
table(amd_seurat_unprocessed$S) # url?
table(amd_seurat_unprocessed$W) # url?
table(amd_seurat_unprocessed$Y) # url?
table(amd_seurat_unprocessed$AC) # url?
table(amd_seurat_unprocessed$AE) # url?
table(amd_seurat_unprocessed$D) # ??
table(amd_seurat_unprocessed$P) # ??
table(amd_seurat_unprocessed$T) # ??
table(amd_seurat_unprocessed$Z) # ??

amd_seurat_unprocessed$cell_type <- Idents(amd_seurat_unprocessed)
amd_seurat_unprocessed$age <- amd_seurat_unprocessed$H
amd_seurat_unprocessed$sex <- amd_seurat_unprocessed$F
amd_seurat_unprocessed$condition <- amd_seurat_unprocessed$N
amd_seurat_unprocessed$location <- amd_seurat_unprocessed$X
amd_seurat_unprocessed$is_rpe <- ifelse(
  amd_seurat_unprocessed$cell_type == "RPE", TRUE, FALSE
)

amd_md_cols <- c(
  "cell_type", "age", "sex", "condition", "location", "is_rpe"
)

# clean gene names

mart_human <- biomaRt::useMart(
  "ensembl",
  #host = "https://nov2020.archive.ensembl.org",
  dataset = "hsapiens_gene_ensembl",
  verbose = TRUE
)
mart_amd_hgnc <- biomaRt::getBM(
  attributes = c(
    "hgnc_symbol",
    "ensembl_gene_id"
  ),
  filters = "hgnc_symbol",
  mart = mart_human,
  values = rownames(amd_seurat_unprocessed)
)
setDT(mart_amd_hgnc)
mart_amd_hgnc <- dcast.data.table(
  mart_amd_hgnc,
  hgnc_symbol ~ .,
  fun.aggregate = function(i) {
    paste(i, collapse = "_")
  },
  value.var = "ensembl_gene_id"
)
setnames(
  mart_amd_hgnc,
  old = ".",
  new = "ensembl_gene_ids"
)
mart_amd_ensembl <- biomaRt::getBM(
  attributes = c(
    "hgnc_symbol",
    "ensembl_gene_id"
  ),
  filters = "ensembl_gene_id",
  mart = mart_human,
  values = rownames(amd_seurat_unprocessed)
)
setDT(mart_amd_ensembl)
mart_amd_ensembl <- mart_amd_ensembl[hgnc_symbol != ""]

amd_genes <- data.table(
  orig_symbol = rownames(amd_seurat_unprocessed)
)
any(duplicated(mart_amd_hgnc$hgnc_symbol))
table(mart_amd_hgnc$hgnc_symbol %in% amd_genes$orig_symbol)
amd_genes[
  mart_amd_hgnc,
  on = "orig_symbol==hgnc_symbol",
  ensembl_gene_ids_1 := i.ensembl_gene_ids
]

any(duplicated(mart_amd_ensembl$hgnc_symbol))
any(duplicated(mart_amd_ensembl$ensembl_gene_id))
amd_genes[
  mart_amd_ensembl,
  on = "orig_symbol==ensembl_gene_id",
  hgnc_symbol_1 := i.hgnc_symbol
]

amd_genes[
  ,
  final_symbol := ifelse(
    !is.na(hgnc_symbol_1) & !hgnc_symbol_1 %in% orig_symbol,
    hgnc_symbol_1,
    orig_symbol
  )
]

any(duplicated(amd_genes$final_symbol))
table(
  amd_genes$orig_symbol %in% mart_amd_hgnc$hgnc_symbol
)

# recreate object

amd_counts <- GetAssayData(amd_seurat_unprocessed, slot = "counts")
identical(rownames(amd_counts), amd_genes$orig_symbol)
rownames(amd_counts) <- amd_genes$final_symbol

amd_seurat <- CreateSeuratObject(
  counts = amd_counts,
  project = "AMD_Seurat",
  meta.data =  amd_seurat_unprocessed[[amd_md_cols]]
)

## Seurat workflow ####

# QC and filtering

# probably already done in the original object
# but we rewrite the code for the shake of clarity

VlnPlot(
  amd_seurat,
  features = c("nFeature_RNA", "nCount_RNA"),
  ncol = 2
)
FeatureScatter(
  amd_seurat,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

amd_seurat <- subset(
  amd_seurat,
  subset = nFeature_RNA > 300 & nFeature_RNA < 7000
)

# Normalization

amd_seurat <- NormalizeData(
  amd_seurat,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# HVG

amd_seurat <- FindVariableFeatures(
  amd_seurat,
  selection.method = "vst",
  nfeatures = 2000
)

hvg_top10 <- head(VariableFeatures(amd_seurat), 10)

VariableFeaturePlot(amd_seurat)

# Scaling and PCA

amd_seurat <- ScaleData(amd_seurat)
amd_seurat <- RunPCA(amd_seurat)

DimPlot(amd_seurat, reduction = "pca")

ElbowPlot(amd_seurat)

# Clustering

amd_seurat <- FindNeighbors(amd_seurat, dims = 1:10)
amd_seurat <- FindClusters(amd_seurat, resolution = 0.05)
amd_seurat <- RunUMAP(amd_seurat, dims = 1:10)

DimPlot(amd_seurat, reduction = "umap")

# As expected clustering correspond to initial cell type annotation
table(amd_seurat$cell_type, amd_seurat$seurat_clusters)

Idents(amd_seurat) <- amd_seurat$cell_type

DimPlot(
  amd_seurat,
  reduction = "umap",
  label = TRUE
) + NoLegend()

# Find cell type markers

amd_markers <- FindAllMarkers(
  amd_seurat,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
amd_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)

amd_features <- c(
  "APOD", "PMEL", "CXCR4", "PLP1", "RPE65", "IGKC", "VWF", "TPSB2"
)

## Save processed object ####

saveRDS(
  amd_seurat,
  paste0(
    path_results,
    "A1_amd_seurat.rds"
  )
)

## Extract informative meta.data in data.table ####

amd_md <- setDT(copy(amd_seurat[[]]))
amd_md_summary <- amd_md[, .N, by = c("age", "sex", "condition", "location")]
amd_md_summary_ct <- amd_md[
  ,
  .N,
  by = c("cell_type", "age", "sex", "condition", "location")
]

## Data for Figure 1A ####

amd_md[, .N, by = c("age", "sex", "condition")]

## Data or Figure 1B ####

amd_fig1B <- amd_md[, .N, by = c("cell_type")][order(-N)]
amd_fig1B[
  ,
  pct := N / sum(N) * 100
]
amd_fig1B

## Figure 1C-D ####

amd_seurat$condition_abbr <- ifelse(
  amd_seurat$condition == "normal",
  "Normal",
  "AMD"
)

amd_fig1CD <- cowplot::plot_grid(
  plotlist = list(
    DimPlot(
      amd_seurat,
      reduction = "umap",
      group.by = "condition_abbr",
      pt.size = 0.5
    ) + ggtitle(""),
        DimPlot(
      amd_seurat,
      label = TRUE
    ) + ggtitle("") + NoLegend()
  ),
  ncol = 1
)

ggsave(
  paste0(
    path_results,
    "images/Fig1CD.png"
  ),
  amd_fig1CD,
  height = 9,
  width = 5
)

## Figure 1E ####

amd_fig1E_list <- VlnPlot(
  amd_seurat,
  features = amd_features,
  ncol = 2,
  pt.size = 0,
  combine = FALSE
)
amd_fig1E_list <- lapply(
  amd_fig1E_list,
  function(i) {
    i + NoLegend() + theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 10),
      plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")
    ) + xlab(
      ""
    ) + ylab(
      "Expression level"
    )
  }
)

amd_fig1E <- cowplot::plot_grid(
  plotlist = amd_fig1E_list,
  ncol = 2
)

ggsave(
  paste0(
    path_results,
    "images/Fig1E.png"
  ),
  amd_fig1E,
  height = 9,
  width = 5
)
