####################################################
##
## Project: amd_aging
##
## D.Dhirachaikulpanich@liverpool.ac.uk
## cyril.lagger@liverpool.ac.uk
##
####################################################
##

## Options ####

#on the server
plan(multicore, workers = 20)

## run scDiffCom analyses ####

icc_scdiffcom <- list(
  icc_detec_all = run_interaction_analysis(
    seurat_object = seurat_object,
    LRI_species = "human",
    seurat_celltype_id = "cell_type",
    seurat_condition_id = NULL,
    scdiffcom_object_name = "icc_detec_all",
    iterations = 10000
  ),
  icc_detec_54 = run_interaction_analysis(
    seurat_object = subset(seurat_object, subset = age == "54 year"),
    LRI_species = "human",
    seurat_celltype_id = "cell_type",
    seurat_condition_id = NULL,
    scdiffcom_object_name = "icc_detec_54",
    iterations = 10000
  ),
  icc_detec_82 = run_interaction_analysis(
    seurat_object = subset(seurat_object, subset = age == "82 year"),
    LRI_species = "human",
    seurat_celltype_id = "cell_type",
    seurat_condition_id = NULL,
    scdiffcom_object_name = "icc_detec_82",
    iterations = 10000
  ),
  icc_detec_79 = run_interaction_analysis(
    seurat_object = subset(seurat_object, subset = age == "79 year"),
    LRI_species = "human",
    seurat_celltype_id = "cell_type",
    seurat_condition_id = NULL,
    scdiffcom_object_name = "icc_detec_79",
    iterations = 10000
  ),
  icc_detec_healthy = run_interaction_analysis(
    seurat_object = subset(
      seurat_object,
      subset = age %in% c("54 year", "82 year")
    ),
    LRI_species = "human",
    seurat_celltype_id = "cell_type",
    seurat_condition_id = NULL,
    scdiffcom_object_name = "icc_detec_healthy",
    iterations = 10000
  ),
  icc_diff_54_vs_79 = run_interaction_analysis(
    seurat_object = subset(
      seurat_object,
      subset = age %in% c("54 year", "79 year")
    ),
    LRI_species = "human",
    seurat_celltype_id = "cell_type",
    seurat_condition_id = list(
      column_name = "age",
      cond1_name = "54 year",
      cond2_name = "79 year"
    ),
    scdiffcom_object_name = "icc_diff_54_vs_79",
    iterations = 10000
  ),
  icc_diff_54_vs_82 = run_interaction_analysis(
    seurat_object = subset(
      seurat_object,
      subset = age %in% c("54 year", "82 year")
    ),
    LRI_species = "human",
    seurat_celltype_id = "cell_type",
    seurat_condition_id = list(
      column_name = "age",
      cond1_name = "54 year",
      cond2_name = "82 year"
    ),
    scdiffcom_object_name = "icc_diff_54_vs_82",
    iterations = 10000
  ),
  icc_diff_82_vs_79 = run_interaction_analysis(
    seurat_object = subset(
      seurat_object,
      subset = age %in% c("82 year", "79 year")
    ),
    LRI_species = "human",
    seurat_celltype_id = "cell_type",
    seurat_condition_id = list(
      column_name = "age",
      cond1_name = "82 year",
      cond2_name = "79 year"
    ),
    scdiffcom_object_name = "icc_diff_82_vs_79",
    iterations = 10000
  ),
  icc_diff_healthy_vs_amd = run_interaction_analysis(
    seurat_object = seurat_object,
    LRI_species = "human",
    seurat_celltype_id = "cell_type",
    seurat_condition_id = list(
      column_name = "condition",
      cond1_name = "normal",
      cond2_name = "wet macular degeneration"
    ),
    scdiffcom_object_name = "icc_diff_healthy_vs_amd",
    iterations = 10000
  ),
  icc_diff_macula_vs_peri = run_interaction_analysis(
    seurat_object = seurat_object,
    LRI_species = "human",
    seurat_celltype_id = "cell_type",
    seurat_condition_id = list(
      column_name = "location",
      cond1_name = "macula lutea",
      cond2_name = "peripheral region of retina"
    ),
    scdiffcom_object_name = "icc_diff_macula_vs_peri",
    iterations = 10000
  )
)

## Save/Read scDiffcom results ####

saveRDS(
  icc_scdiffcom,
  paste0(
    path_results,
    "C1_icc_scdiffcom.rds"
  )
)

icc_scdiffcom <- readRDS(
  paste0(
    path_results,
    "B1_icc_scdiffcom.rds"
  )
)


## upset plot of detected CCIs ####

icc_scd_ccis <- data.table(
  cci = unique(
    c(
      icc_scdiffcom$icc_detec_all@cci_table_detected$CCI,
      icc_scdiffcom$icc_detec_54@cci_table_detected$CCI,
      icc_scdiffcom$icc_detec_82@cci_table_detected$CCI,
      icc_scdiffcom$icc_detec_79@cci_table_detected$CCI
    )
  )
)
icc_scd_ccis[
  ,
  detec_all := ifelse(
    cci %in% icc_scdiffcom$icc_detec_all@cci_table_detected$CCI,
    T,
    F
  )
]
icc_scd_ccis[
  ,
  detec_54 := ifelse(
    cci %in% icc_scdiffcom$icc_detec_54@cci_table_detected$CCI,
    T,
    F
  )
]
icc_scd_ccis[
  ,
  detec_82 := ifelse(
    cci %in% icc_scdiffcom$icc_detec_82@cci_table_detected$CCI,
    T,
    F
  )
]
icc_scd_ccis[
  ,
  detec_79 := ifelse(
    cci %in% icc_scdiffcom$icc_detec_79@cci_table_detected$CCI,
    T,
    F
  )
]

ComplexUpset::upset(
  icc_scd_ccis,
  c("detec_all", "detec_54", "detec_82", "detec_79")
)
ComplexUpset::upset(
  icc_scd_ccis,
  c("detec_54", "detec_82", "detec_79")
)

## Results exploration ####

table(icc_diff_macula_vs_peri@cci_table_detected$REGULATION)
table(icc_diff_healthy_vs_amd@cci_table_detected$REGULATION)
table(icc_diff_54_vs_79@cci_table_detected$REGULATION)
table(icc_diff_54_vs_82@cci_table_detected$REGULATION)
table(icc_diff_82_vs_79@cci_table_detected$REGULATION)

PlotORA(
  icc_scdiffcom$icc_diff_54_vs_79,
  category = "LRI"
)
PlotORA(
  icc_diff_82_vs_79,
  category = "LRI"
)
PlotORA(
  icc_diff_54_vs_79,
  category = "LRI",
  regulation = "DOWN"
)
PlotORA(
  icc_diff_82_vs_79,
  category = "LRI",
  regulation = "DOWN"
)

cci_detec_all <- GetTableCCI(
  icc_detec_all,
  type = "detected",
  simplified = TRUE
)

cci_detec_common <- cci_detec_all[
  CCI %in% icc_all_ccis[detec_all & detec_54 & detec_79 & detec_82]$cci # &
  # grepl("RPE", CCI)
]

all_cci_1 <- GetTableCCI(
  icc_scdiffcom$icc_detec_all, "detected", TRUE
)
all_cci_2 <- all_cci_1[
  ER_CELLTYPES %in% "RPE_endothelial"
]



## cell types based on senesence scores ####

ftable(
  seurat_object$age,
  seurat_object$score_seurat_mayo_top20pct,
  seurat_object$cell_type
)

seurat_object$cell_type_sen <- ifelse(
  # seurat_object$score_seurat_mayo_top50pct == TRUE, # &
  seurat_object$score_seurat_kasit_20pct == TRUE,
  paste0(seurat_object$cell_type, "_sen"),
  paste0(seurat_object$cell_type, "")
)
table(seurat_object$cell_type_sen)

FeaturePlot(
  seurat_object,
  features = "CDKN1A"
)
DimPlot(seurat_object)

## scDiffCom sen Analysis ####

icc_detec_all_sen <- run_interaction_analysis(
  seurat_object = seurat_object,
  LRI_species = "human",
  seurat_celltype_id = "cell_type_sen",
  seurat_condition_id = NULL,
  scdiffcom_object_name = "icc_detec_all_sen",
  iterations = 10000
)

saveRDS(
  icc_detec_all_sen,
  "../results/icc_scdiffcom_sen_kasit_20pct.rds"
)

icc_detec_all_sen

cci_table <- GetTableCCI(icc_detec_all_sen, "detected", FALSE)

test <- cci_table[ER_CELLTYPES == "RPE-sen_endothelial-sen"]
test2 <- cci_table[ER_CELLTYPES == "endothelial-sen_RPE-sen"]
testb <- cci_table[ER_CELLTYPES == "RPE_endothelial"]

test <- cci_table[
  ER_CELLTYPES %in% c(
    "RPE-sen_endothelial-sen",
    "RPE_endothelial"
  )
]
test2 <- dcast.data.table(
  test[, c("ER_CELLTYPES", "LRI", "CCI_SCORE")],
  formula = LRI ~ ER_CELLTYPES,
  value.var = "CCI_SCORE"
)

test2



setdiff(
  test$LRI,
  testb$LRI
)

setdiff(
  testb$LRI,
  test$LRI
)

intersect(
  test$LRI,
  testb$LRI
)

lri_rpe_sen <- cci_table[ER_CELLTYPES == "RPE-sen_RPE-sen"]$LRI
lri_rpe <- cci_table[ER_CELLTYPES == "RPE_RPE"]$LRI

sort(setdiff(
  lri_rpe_sen,
  lri_rpe
))
sort(setdiff(
  lri_rpe,
  lri_rpe_sen
))
sort(intersect(
  lri_rpe,
  lri_rpe_sen
))

library(clusterProfiler)
library(org.Hs.eg.db)

cci_table[, .N, by = ER_CELLTYPES][order(-N)]

t1 <- unique(
  unlist(
    cci_table[ER_CELLTYPES == "RPE-sen_RPE-sen"][, c("LIGAND_1", "RECEPTOR_1")]
  )
)
t1 <- t1[!is.na(t1)]

lri_rpe_sen_go <- enrichGO(
  gene = t1,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "ALL",
  readable = FALSE
)
lri_rpe_sen_go_dt <- setDT(
  lri_rpe_sen_go@result
)


t2 <- unique(
  unlist(
    cci_table[ER_CELLTYPES == "RPE_RPE"][, c("LIGAND_1", "RECEPTOR_1")]
  )
)
t2 <- t2[!is.na(t2)]

lri_rpe_go <- enrichGO(
  gene = t2,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "ALL",
  readable = FALSE
)
lri_rpe_go_dt <- setDT(
  lri_rpe_go@result
)

setdiff(
  lri_rpe_sen_go_dt$Description,
  lri_rpe_go_dt$Description
)[1:20]