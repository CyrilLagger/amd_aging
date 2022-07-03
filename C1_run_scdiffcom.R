####################################################
##
## Project: amd_aging
##
## D.Dhirachaikulpanich@liverpool.ac.uk
## cyril.lagger@liverpool.ac.uk
##
####################################################
##

## run scDiffCom analyses on normal cell types ####

scd_healthy <- run_interaction_analysis(
  seurat_object = subset(
    amd_seurat,
    subset = age %in% c("54 year", "82 year")
  ),
  LRI_species = "human",
  seurat_celltype_id = "cell_type",
  seurat_condition_id = NULL,
  scdiffcom_object_name = "icc_detec_healthy",
  threshold_quantile_score = 0.0,
  iterations = 10000
)

scd_amd <- run_interaction_analysis(
  seurat_object = subset(amd_seurat, subset = age == "79 year"),
  LRI_species = "human",
  seurat_celltype_id = "cell_type",
  seurat_condition_id = NULL,
  scdiffcom_object_name = "icc_detec_amd",
  threshold_quantile_score = 0.0,
  iterations = 10000
)

saveRDS(
  scd_healthy,
  paste0(
    path_results,
    "C1_scd_healthy.rds"
  )
)

saveRDS(
  scd_amd,
  paste0(
    path_results,
    "C1_scd_amd.rds"
  )
)

## run scDiffCom analyses on normal and senescent cell types ####

scd_sen <- run_interaction_analysis(
  seurat_object = amd_seurat,
  LRI_species = "human",
  seurat_celltype_id = "cell_type_senescence_gsea",
  seurat_condition_id = NULL,
  scdiffcom_object_name = "icc_detec_all_sen_gsea",
  iterations = 10000,
  threshold_min_cells = 11,
  threshold_quantile_score = 0.0
)

saveRDS(
  scd_sen,
  paste0(
    path_results,
    "C1_scd_sen.rds"
  )
)
