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

# on the server
plan(multicore, workers = 20)

## run normal scDiffCom analyses ####

icc_scdiffcom <- list(
  icc_detec_all = run_interaction_analysis(
    seurat_object = amd_seurat,
    LRI_species = "human",
    seurat_celltype_id = "cell_type",
    seurat_condition_id = NULL,
    scdiffcom_object_name = "icc_detec_all",
    iterations = 10000
  ),
  icc_detec_54 = run_interaction_analysis(
    seurat_object = subset(amd_seurat, subset = age == "54 year"),
    LRI_species = "human",
    seurat_celltype_id = "cell_type",
    seurat_condition_id = NULL,
    scdiffcom_object_name = "icc_detec_54",
    iterations = 10000
  ),
  icc_detec_82 = run_interaction_analysis(
    seurat_object = subset(amd_seurat, subset = age == "82 year"),
    LRI_species = "human",
    seurat_celltype_id = "cell_type",
    seurat_condition_id = NULL,
    scdiffcom_object_name = "icc_detec_82",
    iterations = 10000
  ),
  icc_detec_79 = run_interaction_analysis(
    seurat_object = subset(amd_seurat, subset = age == "79 year"),
    LRI_species = "human",
    seurat_celltype_id = "cell_type",
    seurat_condition_id = NULL,
    scdiffcom_object_name = "icc_detec_79",
    iterations = 10000
  ),
  icc_detec_healthy = run_interaction_analysis(
    seurat_object = subset(
      amd_seurat,
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
      amd_seurat,
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
      amd_seurat,
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
      amd_seurat,
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
    seurat_object = amd_seurat,
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
    seurat_object = amd_seurat,
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

## Save/Read normal scDiffcom results ####

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
    "C1_icc_scdiffcom.rds"
  )
)

## scDiffCom senescence analysis ####

icc_sen_scdiffcom <- list(
  icc_detec_all = run_interaction_analysis(
    seurat_object = amd_seurat,
    LRI_species = "human",
    seurat_celltype_id = "cell_type_senescence",
    seurat_condition_id = NULL,
    scdiffcom_object_name = "icc_detec_all_sen",
    iterations = 1000
  ),
  icc_detec_54 = run_interaction_analysis(
    seurat_object = subset(amd_seurat, subset = age == "54 year"),
    LRI_species = "human",
    seurat_celltype_id = "cell_type_senescence",
    seurat_condition_id = NULL,
    scdiffcom_object_name = "icc_detec_54_sen",
    iterations = 1000
  ),
  icc_detec_82 = run_interaction_analysis(
    seurat_object = subset(amd_seurat, subset = age == "82 year"),
    LRI_species = "human",
    seurat_celltype_id = "cell_type_senescence",
    seurat_condition_id = NULL,
    scdiffcom_object_name = "icc_detec_82_sen",
    iterations = 1000
  ),
  icc_detec_79 = run_interaction_analysis(
    seurat_object = subset(amd_seurat, subset = age == "79 year"),
    LRI_species = "human",
    seurat_celltype_id = "cell_type_senescence",
    seurat_condition_id = NULL,
    scdiffcom_object_name = "icc_detec_79_sen",
    iterations = 1000
  ),
  icc_detec_healthy = run_interaction_analysis(
    seurat_object = subset(
      amd_seurat,
      subset = age %in% c("54 year", "82 year")
    ),
    LRI_species = "human",
    seurat_celltype_id = "cell_type_senescence",
    seurat_condition_id = NULL,
    scdiffcom_object_name = "icc_detec_healthy_sen",
    iterations = 1000
  ),
  icc_diff_54_vs_79 = run_interaction_analysis(
    seurat_object = subset(
      amd_seurat,
      subset = age %in% c("54 year", "79 year")
    ),
    LRI_species = "human",
    seurat_celltype_id = "cell_type_senescence",
    seurat_condition_id = list(
      column_name = "age",
      cond1_name = "54 year",
      cond2_name = "79 year"
    ),
    scdiffcom_object_name = "icc_diff_54_vs_79_sen",
    iterations = 1000
  ),
  icc_diff_54_vs_82 = run_interaction_analysis(
    seurat_object = subset(
      amd_seurat,
      subset = age %in% c("54 year", "82 year")
    ),
    LRI_species = "human",
    seurat_celltype_id = "cell_type_senescence",
    seurat_condition_id = list(
      column_name = "age",
      cond1_name = "54 year",
      cond2_name = "82 year"
    ),
    scdiffcom_object_name = "icc_diff_54_vs_82_sen",
    iterations = 1000
  ),
  icc_diff_82_vs_79 = run_interaction_analysis(
    seurat_object = subset(
      amd_seurat,
      subset = age %in% c("82 year", "79 year")
    ),
    LRI_species = "human",
    seurat_celltype_id = "cell_type_senescence",
    seurat_condition_id = list(
      column_name = "age",
      cond1_name = "82 year",
      cond2_name = "79 year"
    ),
    scdiffcom_object_name = "icc_diff_82_vs_79_sen",
    iterations = 1000
  ),
  icc_diff_healthy_vs_amd = run_interaction_analysis(
    seurat_object = amd_seurat,
    LRI_species = "human",
    seurat_celltype_id = "cell_type_senescence",
    seurat_condition_id = list(
      column_name = "condition",
      cond1_name = "normal",
      cond2_name = "wet macular degeneration"
    ),
    scdiffcom_object_name = "icc_diff_healthy_vs_amd_sen",
    iterations = 1000
  ),
  icc_diff_macula_vs_peri = run_interaction_analysis(
    seurat_object = amd_seurat,
    LRI_species = "human",
    seurat_celltype_id = "cell_type_senescence",
    seurat_condition_id = list(
      column_name = "location",
      cond1_name = "macula lutea",
      cond2_name = "peripheral region of retina"
    ),
    scdiffcom_object_name = "icc_diff_macula_vs_peri_sen",
    iterations = 1000
  )
)

## Save/Read senescence scDiffcom results ####

saveRDS(
  icc_sen_scdiffcom,
  paste0(
    path_results,
    "C1_icc_sen_scdiffcom.rds"
  )
)

icc_sen_scdiffcom <- readRDS(
  paste0(
    path_results,
    "C1_icc_sen_scdiffcom.rds"
  )
)