####################################################
##
## Project: amd_aging
##
## D.Dhirachaikulpanich@liverpool.ac.uk
## cyril.lagger@liverpool.ac.uk
##
####################################################
##

## GEO query ####

rpe_gset <- getGEO(
  "GSE29801",
  GSEMatrix = TRUE,
  AnnotGPL = TRUE
)
length(rpe_gset)
rpe_gset <- rpe_gset[[1]]
fvarLabels(rpe_gset) <- make.names(fvarLabels(rpe_gset))

## Explore dataset ####

experimentData(rpe_gset)
rownames(exprs(rpe_gset))
colnames(exprs(rpe_gset))
str(exprs(rpe_gset))
sampleNames(rpe_gset)
str(pData(phenoData(rpe_gset)))
varMetadata(phenoData(rpe_gset))
str(pData(featureData(rpe_gset)))
varMetadata(featureData(rpe_gset))

## Keep RPE-choroid and remove Retina samples

rpe_gsms <- paste0(
  "02020202020202022020020202020220202202002020202022",
  "00202022022002002022020020200020002020020202021313",
  "13131313131313131313113131313131313131313131313131",
  "3113131313131313131313131XXXXXXXXXXXXXXXXXXXXXXXXX",
  "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
  "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
)
rpe_sml <- c()
for (i in 1:nchar(rpe_gsms)) {
  rpe_sml[i] <- substr(rpe_gsms, i, i)
}
rpe_sel <- which(rpe_sml != "X")

rpe_gset <- rpe_gset[, rpe_sel]

table(pData(phenoData(rpe_gset))$characteristics_ch1)

## Transform expression values ####

rpe_ex <- exprs(rpe_gset)
rpe_qx <- as.numeric(
  quantile(
    rpe_ex,
    c(0., 0.25, 0.5, 0.75, 0.99, 1.0),
    na.rm = T
  )
)

rpe_LogC <- (
  rpe_qx[5] > 100
) || (
  rpe_qx[6] - rpe_qx[1] > 50 && rpe_qx[2] > 0
) || (
  rpe_qx[2] > 0 && rpe_qx[2] < 1 && rpe_qx[4] > 1 && rpe_qx[4] < 2
)

if (rpe_LogC) {
  rpe_ex[which(rpe_ex <= 0)] <- NaN
  exprs(rpe_gset) <- log2(rpe_ex)
}

## Clean meta.data ####

rpe_gset$age <- as.numeric(rpe_gset$`age (years):ch1`)
rpe_gset$sex <- as.factor(rpe_gset$`gender:ch1`)
rpe_gset$rin <- as.numeric(rpe_gset$`rna integrity number (rin):ch1`)
rpe_gset$mac <- as.factor(rpe_gset$`tissue:ch1`)
rpe_gset$disease <- as.factor(rpe_gset$`ocular disease:ch1`)

## Differential expression analysis on 96 healthy tissue ####

rpe_eel <- which(rpe_sml == "0" | rpe_sml == "2")
rpe_eml <- rpe_sml[rpe_eel]
rpe_eeset <- rpe_gset[, rpe_eel]
table(pData(phenoData(rpe_eeset))$characteristics_ch1)
table(pData(phenoData(rpe_eeset))$disease)

rpe_design <- model.matrix(~ age + sex + mac + 0, rpe_eeset)
rpe_fit <- lmFit(rpe_eeset, rpe_design)
rpe_fit <- eBayes(rpe_fit, 0.01)
rpe_tT <- topTable(
  rpe_fit,
  adjust = "fdr",
  coef = "age",
  sort.by = "B",
  number = 30000
)
setDT(rpe_tT)

rpe_results <- decideTests(rpe_fit)

fwrite(
  rpe_tT,
  paste0(
    path_results,
    "D1_limma_deg.csv"
  )
)

## Load and compare to pre-saved limma results ####

rpe_limma_de <- copy(rpe_tT)

rpe_limma_check <- fread(
  paste0(
    path_data,
    "D1_limma_data.csv"
  )
)
table(rpe_limma_de$ID %in%
rpe_limma_check$ID)
rpe_limma_check[
  rpe_limma_de,
  on = "ID",
  B2 := i.B
]
cor(rpe_limma_check$B, rpe_limma_check$B2)

## Load scDiffCom and cellchat genes ####

scd_genes <- unique(
  unlist(
    LRI_human$LRI_curated[
      ,
      c(
        "LIGAND_1", "LIGAND_2",
        "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3"
      )
    ]
  )
)
scd_genes <- scd_genes[!is.na(scd_genes)]
length(scd_genes)

cc_genes <- unique(
  unlist(
    LRI_human$LRI_curated[
      grepl("CellChat", DATABASE)
      ,
      c(
        "LIGAND_1", "LIGAND_2",
        "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3"
      )
    ]
  )
)
cc_genes <- cc_genes[!is.na(cc_genes)]
length(cc_genes)

## Subset aging genes by ICC, Sup Table 5 ####

table(
  scd_genes %in% unique(rpe_limma_de$Gene.symbol)
)
table(
    cc_genes %in% unique(rpe_limma_de$Gene.symbol)
)

#for scDiffCom
rpe_limma_de_scd <- rpe_limma_de[
  Gene.symbol %in% scd_genes,
  c("Gene.symbol", "logFC", "P.Value", "adj.P.Val")
]
rpe_limma_de_scd[, bh_value := p.adjust(P.Value, method = "BH")]
rpe_limma_de_scd_significant <- rpe_limma_de_scd[
  bh_value <= 0.1 & abs(logFC) >= log2(1.5) / 70
]
rpe_limma_de_scd_significant[
  ,
  pi_score_abs := -log10(bh_value) * abs(logFC)
]
rpe_limma_de_scd_significant <- rpe_limma_de_scd_significant[
  order(-pi_score_abs)
]
rpe_limma_de_scd_significant <- unique(
  rpe_limma_de_scd_significant,
  by = "Gene.symbol"
)
rpe_limma_de_scd_significant[
  ,
  regulation := ifelse(
      logFC < 0, "DOWN", "UP"
  )
]
table(rpe_limma_de_scd_significant$regulation)

fwrite(
    rpe_limma_de_scd_significant,
    paste0(
      path_results,
      "D1_st5.csv"
    )
)

#for CellChat (actually all included in scDiffCom already)
rpe_limma_de_cc <- rpe_limma_de[
  Gene.symbol %in% cc_genes,
  c("Gene.symbol", "logFC", "P.Value", "adj.P.Val")
]
rpe_limma_de_cc[
  ,
  bh_value := p.adjust(P.Value, method = "BH")
]
rpe_limma_de_cc_significant <- rpe_limma_de_cc[
  bh_value <= 0.1 & abs(logFC) >= log2(1.5) / 70
]
rpe_limma_de_cc_significant[
  ,
  pi_score_abs := -log10(bh_value) * abs(logFC)
]
rpe_limma_de_cc_significant <- rpe_limma_de_cc_significant[
  order(-pi_score_abs)
]
rpe_limma_de_cc_significant <- unique(
  rpe_limma_de_cc_significant,
  by = "Gene.symbol"
)
rpe_limma_de_cc_significant[
  ,
  regulation := ifelse(
      logFC < 0, "DOWN", "UP"
  )
]

table(
  rpe_limma_de_cc_significant$Gene.symbol %in%
  rpe_limma_de_scd_significant$Gene.symbol
)

## Find DEG also parts of scDiffCom CCIs, Table 1 ####

scd_cols_genes <- c(
  "LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3"
)
scd_detected_genes <- unique(
  c(
    unlist(
      cci_scd_healthy[, scd_cols_genes, with = FALSE]
    ),
     unlist(
      cci_scd_amd[, scd_cols_genes, with = FALSE]
    )
  )
)
scd_detected_genes <- scd_detected_genes[!is.na(scd_detected_genes)]

rpe_limma_de_icc <- rpe_limma_de_scd_significant[
  Gene.symbol %in% scd_detected_genes
]

#add detected LRIs
scd_lri_melt <- melt.data.table(
  LRI_human$LRI_curated[
    ,
    c("LRI", scd_cols_genes),
    with = FALSE
  ],
  id.vars = "LRI",
  value.name = "gene"
)
scd_lri_melt <- na.omit(scd_lri_melt)
rpe_limma_de_icc[
  dcast.data.table(
    unique(scd_lri_melt[, c("LRI", "gene")])[
      LRI %in% c(cci_scd_healthy$LRI, cci_scd_amd$LRI)
    ],
    formula = gene ~ .,
    value.var = "LRI",
    fun.aggregate = function(i) {
      paste(i, collapse = ",")
    }
  ),
  on = "Gene.symbol==gene",
  LRIs := i..
]

#add CellChat pathways

process_cldb_cellchat <- function(cc) {
  cct <- copy(cc)
  convert_table <- scDiffCom:::CellChat_conversion_human
  genes_to_change <- convert_table[new != "remove"]
  cct[, LIGAND_1 := sub(" - .*", "", interaction_name_2)]
  cct[, temp := sub(".* - ", "", interaction_name_2)]
  cct[, RECEPTOR_1 := ifelse(grepl("+", temp, fixed = TRUE),
                            gsub(".*\\((.+)\\+.*", "\\1", temp), temp)]
  cct[, RECEPTOR_2 := ifelse(grepl("+", temp, fixed = TRUE),
                            gsub(".*\\+(.+)\\).*", "\\1", temp), NA)]
  cct[, temp := NULL]
  cct[, LIGAND_1 := gsub(" ", "", LIGAND_1)]
  cct[, RECEPTOR_1 := gsub(" ", "", RECEPTOR_1)]
  cct[, RECEPTOR_2 := gsub(" ", "", RECEPTOR_2)]
  cct[genes_to_change,
     `:=`(LIGAND_1 = new),
     on = "LIGAND_1==old"
  ][
    genes_to_change,
    `:=`(RECEPTOR_1 = new),
    on = "RECEPTOR_1==old"
  ][
    genes_to_change,
    `:=`(RECEPTOR_2 = new),
    on = "RECEPTOR_2==old"
  ]
  cct[
    ,
    LRI_1 := ifelse(
      is.na(RECEPTOR_2),
      paste(LIGAND_1, RECEPTOR_1, sep = ":"),
      paste(LIGAND_1, paste(RECEPTOR_1, RECEPTOR_2, sep = "_"), sep = ":")
    )
  ]
  cct[
    ,
    LRI_2 := ifelse(
      is.na(RECEPTOR_2),
      paste(LIGAND_1, RECEPTOR_1, sep = ":"),
      paste(LIGAND_1, paste(RECEPTOR_2, RECEPTOR_1, sep = "_"), sep = ":")
    )
  ]
  cct[
    ,
    LRI_3 := ifelse(
      is.na(RECEPTOR_2),
      paste(RECEPTOR_1, LIGAND_1, sep = ":"),
      paste(LIGAND_1, paste(RECEPTOR_2, RECEPTOR_1, sep = "_"), sep = ":")
    )
  ]
  cct[
    ,
    LRI := ifelse(
      LRI_1 %in% LRI_human$LRI_curated$LRI,
      LRI_1,
      ifelse(
        LRI_2 %in% LRI_human$LRI_curated$LRI,
        LRI_2,
        ifelse(
          LRI_3 %in% LRI_human$LRI_curated$LRI,
          LRI_3,
          LRI_3
        )
      )
    )
  ]
  return(cct)
}
cldb_clean <- process_cldb_cellchat(cldb)

rpe_limma_de_pws <- unique(
  rpe_limma_de_icc[, c("Gene.symbol", "LRIs")][
    ,
    c(LRI = strsplit(LRIs, ",")),
    by = "Gene.symbol"
  ][
    cldb_clean,
    on = "LRI",
    pathway := i.pathway_name
  ][, c("Gene.symbol", "pathway")]
)
rpe_limma_de_pws <- na.omit(rpe_limma_de_pws)
rpe_limma_de_pws <- dcast.data.table(
  rpe_limma_de_pws,
  Gene.symbol ~ .,
  value.var = "pathway",
  fun.aggregate = function(i) {
    paste(i, collapse = "_")
  }
)

rpe_limma_de_icc[
  rpe_limma_de_pws,
  on = "Gene.symbol",
  pathways := i..
]
#note: some pathways will be further annotated manually

fwrite(
  rpe_limma_de_icc,
  paste0(
    path_results,
    "D1_t1_rpe_limma_de_icc.csv"
  )
)

## Regression plots for selected genes, Figure 4ACE ####

regp_genes <- c(
    "VEGFA", "KDR", "BMP7", "BMPR2", "TNXB", "SDC4"
)

regp_genes_dt <- as.data.table(
  pData(
    featureData(
      rpe_gset
    )
  )
)[
  Gene.symbol %in% c(
    "VEGFA", "KDR", "BMP7", "BMPR2", "TNXB", "SDC4"
  ),
  c("Gene.symbol", "ID")
]

regp_genes_dt <- regp_genes_dt[
  order(Gene.symbol, ID)
]
regp_genes_dt[
  ,
  uniqueID := paste(
    Gene.symbol,
    ID,
    sep = "_"
  )
]
regp_genes_dt[
  ,
  ymax := c(
    8, 6, 12, 12,
    8, 10, 6,
    10, 6,
    12, 6,
    12, 14, 6, 6, 6, 6,
    12, 8, 8, 10
  )
]

regp_results <- lapply(
  regp_genes_dt$ID,
  function(id) {
    temp_data <- data.table(
      age = rpe_gset$age,
      expr = exprs(rpe_gset)[as.character(id), ]
    )
    temp_lm <- lm(
      expr ~ age,
      data = temp_data
    )
    pval <- summary(temp_lm)$coefficients[2, 4]
    coef <- temp_lm$coefficients[[2]]
    temp_p <- ggplot(
      temp_data,
      aes(
        x = age,
        y = expr
      )
    ) + geom_point(
    ) + geom_smooth(
      method = "lm"
    ) + ylim(
      0,
      12
      #regp_genes_dt[ID == id]$ymax
    ) + ggtitle(
      paste(
        regp_genes_dt[ID == id]$Gene.symbol,
        regp_genes_dt[ID == id]$ID
      )
    )
    list(
      coef = coef,
      pval = pval,
      plot = temp_p
    )
  }
)
names(regp_results) <- regp_genes_dt$uniqueID

# manual visualization and selection

rpe_limma_de_icc[Gene.symbol == "BMP7"]
regp_results$BMP7_4447
regp_results$BMP7_8282
regp_results$BMP7_10756
regp_results$BMP7_15328 #to select

rpe_limma_de_icc[Gene.symbol == "BMPR2"]
regp_results$BMPR2_24626
regp_results$BMPR2_35882 #to select not sure
regp_results$BMPR2_43844

rpe_limma_de_icc[Gene.symbol == "KDR"]
regp_results$KDR_5091
regp_results$KDR_41403 #to select

rpe_limma_de_icc[Gene.symbol == "SDC4"]
regp_results$SDC4_5331 #to select
regp_results$SDC4_24315

rpe_limma_de_icc[Gene.symbol == "TNXB"]
regp_results$TNXB_17215
regp_results$TNXB_21280 #to select
regp_results$TNXB_27052
regp_results$TNXB_34185
regp_results$TNXB_43275
regp_results$TNXB_44221

rpe_limma_de_icc[Gene.symbol == "VEGFA"]
regp_results$VEGFA_10384
regp_results$VEGFA_25811 #to select
regp_results$VEGFA_27120
regp_results$VEGFA_36140

## Figure 4ACE ####

## Figure 4BDF ####

fig_4bdf <- VlnPlot(
  amd_seurat,
  features = regp_genes,
  pt.size = 0,
  ncol = 2
)
ggsave(
  paste0(
    path_results,
    "images/D1_f4bdf.png"
  ),
  fig_4bdf,
  width = 2000,
  height = 3000,
  units = "px"
)
