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

## Differential expression analysis ####

rpe_design <- model.matrix(~ age + sex + mac + disease + 0, rpe_gset)
rpe_fit <- lmFit(rpe_gset, rpe_design)
rpe_fit <- eBayes(rpe_fit, 0.01)
rpe_tT <- topTable(
  rpe_fit,
  adjust = "fdr",
  coef = "age",
  sort.by = "B",
  number = Inf
)

rpe_results <- decideTests(rpe_fit)
vennDiagram(rpe_results, include = c("up", "down"))

## Load pre-saved limma results ####

rpe_limma_de <- fread(
  paste0(
    path_data,
    "D1_limma_data.csv"
  )
)

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

## Subset aging genes by scDiffCom ####

table(
  scd_genes %in% unique(rpe_limma_de$Gene.symbol)
)
rpe_limma_de_scd <- rpe_limma_de[
  Gene.symbol %in% scd_genes,
  c("Gene.symbol", "logFC", "P.Value", "adj.P.Val")
]
rpe_limma_de_scd[, bh_value := p.adjust(P.Value, method = "BH")]
rpe_limma_de_scd_significant <- rpe_limma_de_scd[
  bh_value <= 0.1 & abs(logFC) >= log2(1.5) / 70
]

fwrite(
    rpe_limma_de_scd_significant,
    "../results/D1_bulk_aging_deg_scd.csv"
)

unique(rpe_limma_de_scd_significant$Gene.symbol)

rpe_limma_de_scd_up <- unique(
    rpe_limma_de_scd_significant[logFC >= 0]$Gene.symbol
)
rpe_limma_de_scd_down <- unique(
    rpe_limma_de_scd_significant[logFC <= 0]$Gene.symbol
)

table(
    unique(rpe_limma_de_scd_up) %in% unique(
        unlist(
            cci_scd_healthy[, c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]
        )
    ) | unique(rpe_limma_de_scd_up) %in% unique(
        unlist(
            cci_scd_amd[, c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]
        )
    )
)

table(
    unique(rpe_limma_de_scd_down) %in% unique(
        unlist(
            cci_scd_healthy[, c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]
        )
    ) | unique(rpe_limma_de_scd_down) %in% unique(
        unlist(
            cci_scd_amd[, c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]
        )
    )
)

ba_final_up <- sort(unique(rpe_limma_de_scd_up)[
    unique(rpe_limma_de_scd_up) %in% unique(
        unlist(
            cci_scd_healthy[, c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]
        )
    ) | unique(rpe_limma_de_scd_up) %in% unique(
        unlist(
            cci_scd_amd[, c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]
        )
    )
])

fwrite(
    rpe_limma_de_scd_significant[Gene.symbol %in% ba_final_up][, -4],
    "../../../../../deg_scd_up.csv"
)

ba_final_down <- sort(unique(rpe_limma_de_scd_down)[
    unique(rpe_limma_de_scd_down) %in% unique(
        unlist(
            cci_scd_healthy[, c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]
        )
    ) | unique(rpe_limma_de_scd_down) %in% unique(
        unlist(
            cci_scd_amd[, c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]
        )
    )
])

fwrite(
    rpe_limma_de_scd_significant[Gene.symbol %in% ba_final_down][, -4],
    "../../../../../deg_scd_down.csv"
)

ba_final_up %in% cci_cellchat_healthy$LIGAND_1
ba_final_up %in% cci_cellchat_healthy$RECEPTOR_2

## Subset aging genes by CellChat ####

table(
    cc_genes %in% unique(rpe_limma_de$Gene.symbol)
)
ba_cc <- rpe_limma_de[Gene.symbol %in% cc_genes,
                 c("Gene.symbol", "logFC", "P.Value", "adj.P.Val")]
ba_cc[, bh_value := p.adjust(P.Value, method = "BH")]
ba_cc_significant <- ba_cc[bh_value <= 0.1 & abs(logFC) >= log2(1.5)/70]

fwrite(
    ba_cc_significant,
    "../results/D1_bulk_aging_deg_cc.csv"
)

## Age regulated Ligand-receptor overlaping with scd results ####

icc_ba <- copy(icc_scdiffcom$icc_detec_all@cci_table_detected)
icc_ba[
    ,
    LIGAND_1_UP := ifelse(
        LIGAND_1 %in% ba_genes_up, TRUE, FALSE
    )
]
icc_ba[
    ,
    LIGAND_2_UP := ifelse(
        LIGAND_2 %in% ba_genes_up, TRUE, FALSE
    )
]
icc_ba[
    ,
    RECEPTOR_1_UP := ifelse(
        RECEPTOR_1 %in% ba_genes_up, TRUE, FALSE
    )
]
icc_ba[
    ,
    RECEPTOR_2_UP := ifelse(
        RECEPTOR_2 %in% ba_genes_up, TRUE, FALSE
    )
]
icc_ba[
    ,
    RECEPTOR_3_UP := ifelse(
        RECEPTOR_3 %in% ba_genes_up, TRUE, FALSE
    )
]

icc_ba[
    ,
    LIGAND_1_DOWN := ifelse(
        LIGAND_1 %in% ba_genes_down, TRUE, FALSE
    )
]
icc_ba[
    ,
    LIGAND_2_DOWN := ifelse(
        LIGAND_2 %in% ba_genes_down, TRUE, FALSE
    )
]
icc_ba[
    ,
    RECEPTOR_1_DOWN := ifelse(
        RECEPTOR_1 %in% ba_genes_down, TRUE, FALSE
    )
]
icc_ba[
    ,
    RECEPTOR_2_DOWN := ifelse(
        RECEPTOR_2 %in% ba_genes_down, TRUE, FALSE
    )
]
icc_ba[
    ,
    RECEPTOR_3_DOWN := ifelse(
        RECEPTOR_3 %in% ba_genes_down, TRUE, FALSE
    )
]

icc_ba_relevant <- icc_ba[
    ER_CELLTYPES %in%
        c("endothelial_RPE", "endothelial_endothelial",
          "RPE_RPE", "RPE_endothelial") &
        (LIGAND_1_UP | LIGAND_2_UP | LIGAND_1_DOWN | LIGAND_2_DOWN |
             RECEPTOR_1_UP | RECEPTOR_2_UP | RECEPTOR_3_UP |
             RECEPTOR_1_DOWN | RECEPTOR_2_DOWN | RECEPTOR_3_DOWN)
][, c("LRI", "ER_CELLTYPES", "CCI_SCORE",
      "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", 
      "LIGAND_1_UP", "LIGAND_2_UP" , "LIGAND_1_DOWN" , "LIGAND_2_DOWN" ,
      "RECEPTOR_1_UP" , "RECEPTOR_2_UP" , "RECEPTOR_3_UP" ,
      "RECEPTOR_1_DOWN" , "RECEPTOR_2_DOWN" , "RECEPTOR_3_DOWN"
)
]

fwrite(
    icc_ba_relevant,
    "../results/D1_icc_rpe_endo_bulk_aging.csv"
)

rpe_limma_de_scd_significant

icc_rpe_endo_genes <- unique(
    unlist(
        icc_ba[
            ER_CELLTYPES %in%
                c("endothelial_RPE", "endothelial_endothelial",
                  "RPE_RPE", "RPE_endothelial"), c(
                      "LIGAND_1", "LIGAND_2",
                      "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3"
                  )
        ]
    )
)
icc_rpe_endo_genes <- icc_rpe_endo_genes[!is.na(icc_rpe_endo_genes)]

rpe_limma_de_scd_significant_detected <- rpe_limma_de_scd_significant[Gene.symbol %in% icc_rpe_endo_genes]
fwrite(
    rpe_limma_de_scd_significant_detected,
    "../results/D1_icc_rpe_endo_bulk_aging_detected.csv"
)

## Subset aging genes by senescence ####



table(
    senref_kasit_up$gene %in% unique(rpe_limma_de$Gene.symbol)
)
ba_kasit_up <- rpe_limma_de[Gene.symbol %in% senref_kasit_up$gene,
                       c("Gene.symbol", "logFC", "P.Value", "adj.P.Val")]
ba_kasit_up[, bh_value := p.adjust(P.Value, method = "BH")]