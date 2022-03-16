####################################################
##
## Project: amd_aging
##
## D.Dhirachaikulpanich@liverpool.ac.uk
## cyril.lagger@liverpool.ac.uk
##
####################################################
##


## Load bulk analysis (ba) data ####

ba_data <- fread("../data/D1_limma_data.csv")

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
    scd_genes %in% unique(ba_data$Gene.symbol)
)
ba_scd <- ba_data[Gene.symbol %in% scd_genes,
 c("Gene.symbol", "logFC", "P.Value", "adj.P.Val")]
ba_scd[, bh_value := p.adjust(P.Value, method = "BH")]
ba_scd_significant <- ba_scd[bh_value <= 0.1 & abs(logFC) >= log2(1.5)/70]

fwrite(
    ba_scd_significant,
    "../results/D1_bulk_aging_deg_scd.csv"
)

ba_scd_up <- unique(
    ba_scd_significant[logFC >= 0]$Gene.symbol
)
ba_scd_down <- unique(
    ba_scd_significant[logFC <= 0]$Gene.symbol
)

## Subset aging genes by CellChat ####

table(
    cc_genes %in% unique(ba_data$Gene.symbol)
)
ba_cc <- ba_data[Gene.symbol %in% cc_genes,
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

ba_scd_significant

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

ba_scd_significant_detected <- ba_scd_significant[Gene.symbol %in% icc_rpe_endo_genes]
fwrite(
    ba_scd_significant_detected,
    "../results/D1_icc_rpe_endo_bulk_aging_detected.csv"
)

## Subset aging genes by senescence ####



table(
    senref_kasit_up$gene %in% unique(ba_data$Gene.symbol)
)
ba_kasit_up <- ba_data[Gene.symbol %in% senref_kasit_up$gene,
 c("Gene.symbol", "logFC", "P.Value", "adj.P.Val")]
ba_kasit_up[, bh_value := p.adjust(P.Value, method = "BH")]