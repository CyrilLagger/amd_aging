####################################################
##
## Project: amd_aging
##
## D.Dhirachaikulpanich@liverpool.ac.uk
## cyril.lagger@liverpool.ac.uk
##
####################################################
##

## Upset plot of detected CCIs (scDiffCom) ####

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

## Upset plot of detected CCIs (scDiffCom senescence) ####

icc_scd_ccis_sen <- data.table(
    cci = unique(
        c(
            icc_sen_scdiffcom$icc_detec_all@cci_table_detected$CCI,
            icc_sen_scdiffcom$icc_detec_54@cci_table_detected$CCI,
            icc_sen_scdiffcom$icc_detec_82@cci_table_detected$CCI,
            icc_sen_scdiffcom$icc_detec_79@cci_table_detected$CCI
        )
    )
)
icc_scd_ccis_sen[
    ,
    detec_all := ifelse(
        cci %in% icc_sen_scdiffcom$icc_detec_all@cci_table_detected$CCI,
        T,
        F
    )
]
icc_scd_ccis_sen[
    ,
    detec_54 := ifelse(
        cci %in% icc_sen_scdiffcom$icc_detec_54@cci_table_detected$CCI,
        T,
        F
    )
]
icc_scd_ccis_sen[
    ,
    detec_82 := ifelse(
        cci %in% icc_sen_scdiffcom$icc_detec_82@cci_table_detected$CCI,
        T,
        F
    )
]
icc_scd_ccis_sen[
    ,
    detec_79 := ifelse(
        cci %in% icc_sen_scdiffcom$icc_detec_79@cci_table_detected$CCI,
        T,
        F
    )
]

ComplexUpset::upset(
    icc_scd_ccis_sen,
    c("detec_all", "detec_54", "detec_82", "detec_79")
)
ComplexUpset::upset(
    icc_scd_ccis_sen,
    c("detec_54", "detec_82", "detec_79")
)

icc_scd_ccis_sen_54 <- icc_scd_ccis_sen[
    detec_54 == TRUE & detec_79 == FALSE  & detec_82 == FALSE
]
icc_scd_ccis_sen_54[
    ,
    c("E", "R", "LRI", "X", "Y") := tstrsplit(
        cci, "_"
    )
]
icc_scd_ccis_sen_54[, .N, by = c("E", "R")][order(-N)][1:10]
icc_scd_ccis_sen_54[, .N, by = c("LRI")][order(-N)][1:10]

## Compare scDiffCom to CellChat detection ####

setDT(df_sen_net)
df_sen_net
df_sen_net[
    ,
    ER := do.call(paste, c(.SD, sep = "_")),
    .SDcols = c("source", "target")
]
df_sen_net[
    ,
    LRI := do.call(paste, c(.SD, sep = ":")),
    .SDcols = c("ligand", "receptor")
]
df_sen_net[
    ,
    CCI := do.call(paste, c(.SD, sep = "_")),
    .SDcols = c("ER", "LRI")
]

intersect(
    df_sen_net$CCI,
    icc_sen_scdiffcom$icc_detec_all@cci_table_detected$CCI
)
# does not work because of gene names!!!


## Pathway of interest ####

cldb <- setDT(copy(CellChatDB$interaction))

vegf_pathway <- LRI_human$LRI_curated[
    LIGAND_1 %in% cldb[pathway_name == "VEGF"]$ligand |
    LIGAND_2 %in% cldb[pathway_name == "VEGF"]$ligand |
    RECEPTOR_1 %in% cldb[pathway_name == "VEGF"]$receptor |
    RECEPTOR_2 %in% cldb[pathway_name == "VEGF"]$receptor |
    RECEPTOR_3 %in% cldb[pathway_name == "VEGF"]$receptor 
]

bmp_pathway <- LRI_human$LRI_curated[
    LIGAND_1 %in% cldb[pathway_name == "BMP"]$ligand |
    LIGAND_2 %in% cldb[pathway_name == "BMP"]$ligand |
    RECEPTOR_1 %in% cldb[pathway_name == "BMP"]$receptor |
    RECEPTOR_2 %in% cldb[pathway_name == "BMP"]$receptor |
    RECEPTOR_3 %in% cldb[pathway_name == "BMP"]$receptor
]

tenascin_pathway <- LRI_human$LRI_curated[
    LIGAND_1 %in% cldb[pathway_name == "TENASCIN"]$ligand |
    LIGAND_2 %in% cldb[pathway_name == "TENASCIN"]$ligand #|
    #RECEPTOR_1 %in% cldb[pathway_name == "TENASCIN"]$receptor |
    #RECEPTOR_2 %in% cldb[pathway_name == "TENASCIN"]$receptor |
    #RECEPTOR_3 %in% cldb[pathway_name == "TENASCIN"]$receptor
]


## scDiffCom AMD exploration ####

table(icc_sen_scdiffcom$icc_diff_54_vs_79@cci_table_detected$REGULATION)
intersect(
    icc_scd_ccis_sen_54$cci,
    icc_sen_scdiffcom$icc_diff_54_vs_79@cci_table_detected[
        REGULATION == "DOWN"
    ]$CCI
)

BuildNetwork(
    icc_sen_scdiffcom$icc_diff_54_vs_79
)
PlotORA(
    icc_sen_scdiffcom$icc_diff_54_vs_79,
    category = "LRI",
    regulation = "UP"
)
PlotORA(
    icc_sen_scdiffcom$icc_diff_54_vs_79,
    category = "LRI",
    regulation = "DOWN"
)
PlotORA(
    icc_sen_scdiffcom$icc_diff_54_vs_79,
    category = "GO_TERMS",
    regulation = "UP"
)
PlotORA(
    icc_sen_scdiffcom$icc_diff_54_vs_79,
    category = "GO_TERMS",
    regulation = "DOWN"
)
PlotORA(
    icc_sen_scdiffcom$icc_diff_54_vs_79,
    category = "KEGG_PWS",
    regulation = "UP"
)
PlotORA(
    icc_sen_scdiffcom$icc_diff_54_vs_79,
    category = "KEGG_PWS",
    regulation = "DOWN"
)

icc_sen_scdiffcom$icc_diff_54_vs_79@cci_table_detected[
    LRI %in% vegf_pathway$LRI
][, .N, by = "ER_CELLTYPES"][order(-N)]

table(
    icc_sen_scdiffcom$icc_diff_54_vs_79@cci_table_detected[
    LRI %in% vegf_pathway$LRI
]$REGULATION
)

table(
    icc_sen_scdiffcom$icc_diff_54_vs_79@cci_table_detected$REGULATION
)

test <- icc_sen_scdiffcom$icc_diff_54_vs_79@cci_table_detected[
    LRI %in% vegf_pathway$LRI
][ , .N, by = c("REGULATION", "ER_CELLTYPES")]


#####



groupSize <- as.numeric(table(cellchat_sen@idents))
par(mfrow = c(1,2), xpd=TRUE)
jpeg("cellchat_1.jpg")
netVisual_circle(cellchat_sen@net$count,
 vertex.weight = groupSize, weight.scale = T,
  label.edge= F, title.name = "Number of interactions")
dev.off()

jpeg("cellchat_2.jpg")
netVisual_circle(cellchat_sen@net$weight, 
vertex.weight = groupSize, weight.scale = T, 
label.edge= F, 
title.name = "Interaction weights/strength")
dev.off()

jpeg("cellchat_3.jpg")
netVisual_bubble(cellchat_sen,
 sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
dev.off()

plotGeneExpression(cellchat_sen, signaling = "VEGF")
jpeg("cellchat_4.jpg")
pathways1.show <- c("VEGF") 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat_sen, signaling = pathways1.show,  vertex.receiver = vertex.receiver)
dev.off()

cellchat_sen@net


###

test <- icc_sen_scdiffcom$icc_diff_healthy_vs_amd@cci_table_detected
test <- GetTableCCI(icc_sen_scdiffcom$icc_diff_healthy_vs_amd, "detected", TRUE)

PlotORA(
    icc_sen_scdiffcom$icc_diff_healthy_vs_amd,
    category = "LRI"
)
BuildNetwork(
    icc_sen_scdiffcom$icc_diff_healthy_vs_amd
)

#####





## Explore senescence communication (scDiffCom) ####

cci_sen_all <- copy(icc_sen_scdiffcom$icc_detec_all@cci_table_detected)

vegf_pathway <- LRI_human$LRI_curated[
    LIGAND_1 %in% cldb[pathway_name == "VEGF"]$ligand
]

bmp_pathway <- LRI_human$LRI_curated[
    LIGAND_1 %in% cldb[pathway_name == "BMP"]$ligand
]

tenascin_pathway <- LRI_human$LRI_curated[
    LIGAND_1 %in% cldb[pathway_name == "TENASCIN"]$ligand
]

cldb_pathways <- lapply(
    unique(cldb$pathway_name),
    function(i) {
        LRI_human$LRI_curated[
            LIGAND_1 %in% cldb[pathway_name == i]$ligand
        ]
    }
)


cci_sen_all[LRI %in% vegf_pathway$LRI]

cci_sen_vegf <- cci_sen_all[LRI %in% vegf_pathway$LRI]
cci_sen_bmp <- cci_sen_all[LRI %in% bmp_pathway$LRI]
cci_sen_tenascin <- cci_sen_all[LRI %in% tenascin_pathway$LRI]

test <- cci_sen_tenascin[
    ,
    sum(CCI_SCORE),
    by = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE")
]
setnames(test, old = "V1", new = "weight")

test <- cci_sen_all[LRI %in% vegf_pathway$LRI][
    ,
    .N,
    by = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE")
]
setnames(test, old = "N", new = "weight")

g0 <- graph_from_data_frame(
    test, directed = TRUE
)
coords <- layout_in_circle(g0)
get.adjacency(g0)
E(g0)$weight
is_weighted(g0)

g0_ecentraliy <- eigen_centrality(
    g0,
     directed = FALSE,
     weights = E(g0)$weighted
)
sort(g0_ecentraliy$vector, decreasing = TRUE)


edge_connectivity()

jpeg("g0.jpg")
plot(g0,
  vertex.size=3,
  edge.width = E(g0)$weight,
   #edge.width=ifelse(E(g0)$weight > 8, 0.5*E(g0)$weight ,0.1),
   layout = coords
)
dev.off()




gsub(
    "_", ":",
    cldb[pathway_name == "VEGF"]$interaction_name,
) %in% LRI_human$LRI_curated$LRI

LRI_human$LRI_curated[LIGAND_2 == "VEGFA"][, 1:4]


paste(
    cldb[pathway_name == "VEGF"]$ligand,
    cldb[pathway_name == "VEGF"]$receptor,
    sep = ":"
) %in% LRI_human$LRI_curated$LRI


###

diff_sen <- icc_sen_scdiffcom$icc_diff_54_vs_79@cci_table_detected

diff_sen_up <- diff_sen[REGULATION == "UP"]

g <- PlotORA(
    icc_sen_scdiffcom$icc_diff_82_vs_79,
    category = "LRI",
    regulation = "DOWN"
)
ggsave("scd.png", g)
BuildNetwork(icc_sen_scdiffcom$icc_diff_54_vs_79)

table(
    diff_sen[ER_CELLTYPES == "melanocyte_RPE-sen"]$REGULATION
)

diff_sen[ER_CELLTYPES == "macrophage_macrophage"
 & REGULATION == "UP"][, c("LRI")]

table(diff_sen$REGULATION)

table(
    diff_sen[LRI %in% vegf_pathway$LRI]$REGULATION
)

diff_sen[LRI %in% vegf_pathway$LRI &
 REGULATION == "UP"][, c("CCI")]

diff_sen[LRI %in% tenascin_pathway$LRI &
 REGULATION == "UP"][, c("CCI", "LOGFC")]

###


pathways1.show <- c("VEGF")
pathways2.show <- c("BMP")
pathways3.show <- c("TENASCIN")

par(mfrow=c(1,1))
netVisual_heatmap(
    cellchat_sen,
     signaling = pathways3.show,
      color.heatmap = "Reds",
      font.size = 10
)

cellchat_sen <- netAnalysis_computeCentrality(cellchat_sen, slot.name = "netP")
netAnalysis_signalingRole_network(
  cellchat_sen,
  signaling = pathways1.show, width = 8, height = 2.5, font.size = 10)
