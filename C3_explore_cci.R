####################################################
##
## Project: amd_aging
##
## D.Dhirachaikulpanich@liverpool.ac.uk
## cyril.lagger@liverpool.ac.uk
##
####################################################
##

## Compare scDiffCom to CellChat and healthy to AMD detection ####

LRI_common <- LRI_human$LRI_curated[grepl("CellChat", DATABASE)]$LRI

process_cci_cellchat <- function(cc) {
  convert_table <- scDiffCom:::CellChat_conversion_human
  genes_to_change <- convert_table[new != "remove"]
  cc[, LIGAND_1 := sub(" - .*", "", interaction_name_2)]
  cc[, temp := sub(".* - ", "", interaction_name_2)]
  cc[, RECEPTOR_1 := ifelse(grepl("+", temp, fixed = TRUE),
                            gsub(".*\\((.+)\\+.*", "\\1", temp), temp)]
  cc[, RECEPTOR_2 := ifelse(grepl("+", temp, fixed = TRUE),
                            gsub(".*\\+(.+)\\).*", "\\1", temp), NA)]
  cc[, temp := NULL]
  cc[, LIGAND_1 := gsub(" ", "", LIGAND_1)]
  cc[, RECEPTOR_1 := gsub(" ", "", RECEPTOR_1)]
  cc[, RECEPTOR_2 := gsub(" ", "", RECEPTOR_2)]
  cc[genes_to_change,
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
  cc[
    ,
    LRI_1 := ifelse(
      is.na(RECEPTOR_2),
      paste(LIGAND_1, RECEPTOR_1, sep = ":"),
      paste(LIGAND_1, paste(RECEPTOR_1, RECEPTOR_2, sep = "_"), sep = ":")
    )
  ]
  cc[
    ,
    LRI_2 := ifelse(
      is.na(RECEPTOR_2),
      paste(LIGAND_1, RECEPTOR_1, sep = ":"),
      paste(LIGAND_1, paste(RECEPTOR_2, RECEPTOR_1, sep = "_"), sep = ":")
    )
  ]
  cc[
    ,
    LRI_3 := ifelse(
      is.na(RECEPTOR_2),
      paste(RECEPTOR_1, LIGAND_1, sep = ":"),
      paste(LIGAND_1, paste(RECEPTOR_2, RECEPTOR_1, sep = "_"), sep = ":")
    )
  ]
  cc[
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
  cc[
    ,
    CCI := paste(
      paste(
        source,
        target,
        sep = "_"
      ),
      LRI,
      sep = "_"
    )
  ]
  return(cc)
}

cci_cellchat_healthy_clean <- process_cci_cellchat(cci_cellchat_healthy)
cci_cellchat_amd_clean <- process_cci_cellchat(cci_cellchat_amd)

cci_cellchat_healthy_clean[
  ,
  CCI := ifelse(
    duplicated(CCI),
    paste(CCI, "dup", sep = "_"),
    CCI
  )
]
table(grepl("_dup", cci_cellchat_healthy_clean$CCI))

cci_cellchat_amd_clean[
  ,
  CCI := ifelse(
    duplicated(CCI),
    paste(CCI, "dup", sep = "_"),
    CCI
  )
]
table(grepl("_dup", cci_cellchat_amd_clean$CCI))

## Supplementary Tables 1-4 #####

cci_cellchat_healthy_clean[
  ,
  in_scDiffCom := CCI %in% scd_healthy@cci_table_detected$CCI
]
table(cci_cellchat_healthy_clean$in_scDiffCom)

st1 <- cci_cellchat_healthy_clean[, c(1:11, 20)][order(-prob)]
fwrite(
  st1,
  paste0(
    path_results,
    "C3_st1_cc_healthy.csv"
  )
)

cci_cellchat_amd_clean[
  ,
  in_scDiffCom := CCI %in% scd_amd@cci_table_detected$CCI
]
table(cci_cellchat_amd_clean$in_scDiffCom)
st2 <- cci_cellchat_amd_clean[, c(1:11, 20)][order(-prob)]
fwrite(
  st2,
  paste0(
    path_results,
    "C3_st2_cc_amd.csv"
  )
)

cci_scd_healthy <- scd_healthy@cci_table_detected
cci_scd_amd <- scd_amd@cci_table_detected

cci_scd_healthy[
  ,
  in_CellChat := ifelse(
    CCI %in% cci_cellchat_healthy_clean$CCI,
    TRUE,
    FALSE
    )
]
table(cci_scd_healthy$in_CellChat)
cci_scd_healthy[
  LRI_human$LRI_curated,
  on = "LRI",
  LRI_origin := i.DATABASE
]
st3 <- cci_scd_healthy[, c(3:5, 23:25, 31, 30)][order(-CCI_SCORE)]
fwrite(
  st3,
  paste0(
    path_results,
    "C3_st3_scd_healthy.csv"
  )
)

cci_scd_amd[
  ,
  in_CellChat := ifelse(
    CCI %in% cci_cellchat_amd_clean$CCI, TRUE, FALSE
    )
]
table(cci_scd_amd$in_CellChat)
cci_scd_amd[
  LRI_human$LRI_curated,
  on = "LRI",
  LRI_origin := i.DATABASE
]
st4 <- cci_scd_amd[, c(3:5, 23:25, 31, 30)][order(-CCI_SCORE)]
fwrite(
  st4,
  paste0(
    path_results,
    "C3_st4_scd_amd.csv"
  )
)

## Supplementary Figure 1A-C ####

overlap_venn_healthy <- ggVennDiagram::process_data(
  ggVennDiagram::Venn(
    list(
      "CellChat" = cci_cellchat_healthy_clean$CCI,
      "scDiffCom" = scd_healthy@cci_table_detected[
        LRI %in% LRI_common |
          CCI %in% cci_cellchat_healthy_clean$CCI]$CCI
    )
  )
)
overlap_venn_healthy  <- ggplot() +
  geom_sf(
    aes(fill = count),
    data = venn_region(overlap_venn_healthy),
    show.legend = FALSE
  ) +
  geom_sf(
    aes(color = name),
    data = venn_setedge(overlap_venn_healthy),
    show.legend = FALSE,
    size = 2
  ) +
  geom_sf_text(
    aes(label = name),
    data = venn_setlabel(overlap_venn_healthy),
    size = 6
  ) +
  geom_sf_label(
    aes(
      label = paste0(
        count,
        " (",
        scales::percent(count / sum(count), accuracy = 2),
        ")"
      )
    ),
    data = venn_region(overlap_venn_healthy),
    size = 6
  ) +
  scale_fill_gradient(
    low = "#F4FAFE",
    high = "#4981BF"
  ) +
  scale_color_manual(
    values = c("CellChat" = "grey", "scDiffCom" = "grey")
  ) +
  theme_void() #+
# ggtitle(
#   "Normal"
# ) +
# theme(
#   plot.title = element_text(size = 20)
# )
overlap_venn_healthy

overlap_venn_amd <- ggVennDiagram::process_data(
  ggVennDiagram::Venn(
    list(
      "CellChat" = cci_cellchat_amd_clean$CCI,
      "scDiffCom" = scd_amd@cci_table_detected[
        LRI %in% LRI_common |
          CCI %in% cci_cellchat_amd_clean$CCI]$CCI
    )
  )
)
overlap_venn_amd  <- ggplot() +
  geom_sf(
    aes(fill = count),
    data = venn_region(overlap_venn_amd),
    show.legend = FALSE
  ) +
  geom_sf(
    aes(color = name),
    data = venn_setedge(overlap_venn_amd),
    show.legend = FALSE,
    size = 2
  ) +
  geom_sf_text(
    aes(label = name),
    data = venn_setlabel(overlap_venn_amd),
    size = 6
  ) +
  geom_sf_label(
    aes(
      label = paste0(
        count,
        " (",
        scales::percent(count / sum(count), accuracy = 2),
        ")"
      )
    ),
    data = venn_region(overlap_venn_amd),
    size = 6
  ) +
  scale_fill_gradient(
    low = "#F4FAFE",
    high = "#4981BF"
  ) +
  scale_color_manual(
    values = c("CellChat" = "grey", "scDiffCom" = "grey")
  ) +
  theme_void()# +
# ggtitle(
#   "AMD"
# ) +
# theme(
#   plot.title = element_text(size = 20)
#)
overlap_venn_amd

overlap_venn_combined <- ggVennDiagram::process_data(
  ggVennDiagram::Venn(
    list(
      "CellChat Normal" = cci_cellchat_healthy_clean$CCI,
      "CellChat AMD" = cci_cellchat_amd_clean$CCI,
      "scDiffCom Normal" = scd_healthy@cci_table_detected[
        LRI %in% LRI_common |
          CCI %in% cci_cellchat_healthy_clean$CCI]$CCI,
      "scDiffCom AMD" = scd_amd@cci_table_detected[
        LRI %in% LRI_common |
          CCI %in% cci_cellchat_amd_clean$CCI]$CCI
    )
  )
)
overlap_venn_combined  <- ggplot() +
  geom_sf(
    aes(fill = count),
    data = venn_region(overlap_venn_combined),
    show.legend = FALSE
  ) +
  geom_sf(
    aes(color = name),
    data = venn_setedge(overlap_venn_combined),
    show.legend = FALSE,
    size = 2
  ) +
  geom_sf_text(
    aes(label = name),
    data = venn_setlabel(overlap_venn_combined),
    size = 6,
    nudge_x = c(0.05, 0, 0, -0.05)
  ) +
  geom_sf_label(
    aes(
      label = paste0(
        count,
        " (",
        scales::percent(count / sum(count), accuracy = 2),
        ")"
      )
    ),
    data = venn_region(overlap_venn_combined),
    size = 5
  ) +
  scale_fill_gradient(
    low = "#F4FAFE",
    high = "#4981BF"
  ) +
  scale_color_manual(
    values = c("CellChat" = "grey", "scDiffCom" = "grey")
  ) +
  theme_void()# +
# ggtitle(
#   "Overall"
# ) +
# theme(
#   plot.title = element_text(size = 20)
# )
overlap_venn_combined

overlap_venn_grid <- plot_grid(
  plot_grid(
    overlap_venn_healthy,
    overlap_venn_amd,
    ncol = 1,
    labels = c("a - Normal", "b - AMD")
  ),
  overlap_venn_combined,
  ncol = 2,
  labels = c("", "c - Overall")
)

overlap_venn_title <- ggdraw() +
  draw_label(
    "Overlap of detected CCIs",
    fontface = "bold",
    x = 0,
    hjust = 0,
    size = 18
  ) + theme(
    plot.margin = margin(0,0,0,7)
  )

sf1 <- plot_grid(
  overlap_venn_title,
  overlap_venn_grid,
  ncol = 1,
  rel_heights = c(0.1, 1)
)
ggsave(
  paste0(
    path_results,
    "images/C3_sf1.png"
  ),
  sf1,
  scale = 2
)

## CCIs presents in the 4-overlap ####

cci_common <- cci_cellchat_healthy_clean[
  CCI %in% cci_cellchat_amd_clean$CCI &
    CCI %in% cci_scd_healthy$CCI &
    CCI %in% cci_scd_amd$CCI
]

## Prepare a merged CellChat object ####

object.list <- list(Normal = cellchat_healthy, AMD = cellchat_amd)
cellchatmerge <- mergeCellChat(object.list, add.names = names(object.list))

## Find dominant patterns from/to RPE #####

rpe_sender_pathways <- t(sapply(
  cellchat_healthy@netP$centr,
  function(i) {
    c(i$outdeg[["RPE"]], mean(i$outdeg), sum(i$outdeg > i$outdeg[["RPE"]]))
  }
))

rpe_receiver_pathways <- t(sapply(
  cellchat_healthy@netP$centr,
  function(i) {
    c(i$indeg[["RPE"]], mean(i$indeg), sum(i$indeg > i$indeg[["RPE"]]))
  }
))

rpe_mediator_pathways <- t(sapply(
  cellchat_healthy@netP$centr,
  function(i) {
    c(i$flowbet[[5]], mean(i$flowbet), sum(i$flowbet > i$flowbet[[5]]))
  }
))

rpe_influencer_pathways <- t(sapply(
  cellchat_healthy@netP$centr,
  function(i) {
    c(i$info[[5]], mean(i$info), sum(i$info > i$info[[5]]))
  }
))

netAnalysis_signalingRole_network(
  cellchat_healthy,
  signaling = c("AGRN"),
  width = 8,
  height = 2.5,
  font.size = 10
)

## Secreted pathways heatmaps Figure 2A ####

png(
  paste0(
    path_results,
    "images/C3_f2a.png"
  ),
  width = 600
)
draw(
   netAnalysis_signalingRole_heatmap(
    object.list[[1]],
    signaling = c("VEGF","BMP","PTN","GDF","GRN","HGF"),
    pattern = "all",
    width = 8,
    height = 10,
    title = names(object.list)[1],
    color.heatmap = "OrRd"
  ) + 
    netAnalysis_signalingRole_heatmap(
    object.list[[2]],
    signaling = c("VEGF","BMP","PTN","GDF","GRN","HGF"),
    pattern = "all",
    width = 8,
    height = 10,
    title = names(object.list)[2],
    color.heatmap = "OrRd"
  ),
  ht_gap = unit(0.5, "cm")
)
dev.off()

## VEGF and BMP pathways heatmaps Fig 2B-C ####

png(
  paste0(
    path_results,
    "images/C3_f2b.png"
  ),
  width = 800
)
draw(
  netVisual_heatmap(
    cellchat_healthy,
    signaling = "VEGF",
    color.heatmap = "Reds",
    title.name = paste("VEGF signaling ", names(cc_list)[1]),
    width = 10,
    height = 10,
    font.size = 12
  ) +
    netVisual_heatmap(
      cellchat_amd,
      signaling = "VEGF",
      color.heatmap = "Reds",
      title.name = paste("VEGF signaling ", names(cc_list)[2]),
      width = 10,
      height = 10,
      font.size = 12
    ),
  ht_gap = unit(0.5, "cm")
)
dev.off()

png(
  paste0(
    path_results,
    "images/C3_f2c.png"
  ),
  width = 800
)
draw(
  netVisual_heatmap(
    cellchat_healthy,
    signaling = "BMP",
    color.heatmap = "Reds",
    title.name = paste("BMP signaling ", names(cc_list)[1]),
    width = 10,
    height = 10,
    font.size = 12
  ) +
    netVisual_heatmap(
      cellchat_amd,
      signaling = "BMP",
      color.heatmap = "Reds",
      title.name = paste("BMP signaling ", names(cc_list)[2]),
      width = 10,
      height = 10,
      font.size = 12
    ),
  ht_gap = unit(0.5, "cm")
)
dev.off()

## Supplementary Figures 2A-F ####

sf2_pathways <- data.table(
  id = c(
    "a1", "a2", "b1", "b2", "c1", "c2",
    "d1", "d2", "e1", "e2", "f1"
  ),
  pw = c(
    "VEGF", "VEGF", "BMP", "BMP", "PTN", "PTN",
    "GDF", "GDF", "GRN", "GRN", "HGF"
  ),
  condition = c(
    "healthy", "amd", "healthy", "amd", "healthy", "amd",
    "healthy", "amd", "healthy", "amd", "healthy" 
  )
)

for (i in seq_len(nrow(sf2_pathways))) {
  png(
    paste0(
      path_results,
      "images/C3_sf2",
      sf2_pathways[i]$id,
      ".png"
    ),
    width = 400,
    height = 200
  )
  netAnalysis_signalingRole_network(
    get(paste0("cellchat_", sf2_pathways[i]$condition)),
    signaling = sf2_pathways[i]$pw,
    width = 10,
    height = 3,
    font.size = 14,
    font.size.title = 15
  )
  dev.off()
}

## ECM interactions ####

cci_cellchat_healthy_clean[annotation == "ECM-Receptor"]
cci_cellchat_amd_clean[annotation == "ECM-Receptor"]

cldb[annotation == "ECM-Receptor"]

## ECM heatmaps ####

heat_healthy_ecm <- netAnalysis_signalingRole_heatmap(
  object.list[[1]],
  signaling = c("COLLAGEN","LAMININ","FN1","THBS","TENASCIN","VTN"),
  pattern = "all",
  width = 6,
  height = 7,
  title = names(object.list)[1],
  color.heatmap = "OrRd"
)
heat_amd_ecm <- netAnalysis_signalingRole_heatmap(
  object.list[[2]],
  signaling =  c("COLLAGEN","LAMININ","FN1","THBS","TENASCIN","VTN"),
  pattern = "all",
  width = 6,
  height = 7,
  title = names(object.list)[2],
  color.heatmap = "OrRd"
)
draw(heat_healthy_ecm + heat_amd_ecm, ht_gap = unit(0.5, "cm"))

draw(
  netVisual_heatmap(
    cellchat_healthy,
    signaling = "VEGF",
    color.heatmap = "Reds",
    title.name = paste("VEGF signaling ", names(object.list)[1]),
    width = 10,
    height = 10,
    font.size = 12
  )+
    netVisual_heatmap(
      cellchat_amd,
      signaling = "VEGF",
      color.heatmap = "Reds",
      title.name = paste("VEGF signaling ", names(object.list)[2]),
      width = 10,
      height = 10,
      font.size = 12
    ),
  ht_gap = unit(0.5, "cm")
)

draw(
  netVisual_heatmap(
    cellchat_healthy,
    signaling = "COLLAGEN",
    color.heatmap = "Reds",
    title.name = paste("COLLAGEN signaling ", names(object.list)[1]),
    width = 10,
    height = 10,
    font.size = 12
  )+
    netVisual_heatmap(
      cellchat_amd,
      signaling = "COLLAGEN",
      color.heatmap = "Reds",
      title.name = paste("COLLAGEN signaling ", names(object.list)[2]),
      width = 10,
      height = 10,
      font.size = 12
    ),
  ht_gap = unit(0.5, "cm")
)

draw(
  netVisual_heatmap(
    cellchat_healthy,
    signaling = "THBS",
    color.heatmap = "Reds",
    title.name = paste("THBS signaling ", names(object.list)[1]),
    width = 10,
    height = 10,
    font.size = 12
  )+
    netVisual_heatmap(
      cellchat_amd,
      signaling = "THBS",
      color.heatmap = "Reds",
      title.name = paste("THBS signaling ", names(object.list)[2]),
      width = 10,
      height = 10,
      font.size = 12
    ),
  ht_gap = unit(0.5, "cm")
)

draw(
  netVisual_heatmap(
    cellchat_healthy,
    signaling = "TENASCIN",
    color.heatmap = "Reds",
    title.name = paste("TENASCIN signaling ", names(object.list)[1]),
    width = 10,
    height = 10,
    font.size = 12
  )+
    netVisual_heatmap(
      cellchat_amd,
      signaling = "TENASCIN",
      color.heatmap = "Reds",
      title.name = paste("TENASCIN signaling ", names(object.list)[2]),
      width = 10,
      height = 10,
      font.size = 12
    ),
  ht_gap = unit(0.5, "cm")
)

netVisual_heatmap(
  cellchat_healthy,
  signaling = "COLLAGEN",
  color.heatmap = "Reds",
  title.name = paste("COLLAGEN signaling ", names(object.list)[1])
)

netVisual_heatmap(
  cellchat_amd,
  signaling = "COLLAGEN",
  color.heatmap = "Reds",
  title.name = paste("COLLAGEN signaling ", names(object.list)[2])
)

netAnalysis_signalingRole_network(
  cellchat_healthy,
  signaling = c("COLLAGEN"),
  width = 8,
  height = 2.5,
  font.size = 10
)

netVisual_heatmap(
  cellchat_healthy,
  signaling = "LAMININ",
  color.heatmap = "Reds",
  title.name = paste("LAMININ signaling ", names(object.list)[1])
)

netVisual_heatmap(
  cellchat_amd,
  signaling = "LAMININ",
  color.heatmap = "Reds",
  title.name = paste("LAMININ signaling ", names(object.list)[2])
)

netAnalysis_signalingRole_network(
  cellchat_healthy,
  signaling = c("LAMININ"),
  width = 8,
  height = 2.5,
  font.size = 10
)

netVisual_heatmap(
  cellchat_healthy,
  signaling = "TENASCIN",
  color.heatmap = "Reds",
  title.name = paste("TENASCIN signaling ", names(object.list)[1])
)

netVisual_heatmap(
  cellchat_amd,
  signaling = "LAMININ",
  color.heatmap = "Reds",
  title.name = paste("LAMININ signaling ", names(object.list)[2])
)

## Pathway of interest for scDiffCom ####

cldb <- setDT(copy(CellChatDB$interaction))

vegf_pathway <- LRI_human$LRI_curated[
  LIGAND_1 %in% cldb[pathway_name == "VEGF"]$ligand |
    LIGAND_2 %in% cldb[pathway_name == "VEGF"]$ligand |
    RECEPTOR_1 %in% cldb[pathway_name == "VEGF"]$receptor |
    RECEPTOR_2 %in% cldb[pathway_name == "VEGF"]$receptor |
    RECEPTOR_3 %in% cldb[pathway_name == "VEGF"]$receptor 
]

vegf_pathway2 <- LRI_human$LRI_curated[
  LIGAND_1 %in% cldb[pathway_name == "VEGF"]$ligand |
    LIGAND_2 %in% cldb[pathway_name == "VEGF"]$ligand
]


bmp_pathway <- LRI_human$LRI_curated[
  LIGAND_1 %in% cldb[pathway_name == "BMP"]$ligand |
    LIGAND_2 %in% cldb[pathway_name == "BMP"]$ligand |
    RECEPTOR_1 %in% cldb[pathway_name == "BMP"]$receptor |
    RECEPTOR_2 %in% cldb[pathway_name == "BMP"]$receptor |
    RECEPTOR_3 %in% cldb[pathway_name == "BMP"]$receptor
]

bmp_pathway2 <- LRI_human$LRI_curated[
  LIGAND_1 %in% cldb[pathway_name == "BMP"]$ligand |
    LIGAND_2 %in% cldb[pathway_name == "BMP"]$ligand 
]

collagen_pathway <- LRI_human$LRI_curated[
  LIGAND_1 %in% cldb[pathway_name == "COLLAGEN"]$ligand |
    LIGAND_2 %in% cldb[pathway_name == "COLLAGEN"]$ligand# |
  #RECEPTOR_1 %in% cldb[pathway_name == "COLLAGEN"]$receptor |
  #RECEPTOR_2 %in% cldb[pathway_name == "COLLAGEN"]$receptor |
  #RECEPTOR_3 %in% cldb[pathway_name == "COLLAGEN"]$receptor 
]

thbs_pathway <- LRI_human$LRI_curated[
  LIGAND_1 %in% cldb[pathway_name == "THBS"]$ligand |
    LIGAND_2 %in% cldb[pathway_name == "THBS"]$ligand #|
  #RECEPTOR_1 %in% cldb[pathway_name == "THBS"]$receptor |
  #RECEPTOR_2 %in% cldb[pathway_name == "THBS"]$receptor |
  #RECEPTOR_3 %in% cldb[pathway_name == "THBS"]$receptor 
]


tenascin_pathway <- LRI_human$LRI_curated[
  LIGAND_1 %in% cldb[pathway_name == "TENASCIN"]$ligand |
    LIGAND_2 %in% cldb[pathway_name == "TENASCIN"]$ligand #|
  #RECEPTOR_1 %in% cldb[pathway_name == "TENASCIN"]$receptor |
  #RECEPTOR_2 %in% cldb[pathway_name == "TENASCIN"]$receptor |
  #RECEPTOR_3 %in% cldb[pathway_name == "TENASCIN"]$receptor
]

## scDiffCom heatmaps for specific pathways ####

scd_heatmap <- function(cci_dt, title.name) {
  scd_template <- CJ(
    factor(levels(cellchat_healthy@idents), levels(cellchat_healthy@idents)),
    factor(levels(cellchat_healthy@idents), levels(cellchat_healthy@idents)),
    sorted = FALSE
  )
  setnames(
    scd_template,
    old = colnames(scd_template),
    new = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE")
  )
  dt <- copy(scd_template)[
    cci_dt[
      ,
      sum(CCI_SCORE),
      by = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE")
    ],
    on = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE"),
    score_sum := i.V1
  ]
  dt[is.na(dt)] <- 0
  mat <- reshape2::acast(
    dt,
    EMITTER_CELLTYPE ~ RECEIVER_CELLTYPE,
    value.var = "score_sum"
  )
  
  color.use <- scPalette(ncol(mat))
  names(color.use) <- colnames(mat)
  color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = "Reds")))(100)
  
  df<- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),which = "column",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), which = "row",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
  mat[mat == 0] <- NA
  Heatmap(
    mat,
    col = color.heatmap.use,
    na_col = "white",
    name = "CCI Score Sum",
    bottom_annotation = col_annotation,
    left_annotation = row_annotation,
    top_annotation = ha2,
    right_annotation = ha1,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    row_names_rot = 0,
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    column_title = title.name,
    column_title_gp = gpar(fontsize = 10),
    column_names_rot = 90,
    row_title = "Sources (Sender)",
    row_title_gp = gpar(fontsize = 10),
    row_title_rot = 90,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 8, fontface = "plain"),
      title_position = "leftcenter-rot",
      border = NA, #at = colorbar.break,
      legend_height = unit(20, "mm"),
      labels_gp = gpar(fontsize = 8),
      grid_width = unit(2, "mm")
    )
  )
}

draw(
  scd_heatmap(cci_scd_healthy[LRI %in% vegf_pathway$LRI], "VEGF signaling network - Normal")+
    scd_heatmap(cci_scd_amd[LRI %in% vegf_pathway$LRI], "VEGF signaling network - AMD"),
  ht_gap = unit(0.5, "cm")
)
draw(
  scd_heatmap(cci_scd_healthy[LRI %in% bmp_pathway$LRI], "BMP signaling network - Normal")+
    scd_heatmap(cci_scd_amd[LRI %in% bmp_pathway$LRI], "BMP signaling network - AMD"),
  ht_gap = unit(0.5, "cm")
)

draw(
  scd_heatmap(cci_scd_healthy[LRI %in% collagen_pathway$LRI], "COLLAGEN signaling network - Normal")+
    scd_heatmap(cci_scd_amd[LRI %in% collagen_pathway$LRI], "COLLAGEN signaling network - AMD"),
  ht_gap = unit(0.5, "cm")
)

draw(
  scd_heatmap(cci_scd_healthy[LRI %in% thbs_pathway$LRI], "THBS signaling network - Normal")+
    scd_heatmap(cci_scd_amd[LRI %in% thbs_pathway$LRI], "THBS signaling network - AMD"),
  ht_gap = unit(0.5, "cm")
)

draw(
  scd_heatmap(cci_scd_healthy[LRI %in% tenascin_pathway$LRI], "TENASCIN signaling network - Normal")+
    scd_heatmap(cci_scd_amd[LRI %in% tenascin_pathway$LRI], "TENASCIN signaling network - AMD"),
  ht_gap = unit(0.5, "cm")
)


cci_ER_tenascin <- cci_sen_tenascin[
  ,
  sum(CCI_SCORE),
  by = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE")
][order(-V1)]
setnames(cci_ER_tenascin, old = "V1", new = "weight")
cci_ER_tenascin
cci_ER_tenascin_2 <- dcast.data.table(
  cci_ER_tenascin,
  formula = EMITTER_CELLTYPE ~ RECEIVER_CELLTYPE,
  value.var = "weight",
  fill = 0
)
scd_tenascin_mat <- as.matrix(cci_ER_tenascin_2[,-1])
rownames(scd_tenascin_mat) <- cci_ER_tenascin_2[, 1]$EMITTER_CELLTYPE
heatmap(scd_tenascin_mat, scale = "none")

## Compare scDiffCom to CellChat for senescence detection ####

setDT(cci_cellchat_sen)
scd_sen

cci_cellchat_sen <- process_cci_cellchat(cci_cellchat_sen)

cci_cellchat_sen[
  ,
  CCI := ifelse(
    duplicated(CCI),
    paste(CCI, "dup", sep = "_"),
    CCI
  )
]
table(grepl("_dup", cci_cellchat_sen$CCI))

## Senescence Venn Diagram ####

overlap_venn_sen <- ggVennDiagram::process_data(
  ggVennDiagram::Venn(
    list(
      "CellChat" = cci_cellchat_sen$CCI,
      "scDiffCom" = scd_sen@cci_table_detected[
        LRI %in% LRI_common |
          CCI %in% cci_cellchat_sen$CCI]$CCI
    )
  )
)
overlap_venn_sen  <- ggplot() +
  geom_sf(
    aes(fill = count),
    data = venn_region(overlap_venn_sen),
    show.legend = FALSE
  ) +
  geom_sf(
    aes(color = name),
    data = venn_setedge(overlap_venn_sen),
    show.legend = FALSE,
    size = 2
  ) +
  geom_sf_text(
    aes(label = name),
    data = venn_setlabel(overlap_venn_sen),
    size = 6
  ) +
  geom_sf_label(
    aes(
      label = paste0(
        count,
        " (",
        scales::percent(count / sum(count), accuracy = 2),
        ")"
      )
    ),
    data = venn_region(overlap_venn_sen),
    size = 6
  ) +
  scale_fill_gradient(
    low = "#F4FAFE",
    high = "#4981BF"
  ) +
  scale_color_manual(
    values = c("CellChat" = "grey", "scDiffCom" = "grey")
  ) +
  theme_void() +
  ggtitle(
    "Detected CCIs"
  ) +
  theme(
    plot.title = element_text(size = 20)
  )
overlap_venn_sen

## Find dominant pathways from sen cells #####

rpeSen_sender_pathways <- t(sapply(
  cellchat_sen@netP$centr,
  function(i) {
    c(i$outdeg[["RPE-sen"]], mean(i$outdeg), sum(i$outdeg > i$outdeg[["RPE-sen"]]))
  }
))

rpeNotSen_sender_pathways <- t(sapply(
  cellchat_sen@netP$centr,
  function(i) {
    c(i$outdeg[["RPE"]], mean(i$outdeg), sum(i$outdeg > i$outdeg[["RPE"]]))
  }
))

rpeSen_receiver_pathways <- t(sapply(
  cellchat_sen@netP$centr,
  function(i) {
    c(i$indeg[["RPE-sen"]], mean(i$indeg), sum(i$indeg > i$indeg[["RPE-sen"]]))
  }
))

endoSen_sender_pathways <- t(sapply(
  cellchat_sen@netP$centr,
  function(i) {
    c(i$outdeg[["endothelial-sen"]], mean(i$outdeg), sum(i$outdeg > i$outdeg[["endothelial-sen"]]))
  }
))

endoNoSen_sender_pathways <- t(sapply(
  cellchat_sen@netP$centr,
  function(i) {
    c(i$outdeg[["endothelial"]], mean(i$outdeg), sum(i$outdeg > i$outdeg[["endothelial"]]))
  }
))

endoSen_receiver_pathways <- t(sapply(
  cellchat_sen@netP$centr,
  function(i) {
    c(i$indeg[["endothelial-sen"]], mean(i$indeg), sum(i$indeg > i$indeg[["endothelial-sen"]]))
  }
))

# rpe_mediator_pathways <- t(sapply(
#   cellchat_healthy@netP$centr,
#   function(i) {
#     c(i$flowbet[[5]], mean(i$flowbet), sum(i$flowbet > i$flowbet[[5]]))
#   }
# ))
# 
# rpe_influencer_pathways <- t(sapply(
#   cellchat_healthy@netP$centr,
#   function(i) {
#     c(i$info[[5]], mean(i$info), sum(i$info > i$info[[5]]))
#   }
# ))

## Cellchat Senescence network and centrality ####

netVisual_aggregate(
  cellchat_sen,
  signaling = "VEGF"
)
netVisual_heatmap(
  cellchat_sen,
  signaling = "VEGF",
  color.heatmap = "Reds",
  title.name = "VEGF signaling",
  width = 10,
  height = 10,
  font.size = 12
)
netAnalysis_signalingRole_network(
  cellchat_sen,
  signaling = c("VEGF"),
  width = 10,
  height = 3,
  font.size = 12
)

netVisual_aggregate(
  cellchat_sen,
  signaling = "BMP"
)
netVisual_heatmap(
  cellchat_sen,
  signaling = "BMP",
  color.heatmap = "Reds",
  title.name = "BMP signaling",
  width = 10,
  height = 10,
  font.size = 12
)
netAnalysis_signalingRole_network(
  cellchat_sen,
  signaling = c("BMP"),
  width = 10,
  height = 3,
  font.size = 12
)

netVisual_aggregate(
  cellchat_sen,
  signaling = "TENASCIN"
)
netVisual_heatmap(
  cellchat_sen,
  signaling = "TENASCIN",
  color.heatmap = "Reds",
  title.name = "TENASCIN signaling",
  width = 10,
  height = 10,
  font.size = 12
)
netAnalysis_signalingRole_network(
  cellchat_sen,
  signaling = c("TENASCIN"),
  width = 10,
  height = 3,
  font.size = 12
)

## Prepare supplementary tables C6-7 for CCIs #####

cci_cellchat_sen[
  ,
  in_scDiffCom := CCI %in% scd_sen@cci_table_detected$CCI
]
table(cci_cellchat_sen$in_scDiffCom)
sup_table_c6 <- cci_cellchat_sen[, c(1:11, 20)][order(-prob)]
fwrite(
  sup_table_c6,
  "../results/C3_sup_table_c6_sen.csv"
)

cci_scd_sen <- scd_sen@cci_table_detected

cci_scd_sen[, in_CellChat := ifelse(
  CCI %in% cci_cellchat_sen$CCI, TRUE, FALSE
)]
table(cci_scd_sen$in_CellChat)
cci_scd_sen[
  LRI_human$LRI_curated,
  on = "LRI",
  LRI_origin := i.DATABASE
]
sup_table_c7 <- cci_scd_sen[, c(3:5, 23:25, 31, 30)][order(-CCI_SCORE)]
fwrite(
  sup_table_c7,
  "../results/C3_sup_table_c7_sen.csv"
)


## scDiffCom senescence supplementary #####

scd_heatmap_sen <- function(cci_dt, title.name) {
  scd_template <- CJ(
    factor(levels(cellchat_sen@idents), levels(cellchat_sen@idents)),
    factor(levels(cellchat_sen@idents), levels(cellchat_sen@idents)),
    sorted = FALSE
  )
  setnames(
    scd_template,
    old = colnames(scd_template),
    new = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE")
  )
  dt <- copy(scd_template)[
    cci_dt[
      ,
      sum(CCI_SCORE),
      by = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE")
    ],
    on = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE"),
    score_sum := i.V1
  ]
  dt[is.na(dt)] <- 0
  mat <- reshape2::acast(
    dt,
    EMITTER_CELLTYPE ~ RECEIVER_CELLTYPE,
    value.var = "score_sum"
  )
  color.use <- scPalette(ncol(mat))
  names(color.use) <- colnames(mat)
  color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = "Reds")))(100)
  
  df<- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),which = "column",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), which = "row",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
  mat[mat == 0] <- NA
  Heatmap(
    mat,
    col = color.heatmap.use,
    na_col = "white",
    name = "CCI Score Sum",
    bottom_annotation = col_annotation,
    left_annotation = row_annotation,
    top_annotation = ha2,
    right_annotation = ha1,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    row_names_rot = 0,
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    column_title = title.name,
    column_title_gp = gpar(fontsize = 10),
    column_names_rot = 90,
    row_title = "Sources (Sender)",
    row_title_gp = gpar(fontsize = 10),
    row_title_rot = 90,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 8, fontface = "plain"),
      title_position = "leftcenter-rot",
      border = NA, #at = colorbar.break,
      legend_height = unit(20, "mm"),
      labels_gp = gpar(fontsize = 8),
      grid_width = unit(2, "mm")
    )
  )
}

cci_cellchat_sen[LRI %in% vegf_pathway$LRI]
cci_cellchat_sen[pathway_name == "VEGF"]

ggVennDiagram(
  x = list(
    cci_scd_sen[LRI %in% bmp_pathway$LRI]$CCI,
    cci_cellchat_sen[pathway_name == "BMP"]$CCI
  )
)

sort(unique(cci_cellchat_sen[LRI %in% vegf_pathway$LRI]$LRI))
sort(unique(cci_scd_sen[LRI %in% vegf_pathway$LRI]$LRI))

cci_scd_sen[LRI %in% sort(unique(cci_cellchat_sen[LRI %in% vegf_pathway$LRI]$LRI))]


test <- cci_cellchat_sen[pathway_name == "BMP"]
setnames(
  test,
  old = c("source", "target", "prob"),
  new = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "CCI_SCORE")
)



scd_heatmap_sen(cci_scd_sen[LRI %in% vegf_pathway$LRI], "VEGF signaling network")


scd_heatmap_sen(cci_scd_sen[LRI %in% vegf_pathway2$LRI], "VEGF signaling network")
scd_heatmap_sen(cci_scd_sen[LRI %in% sort(unique(cci_cellchat_sen[LRI %in% vegf_pathway$LRI]$LRI))], "VEGF signaling network")

scd_heatmap_sen(test, "signaling network")


scd_heatmap_sen(cci_scd_sen[LRI %in% bmp_pathway$LRI], "BMP signaling network")

scd_heatmap_sen(cci_scd_sen[LRI %in% tenascin_pathway$LRI], "TENASCIN signaling network")

###########

netVisual_circle(
  cellchat_sen@net$count,
  vertex.weight = groupSize, weight.scale = T,
  label.edge= F, title.name = "Number of interactions"
)
dev.off()

netVisual_circle(cellchat_sen@net$weight, 
                 vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")

netVisual_bubble(cellchat_sen,
                 sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)


plotGeneExpression(cellchat_sen, signaling = "VEGF")


pathways1.show <- c("VEGF") 
vertex.receiver = seq(1,4) # a numeric vector. 
#,  vertex.receiver = vertex.receiver)
dev.off()


## Explore senescence communication ####

cci_sen_all <- copy(scd_sen@cci_table_detected)

cci_sen_vegf <- cci_sen_all[LRI %in% vegf_pathway$LRI]
cci_sen_bmp <- cci_sen_all[LRI %in% bmp_pathway$LRI]
cci_sen_tenascin <- cci_sen_all[LRI %in% tenascin_pathway$LRI]

cci_ER_tenascin <- cci_sen_tenascin[
  ,
  sum(CCI_SCORE),
  by = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE")
][order(-V1)]
setnames(cci_ER_tenascin, old = "V1", new = "weight")
cci_ER_tenascin
cci_ER_tenascin_2 <- dcast.data.table(
  cci_ER_tenascin,
  formula = EMITTER_CELLTYPE ~ RECEIVER_CELLTYPE,
  value.var = "weight",
  fill = 0
)
scd_tenascin_mat <- as.matrix(cci_ER_tenascin_2[,-1])
rownames(scd_tenascin_mat) <- cci_ER_tenascin_2[, 1]$EMITTER_CELLTYPE
heatmap(scd_tenascin_mat, scale = "none")

par(mfrow=c(1,1))
netVisual_heatmap(
  cellchat_sen,
  signaling = pathways_tenascin.show,
  color.heatmap = "Reds",
  font.size = 10
)
netAnalysis_signalingRole_network(
  cellchat_sen,
  signaling = pathways_tenascin.show,
  width = 8,
  height = 2.5,
  font.size = 10
)

cci_ER_vegf <- cci_sen_vegf[
  ,
  sum(CCI_SCORE),
  by = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE")
][order(-V1)]
setnames(cci_ER_vegf, old = "V1", new = "weight")
cci_ER_vegf
cci_ER_vegf_2 <- dcast.data.table(
  cci_ER_vegf,
  formula = EMITTER_CELLTYPE ~ RECEIVER_CELLTYPE,
  value.var = "weight",
  fill = 0
)
scd_vegf_mat <- as.matrix(cci_ER_vegf_2[,-1])
rownames(scd_vegf_mat) <- cci_ER_vegf_2[, 1]$EMITTER_CELLTYPE
heatmap(scd_vegf_mat, scale = "none")

par(mfrow=c(1,1))
netVisual_heatmap(
  cellchat_sen,
  signaling = pathways_vegf.show,
  color.heatmap = "Reds",
  font.size = 10
)
netAnalysis_signalingRole_network(
  cellchat_sen,
  signaling = pathways_vegf.show,
  width = 8,
  height = 2.5,
  font.size = 10
)

cci_ER_bmp <- cci_sen_bmp[
  ,
  sum(CCI_SCORE),
  by = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE")
][order(-V1)]
setnames(cci_ER_bmp, old = "V1", new = "weight")
cci_ER_bmp
cci_ER_bmp_2 <- dcast.data.table(
  cci_ER_bmp,
  formula = EMITTER_CELLTYPE ~ RECEIVER_CELLTYPE,
  value.var = "weight",
  fill = 0
)
scd_bmp_mat <- as.matrix(cci_ER_bmp_2[,-1])
rownames(scd_bmp_mat) <- cci_ER_bmp_2[, 1]$EMITTER_CELLTYPE
heatmap(scd_bmp_mat, scale = "none")

par(mfrow=c(1,1))
netVisual_heatmap(
  cellchat_sen,
  signaling = pathways_bmp.show,
  color.heatmap = "Reds",
  font.size = 10
)
netAnalysis_signalingRole_network(
  cellchat_sen,
  signaling = pathways_bmp.show,
  width = 8,
  height = 2.5,
  font.size = 10
)

## over all pathways ####

cldb_pathways <- lapply(
  unique(cldb$pathway_name),
  function(i) {
    LRI_human$LRI_curated[
      LIGAND_1 %in% cldb[pathway_name == i]$ligand |
        LIGAND_2 %in% cldb[pathway_name == i]$ligand |
        RECEPTOR_1 %in% cldb[pathway_name == i]$receptor |
        RECEPTOR_2 %in% cldb[pathway_name == i]$receptor |
        RECEPTOR_3 %in% cldb[pathway_name == i]$receptor 
    ]
  }
)
names(cldb_pathways) <- unique(cldb$pathway_name)

cci_ER_allpathways <- lapply(
  cldb_pathways,
  function(i) {
    temp <- cci_sen_all[LRI %in% i$LRI]
    temp[
      ,
      sum(CCI_SCORE),
      by = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE")
    ][order(-V1)]
  }
)

cci_ER_allpathways <- rbindlist(
  cci_ER_allpathways,
  idcol = "pathway"
)

cci_ER_tenascin <- cci_sen_tenascin[
  ,
  sum(CCI_SCORE),
  by = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE")
][order(-V1)]


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


test <- icc_sen_scdiffcom$icc_diff_healthy_vs_amd@cci_table_detected
test <- GetTableCCI(icc_sen_scdiffcom$icc_diff_healthy_vs_amd, "detected", TRUE)

PlotORA(
  icc_sen_scdiffcom$icc_diff_healthy_vs_amd,
  category = "LRI"
)
BuildNetwork(
  icc_sen_scdiffcom$icc_diff_healthy_vs_amd
)

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

## Global pathways ####

netAnalysis_signalingRole_heatmap(
  cellchat,
  width = 20,height = 20
)

rankNet(cellchat_healthy, sources.use = "RPE", targets.use = "endothelial", mode = "single", do.stat = TRUE)

colMeans(cellchat@netP[[2]]["RPE", , ]) + colSums(cellchat@netP[[2]][,"RPE" , ])

test <- colMeans(cellchat@netP[[2]]["RPE", , ]) > 
  sapply(1:dim(cellchat@netP[[2]])[[3]], function(i) {mean(cellchat@netP[[2]][, , i])})

test[test]

sort(colSums(cellchat@netP[[2]]["RPE", , ]), decreasing = TRUE)[1:15]
sort(colMeans(cellchat@netP[[2]]["RPE", , ]), decreasing = TRUE)[1:15]

netVisual_heatmap(
  cellchat,
  signaling = "COLLAGEN",
  color.heatmap = "Reds",
  font.size = 10
)

netAnalysis_signalingRole_heatmap(
  cellchat,
  signaling = c("VEGF","BMP","BTLA","GDF","LIGHT","GRN","HGF"),
  pattern = "all",width = 20,height = 16,font.size = 7.5,font.size.title = 10,cluster.rows = FALSE,cluster.cols = FALSE)

rankNet(cellchatmerge, sources.use="RPE",targets.use="endothelial",mode = "comparison", stacked = T, do.stat = T,cutoff.pvalue=0.05)

compareInteractions(cellchatmerge, show.legend = F, group = c(1,2))
rankNet(cellchatmerge, mode = "comparison", stacked = T, do.stat = TRUE)
rankNet(cellchatmerge, sources.use="RPE",targets.use="endothelial",mode = "comparison", stacked = T, do.stat = T,cutoff.pvalue=0.05)

pathways_vegf.show <- c("VEGF")
pathways_bmp.show <- c("BMP")
pathways_tenascin.show <- c("TENASCIN")

## old stuff ######


cellchatmerge <- netAnalysis_computeCentrality(cellchatmerge)

netAnalysis_signalingRole_heatmap(
  cellchatmerge,
  signaling = c("VEGF","BMP","BTLA","GDF","LIGHT","GRN","HGF"),
  pattern = "all",width = 20,height = 16,font.size = 7.5,font.size.title = 10,cluster.rows = FALSE,cluster.cols = FALSE)


netAnalysis_signalingRole_heatmap(cellchat_amd, signaling = c("VEGF","BMP","BTLA","GDF","LIGHT","GRN","HGF"),pattern = "all")
#,width = 20,height = 16,font.size = 7.5,font.size.title = 10,cluster.rows = TRUE,cluster.cols = TRUE)
i <- 1
netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = c("VEGF","BMP","APP","GDF","LIGHT","GRN","HGF"), title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = c("VEGF","BMP","BTLA","GDF","LIGHT","GRN","HGF"), title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))



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