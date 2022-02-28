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

plan(multicore, workers = 20)

## Load senescence reference gene lists ####

senref_mayo <- fread(
  paste0(
    path_data,
    "B1_senref_mayo.csv"
  ),
  header = FALSE
)
colnames(senref_mayo) <- "gene"
senref_mayo$gene[!senref_mayo$gene %in% rownames(amd_seurat)]

senref_kasit_up <- fread(
  paste0(
    path_data,
    "B1_senref_kasit_up.txt"
  ),
  header = FALSE
)
colnames(senref_kasit_up) <- "gene"
senref_kasit_up$gene[!senref_kasit_up$gene %in% rownames(amd_seurat)]

senref_kasit_down <- fread(
  paste0(
    path_data,
    "B1_senref_kasit_down.txt"
  ),
  header = FALSE
)
colnames(senref_kasit_down) <- "gene"
senref_kasit_down$gene[!senref_kasit_down$gene %in% rownames(amd_seurat)]

senref_reactome <- getGeneSets(
  library = "C2",
  gene.sets = "REACTOME_CELLULAR_SENESCENCE"
)
senref_reactome <- GSEABase::geneIds(
  senref_reactome
)$REACTOME_CELLULAR_SENESCENCE
senref_reactome <- data.table(
  gene = senref_reactome
)

## Overlap senescence gene lists #####

senref_overlap <- ggVennDiagram::process_data(
  ggVennDiagram::Venn(
    list(
      "MayoClinic" = senref_mayo$gene,
      "Kasit_UP" = senref_kasit_up$gene,
      # "Kasit_DOWN" = senref_kasit_down$gene,
      "REACTOME" = senref_reactome$gene,
      "Seurat" = rownames(amd_seurat)
    )
  )
)

senref_overlap_plot <- ggplot() +
  geom_sf(
    aes(fill = count),
    data = venn_region(senref_overlap),
    show.legend = FALSE
  ) +
  geom_sf(
    aes(color = name),
    data = venn_setedge(senref_overlap),
    show.legend = FALSE,
    size = 2
  ) +
  geom_sf_text(
    aes(label = name),
    data = venn_setlabel(senref_overlap),
    size = 5
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
    data = venn_region(senref_overlap),
    size = 5
  ) +
  scale_fill_gradient(
    low = "#F4FAFE",
    high = "#4981BF"
  ) +
  theme_void() +
  ggtitle(
    "Overlap Senescence Reference"
  ) +
  theme(
    plot.title = element_text(size = 20)
  )
senref_overlap_plot

ggsave(
  paste0(
    path_results,
    "images/B1_senref_overlap.png"
  ),
  senref_overlap_plot,
  scale = 2
)

## Create a combined senescence gene list ####

senref_combined <- unique(
  rbindlist(
    l = list(
      senref_mayo,
      senref_kasit_up,
      senref_reactome
    )
  )
)

## Score with Seurat ####

amd_seurat <- AddModuleScore(
  amd_seurat,
  features = list(sen = senref_mayo$gene),
  name = "score_seurat_mayo"
)
amd_seurat <- AddModuleScore(
  amd_seurat,
  features = list(sen = senref_kasit_up$gene),
  name = "score_seurat_kasit_up"
)
amd_seurat <- AddModuleScore(
  amd_seurat,
  features = list(sen = senref_kasit_down$gene),
  name = "score_seurat_kasit_down"
)
amd_seurat <- AddModuleScore(
  amd_seurat,
  features = list(sen = senref_reactome$gene),
  name = "score_seurat_reactome"
)
amd_seurat <- AddModuleScore(
  amd_seurat,
  features = list(sen = senref_combined$gene),
  name = "score_seurat_combined"
)
amd_seurat$score_seurat_mayo1
amd_seurat$score_seurat_kasit_up1
amd_seurat$score_seurat_kasit_down1
amd_seurat$score_seurat_reactome1
amd_seurat$score_seurat_combined1

amd_seurat$score_seurat_kasit_updown <- amd_seurat$score_seurat_kasit_up1 -
  amd_seurat$score_seurat_kasit_down1

## Score with enrichIt ####

sen_enrichit <- enrichIt(
  obj = amd_seurat,
  gene.sets = list(
    score_it_mayo = senref_mayo$gene,
    score_it_kasit_up = senref_kasit_up$gene,
    score_it_kasit_down = senref_kasit_down$gene,
    score_it_reactome = senref_reactome$gene,
    score_it_combined = senref_combined$gene
  ),
  groups = 1000,
  cores = 20
)
# TODO
# amd_seurat$score_gsea_mayo <- sen_nes

## Correlations between scores ####

cor(
  amd_seurat$score_seurat_mayo1,
  amd_seurat$score_seurat_kasit_up1
)
cor(
  amd_seurat$score_seurat_mayo1,
  amd_seurat$score_seurat_reactome1
)
cor(
  amd_seurat$score_seurat_reactome1,
  amd_seurat$score_seurat_kasit_up1
)

## Plot scores ####

ggplot(
  amd_seurat[[]],
  aes(
    x = score_seurat_mayo1,
    color = cell_type,
    fill = cell_type
  )
) +
  geom_histogram(
    bins = 50,
    position = "identity",
    alpha = 0.3
  )

ggplot(
  amd_seurat[[]],
  aes(
    x = score_seurat_kasit_up1,
    color = cell_type,
    fill = cell_type
  )
) +
  geom_histogram(
    bins = 50,
    position = "identity",
    alpha = 0.3
  )

ggplot(
  amd_seurat[[]],
  aes(
    x = score_seurat_combined1,
    color = cell_type,
    fill = cell_type
  )
) +
  geom_histogram(
    bins = 50,
    position = "identity",
    alpha = 0.3
  )

sen_score_km_plot <- ggplot(
  amd_seurat[[]],
  aes(
    x = score_seurat_kasit_up1,
    y = score_seurat_mayo1,
    color = is_rpe
  )
) +
  geom_point()
sen_score_km_plot <- ggMarginal(
  sen_score_km_plot,
  groupColour = TRUE,
  groupFill = TRUE
)
ggsave(
  paste0(
    path_results,
    "images/B1_sen_score_km.png"
  ),
  sen_score_km_plot
)

sen_score_kc_plot <- ggplot(
  amd_seurat[[]],
  aes(
    x = score_seurat_kasit_up1,
    y = score_seurat_combined1,
    color = is_rpe
  )
) +
  geom_point()
sen_score_kc_plot <- ggMarginal(
  sen_score_kc_plot,
  groupColour = TRUE,
  groupFill = TRUE
)
ggsave(
  paste0(
    path_results,
    "images/B1_sen_score_kc.png"
  ),
  sen_score_kc_plot
)

## Find appropriate cutoff ####

table(
  amd_seurat$cell_type,
  amd_seurat$score_seurat_mayo1 >=
    quantile(amd_seurat$score_seurat_mayo1, 0.8)
)
table(
  amd_seurat$cell_type,
  amd_seurat$score_seurat_kasit_up1 >=
    quantile(amd_seurat$score_seurat_kasit_up1, 0.8)
)

table(
  amd_seurat$cell_type,
  amd_seurat$score_seurat_kasit_updown >=
    quantile(amd_seurat$score_seurat_kasit_updown, 0.8)
)

table(
  amd_seurat$cell_type,
  amd_seurat$score_seurat_combined1 >=
    quantile(amd_seurat$score_seurat_combined1, 0.8)
)

amd_seurat$sen_kasit_up_20pct <- ifelse(
  amd_seurat$score_seurat_kasit_up1 >=
    quantile(amd_seurat$score_seurat_kasit_up1, 0.8),
  TRUE,
  FALSE
)
table(amd_seurat$sen_kasit_up_20pct)

amd_seurat$sen_kasit_updown_20pct <- ifelse(
  amd_seurat$score_seurat_kasit_updown >=
    quantile(amd_seurat$score_seurat_kasit_updown, 0.8),
  TRUE,
  FALSE
)
table(amd_seurat$sen_kasit_updown_20pct)
ftable(
  amd_seurat$age,
  amd_seurat$sen_kasit_updown_20pct,
  amd_seurat$cell_type
)

## Rename cell types based on senesence scores ####

amd_seurat$cell_type_senescence <- ifelse(
  amd_seurat$sen_kasit_updown_20pct == TRUE,
  paste0(amd_seurat$cell_type, "_sen"),
  paste0(amd_seurat$cell_type, "")
)
table(amd_seurat$cell_type_senescence)

## Plots for article ####

sen_score_k_updown_plot <- ggplot(
  amd_seurat[[]],
  aes(
    x = score_seurat_kasit_up1,
    y = score_seurat_kasit_down1
  )
) + geom_point(
  ) + xlab(
    "Senescence score (up-regulated signatures)"
  ) + ylab(
    "Senescence score (down-regulated signatures)"
  )
sen_score_k_updown_plot <- ggMarginal(
  sen_score_k_updown_plot
)
ggsave(
  paste0(
    path_results,
    "images/B1_sen_score_k_updown.png"
  ),
  sen_score_k_updown_plot
)

sen_score_k_updown_plot_ct <- ggplot(
  amd_seurat[[]],
  aes(
    x = score_seurat_kasit_up1,
    y = score_seurat_kasit_down1,
    color = cell_type
  )
) +
  geom_point(
  ) + xlab(
    "Senescence score (up-regulated signatures)"
  ) + ylab(
    "Senescence score (down-regulated signatures)"
  ) + labs(
    color = "Cell type"
  )
sen_score_k_updown_plot_ct <- ggMarginal(
  sen_score_k_updown_plot_ct,
  groupColour = TRUE,
  groupFill = TRUE
)
ggsave(
  paste0(
    path_results,
    "images/B1_sen_score_k_updown_ct.png"
  ),
  sen_score_k_updown_plot_ct
)

sen_score_k_updown_plot_rpe <- ggplot(
  amd_seurat[[]],
  aes(
    x = score_seurat_kasit_up1,
    y = score_seurat_kasit_down1,
    color = is_rpe
  )
) +
  geom_point(
  ) + xlab(
    "Senescence score (up-regulated signatures)"
  ) + ylab(
    "Senescence score (down-regulated signatures)"
  ) + labs(
    color = "RPE"
  )
sen_score_k_updown_plot_rpe <- ggMarginal(
  sen_score_k_updown_plot_rpe,
  groupColour = TRUE,
  groupFill = TRUE
)
ggsave(
  paste0(
    path_results,
    "images/B1_sen_score_k_updown_rpe.png"
  ),
  sen_score_k_updown_plot_rpe
)

sen_score_k_updown_hist <- ggplot(
  amd_seurat[[]],
  aes(
    x = score_seurat_kasit_updown
  )
) + geom_histogram(
  bins = 50,
  alpha = 0.5
) + geom_vline(
  xintercept = quantile(amd_seurat$score_seurat_kasit_updown, 0.8)
) + xlab(
  "Senescence score (up - down signatures)"
)
ggsave(
  paste0(
    path_results,
    "images/B1_sen_score_k_updown_hist.png"
  ),
  sen_score_k_updown_hist
)

sen_score_k_updown_hist_ct <- ggplot(
  amd_seurat[[]],
  aes(
    x = score_seurat_kasit_updown,
    color = cell_type,
    fill = cell_type
  )
) + geom_histogram(
  bins = 50,
  alpha = 0.5, 
  position = "identity"
) + geom_vline(
  xintercept = quantile(amd_seurat$score_seurat_kasit_updown, 0.8)
) + xlab(
  "Senescence score (up - down signatures)"
)
ggsave(
  paste0(
    path_results,
    "images/B1_sen_score_k_updown_hist_ct.png"
  ),
  sen_score_k_updown_hist_ct
)

sen_score_k_updown_feature <- FeaturePlot(
  amd_seurat,
  features = "score_seurat_kasit_updown",
  cols = c("yellow", "blue")
)
ggsave(
  paste0(
    path_results,
    "images/B1_sen_score_k_updown_feature.png"
  ),
  sen_score_k_updown_feature
)

amd_seurat$senescence_score <- ifelse(
  amd_seurat$sen_kasit_updown_20pct, "High", "Low"
)

sen_cut_k_updown_plot <- DimPlot(
  amd_seurat,
  reduction = "umap",
  group.by = "senescence_score",
  pt.size = 1,
  order = "High"
)
ggsave(
  paste0(
    path_results,
    "images/B1_sen_cut_k_updown.png"
  ),
  sen_cut_k_updown_plot
)

