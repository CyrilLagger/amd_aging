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
      #"Kasit_DOWN" = senref_kasit_down$gene,
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
amd_seurat$score_seurat_mayo1
amd_seurat$score_seurat_kasit1
amd_seurat$score_seurat_reactome1

## Score with enrichIt ####

sen_enrichit <- enrichIt(
  obj = amd_seurat,
  gene.sets = list(sen = senref_mayo$gene),
  groups = 1000,
  cores = 20
)
sen_nes
amd_seurat$score_gsea_mayo <- sen_nes
amd_seurat$score_gsea_mayo_top20pct <- ifelse(
  amd_seurat$score_gsea_mayo >= quantile(amd_seurat$score_gsea_mayo, 0.5),
  TRUE,
  FALSE
)



cor(
  amd_seurat$score_seurat_mayo1,
  amd_seurat$score_seurat_kasit1
)
cor(
  amd_seurat$score_seurat_mayo1,
  amd_seurat$score_seurat_reactome1
)

cor(
  amd_seurat$score_seurat_reactome1,
  amd_seurat$score_seurat_kasit1
)

amd_seurat$is_rpe <- ifelse(
  amd_seurat$cell_type == "RPE", TRUE, FALSE
)
table(amd_seurat$is_rpe)

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
    x = score_seurat_kasit1,
    color = cell_type,
    fill = cell_type
  )
) +
  geom_histogram(
    bins = 50,
    position = "identity",
    alpha = 0.3
  )

p <- ggplot(
  amd_seurat[[]],
  aes(
    x = score_seurat_kasit1,
    y = score_seurat_mayo1,
    color = is_rpe
  )
) +
  geom_point()
p <- ggMarginal(
  p,
  groupColour = TRUE,
  groupFill = TRUE
)
ggsave(
  "../results/score_mayo_vs_kasit_rpe.png",
  p
)


table(
  amd_seurat$cell_type,
  amd_seurat$score_seurat_mayo1 >=
    quantile(amd_seurat$score_seurat_mayo1, 0.5)
)
table(
  amd_seurat$cell_type,
  amd_seurat$score_seurat_kasit1 >=
    quantile(amd_seurat$score_seurat_kasit1, 0.8)
)

amd_seurat$score_seurat_kasit_20pct <- ifelse(
  amd_seurat$score_seurat_kasit1 >= quantile(amd_seurat$score_seurat_kasit1, 0.8),
  TRUE,
  FALSE
)

amd_seurat$score_seurat_mayo_top50pct <- ifelse(
  amd_seurat$score_seurat_mayo1 >= quantile(amd_seurat$score_seurat_mayo1, 0.5),
  TRUE,
  FALSE
)


amd_seurat$score_seurat_kasit1
amd_seurat$score_seurat_kasit_top20pct <- ifelse(
  amd_seurat$score_seurat_kasit1 >= quantile(amd_seurat$score_seurat_kasit1, 0.5),
  TRUE,
  FALSE
)

amd_seurat <- AddModuleScore(
  amd_seurat,
  features = list(sen = senref_kasit_down$gene),
  name = "score_seurat_senref_kasit_down"
)


amd_seurat$score_seurat_reactome1
amd_seurat$score_seurat_reactome_top20pct <- ifelse(
  amd_seurat$score_seurat_reactome1 >= quantile(amd_seurat$score_seurat_reactome1, 0.5),
  TRUE,
  FALSE
)

table(
  amd_seurat$score_seurat_mayo_top20pct,
  amd_seurat$score_gsea_mayo_top20pct
)

table(
  amd_seurat$score_seurat_mayo_top20pct,
  amd_seurat$score_seurat_reactome_top20pct
)

ftable(
  amd_seurat$cell_type,
  amd_seurat$score_seurat_mayo_top20pct,
  amd_seurat$score_seurat_reactome_top20pct,
  amd_seurat$score_seurat_kasit_top20pct
)

ftable(
  amd_seurat$cell_type,
  amd_seurat$score_seurat_mayo_top20pct,
  amd_seurat$score_seurat_kasit_top20pct
)

ftable(
  amd_seurat$cell_type,
  amd_seurat$score_seurat_kasit_top20pct
)

cor(
  amd_seurat$score_gsea_mayo,
  amd_seurat$score_seurat_mayo1,
  method = "spearman"
)
cor(
  amd_seurat$score_seurat_kasit1,
  amd_seurat$score_seurat_mayo1,
  method = "spearman"
)
cor(
  amd_seurat$score_seurat_kasit1,
  amd_seurat$score_gsea_mayo,
  method = "spearman"
)
cor(
  amd_seurat$score_seurat_mayo1,
  amd_seurat$score_seurat_reactome1,
  method = "spearman"
)
cor(
  amd_seurat$score_seurat_kasit1,
  amd_seurat$score_seurat_reactome1,
  method = "spearman"
)


ggplot(
  amd_seurat[[]],
  aes(
    x = score_gsea_mayo,
    fill = cell_type,
    color = cell_type
  )
) +
  geom_histogram(
    aes(y = ..density..),
    position = "identity",
    alpha = 0.4,
    bins = 50
  )




ggplot(
  amd_seurat[[]],
  aes(
    x = score_gsea_mayo,
    y = score_seurat_mayo1
  )
) +
  geom_point()

ggplot(
  amd_seurat[[]],
  aes(
    x = score_seurat_kasit1,
    y = score_seurat_senref_kasit_down1
  )
) +
  geom_point()


ggplot(
  amd_seurat[[]],
  aes(
    x = score_seurat_reactome1,
    y = score_seurat_mayo1
  )
) +
  geom_point()


ggplot(
  amd_seurat[[]],
  aes(
    x = score_seurat_reactome1,
    y = score_seurat_kasit1,
    color = is_rpe
  )
) +
  geom_point()


amd_seurat$is_rpe <- ifelse(
  amd_seurat$cell_type == "RPE", TRUE, FALSE
)
table(amd_seurat$is_rpe)

ggplot(
  amd_seurat[[]],
  aes(
    x = score_seurat_mayo1,
    y = score_seurat_kasit1
  )
) +
  geom_point()

ggplot(
  amd_seurat[[]],
  aes(
    x = score_gsea_mayo,
    fill = age,
    color = age
  )
) +
  geom_histogram(
    aes(y = ..density..),
    position = "identity",
    alpha = 0.4,
    bins = 50
  )