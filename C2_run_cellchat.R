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

plan(multicore, workers = 24)

## Prepare CellChat workflow ####

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB

## Cellchat on standard cells ####

cellchat_all <- createCellChat(
    object = amd_seurat,
    group.by = "cell_type"
)
levels(cellchat_all@idents)
cellchat_all@DB <- CellChatDB.use
cellchat_all <- subsetData(cellchat_all)
cellchat_all <- identifyOverExpressedGenes(cellchat_all)
cellchat_all <- identifyOverExpressedInteractions(cellchat_all)
cellchat_all <- projectData(cellchat_all, PPI.human)
cellchat_all <- computeCommunProb(cellchat_all, type = "truncatedMean", trim = 0.1)
cellchat_all <- filterCommunication(cellchat_all, min.cells = 10)
cellchat_all <- computeCommunProbPathway(cellchat_all)
cellchat_all <- aggregateNet(cellchat_all)
cellchat_all <- netAnalysis_computeCentrality(cellchat_all)
cci_cellchat_all <- subsetCommunication(cellchat_all)

cellchat_healthy <- createCellChat(
  object = subset(
    amd_seurat,
    subset = age %in% c("54 year", "82 year")
  ),
  group.by = "cell_type"
)
levels(cellchat_healthy@idents)
cellchat_healthy@DB <- CellChatDB.use
cellchat_healthy <- subsetData(cellchat_healthy)
cellchat_healthy <- identifyOverExpressedGenes(cellchat_healthy)
cellchat_healthy <- identifyOverExpressedInteractions(cellchat_healthy)
cellchat_healthy <- projectData(cellchat_healthy, PPI.human)
cellchat_healthy <- computeCommunProb(cellchat_healthy, type = "truncatedMean", trim = 0.1)
cellchat_healthy <- filterCommunication(cellchat_healthy, min.cells = 10)
cellchat_healthy <- computeCommunProbPathway(cellchat_healthy)
cellchat_healthy <- aggregateNet(cellchat_healthy)
cellchat_healthy <- netAnalysis_computeCentrality(cellchat_healthy)
cci_cellchat_healthy <- subsetCommunication(cellchat_healthy)


cellchat_amd <- createCellChat(
  object = subset(amd_seurat, subset = age == "79 year"),
  group.by = "cell_type"
)
levels(cellchat_amd@idents)
cellchat_amd@DB <- CellChatDB.use
cellchat_amd <- subsetData(cellchat_amd)
cellchat_amd <- identifyOverExpressedGenes(cellchat_amd)
cellchat_amd <- identifyOverExpressedInteractions(cellchat_amd)
cellchat_amd <- projectData(cellchat_amd, PPI.human)
cellchat_amd <- computeCommunProb(cellchat_amd, type = "truncatedMean", trim = 0.1)
cellchat_amd <- filterCommunication(cellchat_amd, min.cells = 10)
cellchat_amd <- computeCommunProbPathway(cellchat_amd)
cellchat_amd <- aggregateNet(cellchat_amd)
cellchat_amd <- netAnalysis_computeCentrality(cellchat_amd)
cci_cellchat_amd <- subsetCommunication(cellchat_amd)

## Senescence CellChat workflow ####

cellchat_sen <- createCellChat(
    object = amd_seurat,
    group.by = "cell_type_senescence_gsea"
)
levels(cellchat_sen@idents)
cellchat_sen@DB <- CellChatDB.use
cellchat_sen <- subsetData(cellchat_sen)
cellchat_sen <- identifyOverExpressedGenes(cellchat_sen)
cellchat_sen <- identifyOverExpressedInteractions(cellchat_sen)
cellchat_sen <- projectData(cellchat_sen, PPI.human)
cellchat_sen <- computeCommunProb(cellchat_sen, type = "truncatedMean", trim = 0.1)
cellchat_sen <- filterCommunication(cellchat_sen, min.cells = 10)
cellchat_sen <- computeCommunProbPathway(cellchat_sen)
cellchat_sen <- aggregateNet(cellchat_sen)
cellchat_sen <- netAnalysis_computeCentrality(cellchat_sen)
cci_cellchat_sen <- subsetCommunication(cellchat_sen)

#saveRDS(cellchat_sen, "../../../../../cellchat_sen.rds")

## Save senescence results ####

saveRDS(
    cellchat_sen,
    paste0(
        path_results,
        "C2_cellchat_sen.rds"
  )
)