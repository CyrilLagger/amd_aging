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

cellchat <- createCellChat(
    object = amd_seurat,
    group.by = "cell_type"
)
levels(cellchat@idents)
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 5)

df_net <- subsetCommunication(cellchat)

## Senescence CellChat workflow ####

cellchat_sen <- createCellChat(
    object = amd_seurat,
    group.by = "cell_type_senescence"
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

df_sen_net <- subsetCommunication(cellchat_sen)

## Save senescence results ####

saveRDS(
    cellchat_sen,
    paste0(
        path_results,
        "C2_cellchat_sen.rds"
  )
)