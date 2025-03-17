fetal <- RenameIdents(fetal, new_idents)

DimPlot(fetal, raster = F, label = T, repel = T, label.box = T, group.by = 'timepoint')

fetal$EK_PB_annov1 <- fetal@active.ident

DimPlot(fetal, raster = F, label = T, repel = T, label.box = T, group.by = 'EK_PB_annov1')

SaveH5Seurat(fetal, 'C://Bioinf/HUMAN_FETAL_RETINA/COMBINED_EKPB_v1.h5Seurat')

fetal_day105 <- subset(fetal, subset = timepoint == c('Day 105'))

fetal_W27 <- subset(fetal, subset = timepoint == c('Week 27'))

AverageExpression(fetal, features = c('BDNF'), group.by = c('EK_PB_annov1', 'timepoint'), assays = 'RNA')

performCellChatAnalysis <- function(cellchat, cellChatDB) {
  # Set the CellChat database
  cellchat@DB <- cellChatDB
  
  # Subset the data
  cellchat <- subsetData(cellchat)
  
  # Identify over-expressed genes
  cellchat <- identifyOverExpressedGenes(cellchat)
  
  # Identify over-expressed interactions
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- projectData(cellchat, PPI.human)
  # Compute communication probability
  cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, raw.use = TRUE, population.size = FALSE)
  
  # Filter communication based on minimum cells
  cellchat <- filterCommunication(cellchat, min.cells = 5)
  
  # Compute communication probability for pathways
  cellchat <- computeCommunProbPathway(cellchat)
  
  # Aggregate networks
  cellchat <- aggregateNet(cellchat)
  
  # Compute centrality for the aggregated network
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  # Return the processed CellChat object
  return(cellchat)
}

# Example usage:
# Assuming 'SV_zebrafish_Adult' is your initial CellChat object and 'CellChatDB.zebrafish' is the database
cc_fetal_d105 <- createCellChat(object = fetal_day105, group.by = "EK_PB_annov1", assay = 'RNA')
CellChatDB.use <- CellChatDB.human

# Perform the analysis
cc_fetal_d105 <- performCellChatAnalysis(cc_fetal_d105, CellChatDB.use)

fetal_day125 <- subset(fetal, subset = timepoint == c('Day 125'))
cc_fetal_d125 <- createCellChat(object = fetal_day125, group.by = "EK_PB_annov1", assay = 'RNA')
CellChatDB.use <- CellChatDB.human

# Perform the analysis
cc_fetal_d125 <- performCellChatAnalysis(cc_fetal_d125, CellChatDB.use)

fetal_day59 <- subset(fetal, subset = timepoint == c('Day 59'))
cc_fetal_d59 <- createCellChat(object = fetal_day59, group.by = "EK_PB_annov1", assay = 'RNA')
CellChatDB.use <- CellChatDB.human

# Perform the analysis
cc_fetal_d59 <- performCellChatAnalysis(cc_fetal_d59, CellChatDB.use)


performCellChatForDay <- function(fetal_data, day, group_by, assay) {
  # Subset fetal data for the specific day
  fetal_subset <- subset(fetal_data, subset = timepoint == day)
  fetal_subset@active.ident <- droplevels(fetal_subset@active.ident, exclude = setdiff(levels(fetal_subset@active.ident), unique(fetal_subset@active.ident)))
  fetal_subset$EK_PB_annov1 <- droplevels(fetal_subset$EK_PB_annov1, exclude = setdiff(levels(fetal_subset$EK_PB_annov1), unique(fetal_subset$EK_PB_annov1)))
  
  
  # Create CellChat object
  cc_fetal <- createCellChat(object = fetal_subset, group.by = group_by, assay = assay)
  
  # Specify the CellChat database to use
  CellChatDB.use <- CellChatDB.human
  
  # Perform CellChat analysis
  cc_fetal <- performCellChatAnalysis(cc_fetal, CellChatDB.use)
  
  return(cc_fetal)
}

cc_fetal_d59 <- performCellChatForDay(fetal, 'Day 59', 'EK_PB_annov1', 'RNA')
cc_fetal_d80 <- performCellChatForDay(fetal, 'Day 80', 'EK_PB_annov1', 'RNA')
cc_fetal_d82 <- performCellChatForDay(fetal, 'Day 82', 'EK_PB_annov1', 'RNA')
cc_fetal_w11 <- performCellChatForDay(fetal, 'Week 11', 'EK_PB_annov1', 'RNA')
cc_fetal_w12 <- performCellChatForDay(fetal, 'Week 12', 'EK_PB_annov1', 'RNA')
cc_fetal_w13 <- performCellChatForDay(fetal, 'Week 13', 'EK_PB_annov1', 'RNA')
cc_fetal_w14 <- performCellChatForDay(fetal, 'Week 14', 'EK_PB_annov1', 'RNA')
cc_fetal_w15 <- performCellChatForDay(fetal, 'Week 15', 'EK_PB_annov1', 'RNA')
cc_fetal_w16 <- performCellChatForDay(fetal, 'Week 16', 'EK_PB_annov1', 'RNA')
cc_fetal_w17 <- performCellChatForDay(fetal, 'Week 17', 'EK_PB_annov1', 'RNA')
cc_fetal_w18 <- performCellChatForDay(fetal, 'Week 18', 'EK_PB_annov1', 'RNA')

cc_fetal_w19 <- performCellChatForDay(fetal, 'Week 19', 'EK_PB_annov1', 'RNA')
cc_fetal_w20 <- performCellChatForDay(fetal, 'Week 20', 'EK_PB_annov1', 'RNA')
cc_fetal_w22 <- performCellChatForDay(fetal, 'Week 22', 'EK_PB_annov1', 'RNA')
cc_fetal_w24 <- performCellChatForDay(fetal, 'Week 24', 'EK_PB_annov1', 'RNA')
cc_fetal_w27 <- performCellChatForDay(fetal, 'Week 27', 'EK_PB_annov1', 'RNA')
cc_fetal_w9 <- performCellChatForDay(fetal, 'Week 9', 'EK_PB_annov1', 'RNA')

object.list <- list(FD59 = cc_fetal_d59, FD80 = cc_fetal_d80, FD82 = cc_fetal_d82,
                    FD105 = cc_fetal_d105, FD125 = cc_fetal_d125,
                      W9 = cc_fetal_w9, W11 = cc_fetal_w11,
                    W12 = cc_fetal_w12, W13 = cc_fetal_w13,
                    W14 = cc_fetal_w14, W15 = cc_fetal_w15,
                    W16 = cc_fetal_w16, W17 = cc_fetal_w17,
                    W18 = cc_fetal_w18, W19 = cc_fetal_w19,
                    W20 = cc_fetal_w20, W22 = cc_fetal_w22,
                    W24 = cc_fetal_w24, W27 = cc_fetal_w27)

MERGED_cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = T)

gg1 <- compareInteractions(MERGED_cellchat, show.legend = F, group = c(1:19))
gg2 <- compareInteractions(MERGED_cellchat, show.legend = F, group = c(1:19), measure = "weight")
gg1 + gg2
