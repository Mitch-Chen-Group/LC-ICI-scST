## =====================================================
## Libraries
library(Seurat)
library(CellChat)
library(patchwork)
library(circlize)
library(ComplexHeatmap)
library(presto)
library(dplyr)

## =====================================================
## Load data
base_dir <- "/path/to/folder"
seu <- readRDS(file.path(base_dir, "ALL_Seurat_annot_combined_xP7.rds"))

## =====================================================
## Define groups + Idents
Idents(seu) <- seu$annot_combined
seu <- subset(seu, subset = annot_combined != "Unassigned")
seu$annot_combined <- droplevels(seu$annot_combined)

## =====================================================
## Split Seurat
seu_PR <- subset(seu, subset = response_lbl == "PR")
seu_PD <- subset(seu, subset = response_lbl == "PD")
stopifnot(ncol(seu_PR) > 0, ncol(seu_PD) > 0)

CellChatDB <- CellChatDB.human

saveRDS(seu_PD, file.path(base_dir, "seu_PD_xP7.rds"))
saveRDS(seu_PR, file.path(base_dir, "seu_PR_xP7.rds"))

## =====================================================
## PR CellChat (spatial mode)
data.input <- Seurat::GetAssayData(seu_PR, slot = "data", assay = "RNA")

meta <- data.frame(
  labels = Seurat::Idents(seu_PR),
  samples = seu_PR$PatientID,
  row.names = names(Seurat::Idents(seu_PR))
)
meta$samples <- factor(meta$samples)

spatial.locs <- data.frame(
  x = seu_PR@meta.data$CenterY_global_px,
  y = seu_PR@meta.data$CenterX_global_px
)
rownames(spatial.locs) <- rownames(seu_PR@meta.data)
spatial.locs <- as.matrix(spatial.locs)

conversion.factor <- 0.12028
d <- computeCellDistance(spatial.locs) # memory intensive
spot.size <- min(d) * conversion.factor
spatial.factors <- data.frame(ratio = conversion.factor, tol = spot.size/2)

cellchat_PR <- createCellChat(
  object = data.input, meta = meta, group.by = "labels",
  datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors
)

CellChatDB.use <- subsetDB(CellChatDB)
cellchat_PR@DB <- CellChatDB.use
cellchat_PR <- subsetData(cellchat_PR)
cellchat_PR <- identifyOverExpressedGenes(cellchat_PR)
cellchat_PR <- identifyOverExpressedInteractions(cellchat_PR, variable.both = FALSE)

cellchat_PR <- computeCommunProb(
  cellchat_PR,
  type = "truncatedMean", trim = 0.1,
  distance.use = TRUE, interaction.range = 250,
  scale.distance = 4.5,
  contact.dependent = TRUE, contact.range = 20
)

cellchat_PR <- computeCommunProbPathway(cellchat_PR)
cellchat_PR <- aggregateNet(cellchat_PR)
cellchat_PR <- netAnalysis_computeCentrality(cellchat_PR, slot.name = "netP")
saveRDS(cellchat_PR, file = file.path(base_dir, "cellchat_PR_spatial.rds"))

## =====================================================
## PD CellChat (spatial mode)
data.input <- Seurat::GetAssayData(seu_PD, layer = "data", assay = "RNA")

meta <- data.frame(
  labels = Seurat::Idents(seu_PD),
  samples = seu_PD$PatientID,
  row.names = names(Seurat::Idents(seu_PD))
)
meta$samples <- factor(meta$samples)

spatial.locs <- data.frame(
  x = seu_PD@meta.data$CenterY_global_px,
  y = seu_PD@meta.data$CenterX_global_px
)
rownames(spatial.locs) <- rownames(seu_PD@meta.data)
spatial.locs <- as.matrix(spatial.locs)

conversion.factor <- 0.12028
d <- computeCellDistance(spatial.locs) # memory intensive
spot.size <- min(d) * conversion.factor
spatial.factors <- data.frame(ratio = conversion.factor, tol = spot.size/2)

cellchat_PD <- createCellChat(
  object = data.input, meta = meta, group.by = "labels",
  datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors
)

CellChatDB.use <- subsetDB(CellChatDB)
cellchat_PD@DB <- CellChatDB.use
cellchat_PD <- subsetData(cellchat_PD)
cellchat_PD <- identifyOverExpressedGenes(cellchat_PD)
cellchat_PD <- identifyOverExpressedInteractions(cellchat_PD, variable.both = FALSE)

cellchat_PD <- computeCommunProb(
  cellchat_PD,
  type = "truncatedMean", trim = 0.1,
  distance.use = TRUE, interaction.range = 250,
  scale.distance = 4.5,
  contact.dependent = TRUE, contact.range = 20
)

cellchat_PD <- computeCommunProbPathway(cellchat_PD)
cellchat_PD <- aggregateNet(cellchat_PD)
cellchat_PD <- netAnalysis_computeCentrality(cellchat_PD, slot.name = "netP")
saveRDS(cellchat_PD, file = file.path(base_dir, "cellchat_PD_spatial.rds"))

## =====================================================
## Merge PR vs PD
object.list <- list(PR = cellchat_PR, PD = cellchat_PD)
cellchat_PRvsPD <- mergeCellChat(object.list, add.names = names(object.list))
saveRDS(cellchat_PRvsPD, file = file.path(base_dir, "cellchat_PRvsPD_spatial_merged.rds"))

## =====================================================
cellchat_PR <- readRDS(file = file.path(base_dir, "cellchat_PR_xP7_spatial.rds"))
cellchat_PD <- readRDS(file = file.path(base_dir, "cellchat_PD_xP7_spatial.rds"))
cellchat_PRvsPD <- readRDS(file = file.path(base_dir, "cellchat_PRvsPD_xP7_spatial_merged.rds"))

## =====================================================
## overall interactions (count + strength)
p_count  <- compareInteractions(cellchat_PRvsPD, show.legend = TRUE, group = names(object.list))
p_weight <- compareInteractions(cellchat_PRvsPD, show.legend = TRUE, measure = "weight", group = names(object.list))
p_count + p_weight

## =====================================================
## incoming vs outgoing scatter (PR and PD)
num.link <- sapply(object.list, function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})
weight.MinMax <- c(min(num.link), max(num.link))

gg <- list()
for (i in seq_along(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(
    object.list[[i]],
    title = names(object.list)[i],
    weight.MinMax = weight.MinMax
  )
}
patchwork::wrap_plots(plots = gg)


# Identify altered signaling with distinct interaction strength
# Compare the overall information flow of each signaling pathway or ligand-receptor pair
gg1 <- rankNet(cellchat_PRvsPD, slot.name = "netP", mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat_PRvsPD, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
gg1 + gg2

## =====================================================
## bubble plots (SPP1_Mac <-> Epi), increased/decreased in PD
# Note: object.list order is PR (1) then PD (2).
# Increased in PD => max.dataset = 2, Decreased in PD => max.dataset = 1.

# SPP1_Mac -> Epi
sources <- c("SPP1_APOE_Mac")
targets <- c("Epi")
ggC_inc <- netVisual_bubble(
  cellchat_PRvsPD,
  sources.use = sources, targets.use = targets,
  comparison = c(1, 2),
  max.dataset = 2,
  title.name = "Increased signaling in PD",
  angle.x = 45, remove.isolate = TRUE
)
ggC_dec <- netVisual_bubble(
  cellchat_PRvsPD,
  sources.use = sources, targets.use = targets,
  comparison = c(1, 2),
  max.dataset = 1,
  title.name = "Decreased signaling in PD",
  angle.x = 45, remove.isolate = TRUE
)
ggC_inc + ggC_dec

# Epi -> SPP1_Mac
sources <- c("Epi")
targets <- c("SPP1_APOE_Mac")
ggD_inc <- netVisual_bubble(
  cellchat_PRvsPD,
  sources.use = sources, targets.use = targets,
  comparison = c(1, 2),
  max.dataset = 2,
  title.name = "Increased signaling in PD",
  angle.x = 45, remove.isolate = TRUE
)
ggD_dec <- netVisual_bubble(
  cellchat_PRvsPD,
  sources.use = sources, targets.use = targets,
  comparison = c(1, 2),
  max.dataset = 1,
  title.name = "Decreased signaling in PD",
  angle.x = 45, remove.isolate = TRUE
)
ggD_inc + ggD_dec

## =====================================================
## chord diagrams for SPP1 pathway in PD vs PR
pathways.show <- "SPP1"
par(mfrow = c(1,2), xpd = TRUE)
netVisual_aggregate(cellchat_PD, signaling = pathways.show, layout = "chord",
                    signaling.name = "SPP1 Signalling (PD)", vertex.weight = NULL)
netVisual_aggregate(cellchat_PR, signaling = pathways.show, layout = "chord",
                    signaling.name = "SPP1 Signalling (PR)", vertex.weight = NULL)
par(mfrow = c(1,1))
