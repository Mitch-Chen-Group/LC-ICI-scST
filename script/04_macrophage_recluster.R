## =====================================================
## libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)

## =====================================================
## Paths + input
dir <- "/path/to/folder"
seu <- readRDS(file.path(dir, "ALL_Seurat_major_annot.rds"))

# subset Mono/Mac
mac_seu <- subset(seu, subset = major_annot %in% "Mono/Mac")

## =====================================================
# re-normalisation and processing
DefaultAssay(mac_seu) = "RNA"
median_nFeature <- median(mac_seu$nFeature_RNA)
mac_seu <- NormalizeData(mac_seu, normalization.method = "RC", scale.factor = median_nFeature)
mac_seu <- ScaleData(mac_seu, features = rownames(mac_seu))

# RunPCA
DefaultAssay(mac_seu)
mac_seu = RunPCA(mac_seu, features=rownames(mac_seu))
# HARMONY
mac_seu = RunHarmony(mac_seu, group.by.vars = c("runID","slideID"), plot_convergence = TRUE)

npcs = 15
ElbowPlot(mac_seu, ndims = 50) +
  geom_vline(xintercept = npcs, linetype="dotted", color = "red", size=1) + 
  geom_hline(yintercept = mac_seu@reductions$pca@stdev[npcs], linetype="dotted", color = "red", size=1)

# UMAP and reclustering
set.seed(1234)
DefaultAssay(mac_seu) <- "RNA"
mac_seu <- FindNeighbors(mac_seu, reduction = "harmony", dims = 1:npcs)
mac_seu <- RunUMAP(mac_seu, reduction = "harmony", dims = 1:npcs)
mac_seu <- FindClusters(mac_seu, resolution = 0.5, algorithm = 4)
DimPlot(mac_seu, reduction = "umap", label=T, cols="alphabet", raster=FALSE, group.by = "RNA_snn_res.0.5")

# find markers for each cluster
DefaultAssay(mac_seu) <- "RNA"
Idents(mac_seu) = mac_seu$RNA_snn_res.0.5
mac_seu.markers <- FindAllMarkers(mac_seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
mac_seu.markers_df <- mac_seu.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC) %>%
  print(n=30)

## =====================================================
## annotation
mac_seu$seurat_clusters <- mac_seu$RNA_snn_res.0.5
Idents(mac_seu) <- mac_seu$seurat_clusters
mac_clusters <- c("1", "4", "5", "6")
mac_seu$mac_subtype <- factor(ifelse(mac_seu$seurat_clusters %in% mac_clusters, "Mac", "Unassigned"),
                              levels = c("Mac", "Unassigned"))

saveRDS(mac_seu, file.path(dir, "ALL_Seurat_Mac_Recluster_1.rds"))

## =====================================================
## Re-clustering
## removing "Unassigned"
mac_seu <- readRDS(file.path(dir, "ALL_Seurat_Mac_Recluster_1.rds"))
mac_only <- subset(mac_seu, subset = mac_subtype %in% "Mac")

## Re-normalize
DefaultAssay(mac_only) <- "RNA"
median_nFeature <- median(mac_only$nFeature_RNA)
mac_only <- NormalizeData(mac_only, normalization.method = "RC", scale.factor = median_nFeature)
mac_only <- ScaleData(mac_only, features = rownames(mac_only))

# RunPCA
DefaultAssay(mac_only)
mac_only = RunPCA(mac_only, features=rownames(mac_only))
# HARMONY
mac_only = RunHarmony(mac_only, group.by.vars = c("runID","slideID"), plot_convergence = TRUE)

npcs = 15
ElbowPlot(mac_only, ndims = 50) +
  geom_vline(xintercept = npcs, linetype="dotted", color = "red", size=1) + 
  geom_hline(yintercept = mac_only@reductions$pca@stdev[npcs], linetype="dotted", color = "red", size=1)

# UMAP and reclustering
set.seed(1234)
DefaultAssay(mac_only) <- "RNA"
mac_only <- FindNeighbors(mac_only, reduction = "harmony", dims = 1:npcs)
mac_only <- RunUMAP(mac_only, reduction = "harmony", dims = 1:npcs)
mac_only <- FindClusters(mac_only, resolution = 0.3, algorithm = 4)
DimPlot(mac_only, reduction = "umap", label=T, cols="alphabet", raster=FALSE, group.by = "RNA_snn_res.0.3")

# find markers for each cluster
DefaultAssay(mac_only) <- "RNA"
Idents(mac_only) = mac_only$RNA_snn_res.0.3
mac_only.markers <- FindAllMarkers(mac_only, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
mac_only.markers_df <- mac_only.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC) %>%
  print(n=30)

## =====================================================
## Re-annotation
mac_only$seurat_clusters <- mac_only$RNA_snn_res.0.3
Idents(mac_only) <- mac_only$seurat_clusters
new_cluster_ids <- c(
  "1" = "MMP19_Mac",
  "2" = "SPP1_Mac",
  "3" = "CXCL9_Mac",
  "4" = "MARCO_Mac"
)

mac_only$mac_subtype <- plyr::mapvalues(
  mac_only$seurat_clusters,
  from = names(new_cluster_ids),
  to = new_cluster_ids
)
mac_only$mac_subtype <- factor(mac_only$mac_subtype, levels = c("MMP19_Mac", "SPP1_Mac", "MARCO_Mac", "CXCL9_Mac"))
DimPlot(mac_only, reduction = "umap", label=F, raster=FALSE, group.by = "mac_subtype", pt.size = 0.3)

# find markers for each cluster
DefaultAssay(mac_only) <- "RNA"
Idents(mac_only) = mac_only$mac_subtype
mac_only.markers <- FindAllMarkers(mac_only, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
mac_only.markers_df <- mac_only.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC) %>%
  print(n=30)
write.csv(mac_only.markers, file=file.path(dir, "mac_only_markers_annot.csv"))

## =====================================================
## Dotplot to show relevant markers
custom_markers <- list(
  MMP19_Mac       = c("MMP19", "ZFP36", "FOS", "RGS2", "DUSP1"),
  SPP1_Mac       = c("SPP1", "CHI3L1", "APOE", "APOC1", "GPNMB"),
  MARCO_Mac     = c("FABP4", "MARCO", "FN1", "S100A4", "ALOX5AP"),
  CXCL9_Mac     = c("CCL13", "CXCL9", "HLA.DQA1", "SLC40A1", "SELENOP")
)

features_final <- unique(unlist(custom_markers))
DotPlot(mac_only, features = features_final, group.by = "mac_subtype", cols = "RdYlBu") +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12),
        axis.title = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

saveRDS(mac_only, file=file.path(dir, "ALL_Seurat_Mac_Recluster_2.rds"))
