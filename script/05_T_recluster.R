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

# subset T
T_seu <- subset(seu, subset = major_annot %in% "T")

## =====================================================
# re-normalisation and processing
DefaultAssay(T_seu) = "RNA"
median_nFeature <- median(T_seu$nFeature_RNA)
T_seu <- NormalizeData(T_seu, normalization.method = "RC", scale.factor = median_nFeature)
T_seu <- ScaleData(T_seu, features = rownames(T_seu))

# RunPCA
DefaultAssay(T_seu)
T_seu = RunPCA(T_seu, features=rownames(T_seu))
# HARMONY
T_seu = RunHarmony(T_seu, group.by.vars = c("runID","slideID"), plot_convergence = TRUE)

npcs = 15
ElbowPlot(T_seu, ndims = 50) +
  geom_vline(xintercept = npcs, linetype="dotted", color = "red", size=1) + 
  geom_hline(yintercept = T_seu@reductions$pca@stdev[npcs], linetype="dotted", color = "red", size=1)

# UMAP and reclustering
set.seed(1234)
DefaultAssay(T_seu) <- "RNA"
T_seu <- FindNeighbors(T_seu, reduction = "harmony", dims = 1:npcs)
T_seu <- RunUMAP(T_seu, reduction = "harmony", dims = 1:npcs)
T_seu <- FindClusters(T_seu, resolution = 0.5, algorithm = 4)
DimPlot(T_seu, reduction = "umap", label=T, cols="alphabet", raster=FALSE, group.by = "RNA_snn_res.0.5")

# find markers for each cluster
DefaultAssay(T_seu) <- "RNA"
Idents(T_seu) = T_seu$RNA_snn_res.0.5
T_seu.markers <- FindAllMarkers(T_seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
T_seu.markers_df <- T_seu.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC) %>%
  print(n=30)

## =====================================================
## annotation
T_seu$seurat_clusters <- T_seu$RNA_snn_res.0.5
Idents(T_seu) <- T_seu$seurat_clusters
T_clusters <- c("1", "4")
T_seu$T_subtype <- factor(ifelse(T_seu$seurat_clusters %in% T_clusters, "T", "Unassigned"),
                          levels = c("T", "Unassigned"))

saveRDS(T_seu, file.path(dir, "ALL_Seurat_T_Recluster_1.rds"))

## =====================================================
## Re-clustering
## removing "Unassigned"
T_only <- subset(T_seu, subset = T_subtype %in% "T")

## Re-normalize
DefaultAssay(T_only) <- "RNA"
median_nFeature <- median(T_only$nFeature_RNA)
T_only <- NormalizeData(T_only, normalization.method = "RC", scale.factor = median_nFeature)
T_only <- ScaleData(T_only, features = rownames(T_only))

# RunPCA
DefaultAssay(T_only)
T_only = RunPCA(T_only, features=rownames(T_only))
# HARMONY
T_only = RunHarmony(T_only, group.by.vars = c("runID","slideID"), plot_convergence = TRUE)

npcs = 10
ElbowPlot(T_only, ndims = 50) + 
  geom_vline(xintercept = npcs, linetype="dotted", color = "red", size=1) + 
  geom_hline(yintercept = T_only@reductions$pca@stdev[npcs], linetype="dotted", color = "red", size=1)

# UMAP and reclustering
set.seed(1234)
DefaultAssay(T_only) <- "RNA"
T_only <- FindNeighbors(T_only, reduction = "harmony", dims = 1:npcs)
T_only <- RunUMAP(T_only, reduction = "harmony", dims = 1:npcs)
T_only <- FindClusters(T_only, resolution = 0.4, algorithm = 4)
DimPlot(T_only, reduction = "umap", label=T, cols="alphabet", raster=FALSE, group.by = "RNA_snn_res.0.4")

# find markers for each cluster
DefaultAssay(T_only) <- "RNA"
Idents(T_only) = T_only$RNA_snn_res.0.4
T_only.markers <- FindAllMarkers(T_only, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
T_only.markers_df <- T_only.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC) %>%
  print(n=30)

# get top 5 markers per cluster
top5 <- T_only.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  pull(gene)

# add canonical markers
genes_use <- unique(c(top5, "CD4", "CD8A", "CD8B"))

DotPlot(T_only,features = genes_use,group.by = "RNA_snn_res.0.4",cols = "RdYlBu") +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12))

## =====================================================
## Re-annotation
T_only$seurat_clusters <- T_only$RNA_snn_res.0.4
Idents(T_only) <- T_only$seurat_clusters
new_cluster_ids <- c(
  "1" = "Activated_CD4T",
  "2" = "Cytotoxic_CD8T",
  "3" = "Cycling_like_T",
  "4" = "Treg"
)

T_only$T_subtype <- plyr::mapvalues(
  T_only$seurat_clusters,
  from = names(new_cluster_ids),
  to = new_cluster_ids
)
T_only$T_subtype <- factor(T_only$T_subtype, levels = c("Cytotoxic_CD8T", "Activated_CD4T", "Treg", "Cycling_like_T"))
palette.pbmc <- c("#b30024",
                  # "#e31a1c",
                  # "#ff7f00",
                  "#08519c",
                  # "#762a83",
                  "#f28ac6", 
                  "#20b2aa")
DimPlot(T_only, reduction = "umap", label=F, raster=FALSE, group.by = "T_subtype", pt.size = 0.4, cols = palette.pbmc)

# find markers for each cluster
DefaultAssay(T_only) <- "RNA"
Idents(T_only) = T_only$T_subtype
T_only.markers <- FindAllMarkers(T_only, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
T_only.markers_df <- T_only.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC) %>%
  print(n=30)
write.csv(T_only.markers, file=file.path(dir, "T_only_markers_annot.csv"))

## =====================================================
## Dotplot to show relevant markers
custom_markers <- list(
  Cytotoxic_CD8T       = c("GNLY", "GZMB", "PRF1", "NKG7", "GZMA"),
  Activated_CD4T       = c("TNFSF8", "FOS", "DUSP1", "CCR7", "CXCR4"),
  Treg                 = c("IL2RA", "TNFRSF4", "CTLA4", "TNFRSF18", "ICA1"),
  Cycling_like_T       = c("RPL22", "RPL32", "RPL34", "TPT1", "NACA")
)

features_final <- c("CD4", "CD8A", "CD8B", unique(unlist(custom_markers)))
DotPlot(T_only, features = features_final, group.by = "T_subtype", cols = "RdYlBu") +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12),
        axis.title = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

DotPlot(T_only, features = features_final, group.by = "T_subtype", cols = "RdYlBu") +
  RotatedAxis() +
  scale_y_discrete(limits = rev(levels(T_only$T_subtype))) +
  guides(
    colour = guide_colourbar(title = "Average Expression", title.position = "top",
                             title.hjust = 0.5, barwidth = 10, barheight = 0.6),
    size = guide_legend(title = "Percent Expressed", title.position = "top",
                        title.hjust = 0.5)
  ) +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 12),
    axis.title = element_blank(),
    legend.position = "bottom",        # move legends below
    legend.direction = "horizontal",   # make them horizontal
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

saveRDS(T_only, file=file.path(dir, "ALL_Seurat_T_Recluster_2.rds"))
