## =====================================================
## libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(Matrix)
library(SeuratWrappers)
library(Banksy)
library(harmony)
library(cowplot)
library(rlang)
library(scales)
library(plyr)

## =====================================================
## Paths + input
in_dir  <- "/path/to/in_folder"
out_dir <- "/path/to/out_folder"
seu <- readRDS(file.path(in_dir, "ALL_Seurat_after_qc.rds"))

## =====================================================
## Normalization/transformation
options(future.globals.maxSize = 16000 * 1024^2)
# traditional way of normalizing data (using RC as normalization method and scale by median_nFeature)
median_nFeature <- median(seu$nFeature_RNA) 
seu <- NormalizeData(seu, normalization.method = "RC", scale.factor = median_nFeature)
all.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.genes)

## =====================================================
# add coordinates
# Adding CosMx coordinates as a global FOV 
cell_centroid_df = data.frame(y = seu$CenterX_global_px, x = seu$CenterY_global_px, cell = colnames(seu))
centroid_data = list("centroids" = CreateCentroids(cell_centroid_df))
coords = CreateFOV(coords = centroid_data, type = c("centroids"), assay = "RNA")
seu[["global"]] = coords
ImageDimPlot(seu, fov = "global", axes = TRUE, cols = "glasbey", size = 0.75, group.by = "PatientID")

# Adding CosMx coordinates as "reduction"  
coordinates = as.matrix(data.frame(CosMx_1 = seu$CenterX_global_px, CosMx_2 = seu$CenterY_global_px))
seu[['CosMx']] <- CreateDimReducObject(embeddings = coordinates, key = "CosMx_", global = T, assay = "RNA")

# save data
saveRDS(seu, file=file.path(out_dir, "ALL_Seurat_Norm.rds"))

## =====================================================
# neighbourhood expression construction
# Run Banksy (combines individual spot/cell expression and spatial-neighbourhood summaries into a seurat object)
# [this part wil replace your original scale.data] do not call scale.data after running banksy
# using global_fovID as group
# memory intensive
seu <- RunBanksy(seu, lambda = 0.2, # 0.2 for weaker smoothing, closer to 1 is for stronger spatial smoothing
                 verbose=TRUE, 
                 dimx = "CenterY_global_px",
                 dimy = "CenterX_global_px",
                 assay = "RNA", slot = "data", features = "all",
                 group = "global_fovID", # must specify otherwise the coordinates will overlap
                 split.scale = FALSE, # TRUE if you have a strong batch effect, but PCA won't run... don't know why
                 k_geom = 5) # how large a region you want to smooth: recommend 5-20
saveRDS(seu, file=file.path(out_dir, "ALL_Seurat_Banksy.rds"))

## =====================================================
# PCA
# using all genes
seu = RunPCA(seu, assay="BANKSY", features=rownames(seu))

DimPlot(seu, reduction = "pca", raster=F)
DimPlot(seu, reduction = "pca", group.by="runID")
DimPlot(seu, reduction = "pca", group.by="slideID")
DimPlot(seu, reduction = "pca", group.by="response")

DimPlot(seu, reduction = "pca", split.by="runID")
DimPlot(seu, reduction = "pca", split.by="slideID")
DimPlot(seu, reduction = "pca", split.by="response")

#Determining the number of PCs for UMAP
npcs = 25
ElbowPlot(seu, ndims = 50) + 
  geom_vline(xintercept = npcs, linetype="dotted", color = "red", size=1) + 
  geom_hline(yintercept = seu@reductions$pca@stdev[npcs], linetype="dotted", color = "red", size=1)

## =====================================================
# batch correction
# HARMONY
seu_harmony = RunHarmony(seu, group.by.vars = c("runID", "slideID"), plot_convergence = TRUE)

DimPlot(seu_harmony, reduction = "harmony")
DimPlot(seu_harmony, reduction = "harmony", group.by = "runID")
DimPlot(seu_harmony, reduction = "harmony", split.by = "runID")
DimPlot(seu_harmony, reduction = "harmony", split.by = "slideID")

# if you are happy with harmony
seu = seu_harmony
rm(seu_harmony)

saveRDS(seu, file=file.path(out_dir, "ALL_Seurat_Harmony.rds"))

## =====================================================
# UMAP and Clustering
set.seed(1234)
seu <- FindNeighbors(seu, reduction = "harmony", dims=1:20)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:20)

DefaultAssay(seu) <- "BANKSY"
seu <- FindClusters(seu, resolution = 0.5, algorithm = 4)

DimPlot(seu, reduction = "umap", label=T, cols="alphabet", raster=FALSE, group.by = "BANKSY_snn_res.0.5")
DimPlot(seu, reduction = "umap", label=T, cols="alphabet", raster=FALSE, split.by = "BANKSY_snn_res.0.5")

seu$BANKSY_snn_res.0.5 = factor(seu$BANKSY_snn_res.0.5, levels = c(1:13))

## =====================================================
## Find markers
Idents(seu) = seu$BANKSY_snn_res.0.5

DefaultAssay(seu) = "RNA"
seu.markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
seu.markers_df <- seu.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC) %>%
  print(n=10)
write.csv(seu.markers_df, file=file.path(out_dir, "ALL_Seurat_markers_r0.5.csv"))

## =====================================================
## Cell type annotation based on markers, and visualisation
seu$seurat_clusters <- seu$BANKSY_snn_res.0.5
Idents(seu) <- seu$seurat_clusters

# manual annotation
new_cluster_ids <- c(
  "1" = "Fib",
  "2" = "Mono/Mac",
  "3" = "B/Plasma", 
  "4" = "Epi",
  "5" = "Unassigned",
  "6" = "Epi",
  "7" = "Epi",
  "8" = "T",
  "9" = "Neu",
  "10" = "Epi",
  "11" = "Endo",
  "12" = "Epi",
  "13" = "B/Plasma"
)

# Re-annotate remaining clusters
seu$major_annot <- plyr::mapvalues(
  seu$seurat_clusters,
  from = names(new_cluster_ids),
  to = new_cluster_ids
)

# for light background
light_cluster_colors <- c(
  "#E41A1C",  # Red
  "#377EB8",  # Blue
  "#4DAF4A",  # Green
  "#984EA3",  # Purple
  "#FF7F00",  # Orange
  "#A65628",  # Brown
  "#F781BF",  # Pink
  "#999999",  # Grey
  "#66C2A5",  # Turquoise
  "#FC8D62",  # Coral
  "#8DA0CB",  # Steel Blue
  "#E78AC3",  # Light Pink
  "#A6D854",  # Lime Green
  "#FFD92F",  # Yellow
  "#E5C494",  # Tan
  "#B3B3B3",  # Light Grey
  "#1B9E77",  # Teal
  "#D95F02"   # Dark Orange
)

lvl <- c("Epi","Endo","Fib","Mono/Mac","Neu","T","B/Plasma","Unassigned")
seu$major_annot <- factor(seu$major_annot, levels = lvl)
DimPlot(seu, reduction = "umap", cols = light_cluster_colors,
        label=F, raster=FALSE, group.by = "major_annot", pt.size = 0.3)

# Dot plot
Idents(seu) = seu$major_annot
unique(seu$major_annot)

# Find markers
seu.markers_major_annot <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(seu.markers_major_annot, file=file.path(out_dir, "ALL_Seurat_markers_major_annot.csv"))

# top 5 genes according to log2FC
seu.markers_major_annot %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
DotPlot(seu, features = top5$gene, group.by = "major_annot") +
  scale_y_discrete(limits = rev(levels(seu$major_annot))) +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 12),
    axis.title = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10))

# custom gene list
custom_top5 <- list(
  "Epi" = c("KRT19", "EPCAM", "FGG", "KRT5", "CEACAM6"),
  "Endo" = c("PECAM1", "SPARCL1", "RGS5", "FLT1", "VWF"),
  "Fib" = c("COL1A1", "COL3A1", "COL1A2", "COL6A3", "THBS2"),
  "Mono/Mac" = c("LYZ", "C1QC", "APOE", "C1QB", "C1QA"),
  "Neu" = c("CXCL8", "HCAR2.3", "G0S2", "CCL3.L1.L3", "ITGAX"),
  "T" = c("LTB", "CD2", "CD3E", "SPOCK2", "PTPRCAP"),
  "B/Plasma" = c("IGHG2", "IGHG1", "IGKC", "IGHA1", "JCHAIN"),
  "Unassigned" = c("ACKR3", "TNFSF14", "ADGRV1", "THSD4", "CLOCK")
)

features_use <- unlist(custom_top5, use.names = FALSE)
DotPlot(seu, features = features_use, group.by = "major_annot", cols = "RdYlBu") +
  scale_y_discrete(limits = rev(levels(seu$major_annot))) +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 12),
    axis.title = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# without Unassigned
seu_plot <- subset(seu, subset = major_annot != "Unassigned")
seu_plot$major_annot <- droplevels(seu_plot$major_annot)

DotPlot(seu_plot, features = features_use, group.by = "major_annot", cols = "RdYlBu") +
  scale_y_discrete(limits = rev(levels(seu_plot$major_annot))) +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 12),
    axis.title = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

saveRDS(seu, file=file.path(out_dir, "ALL_Seurat_major_annot.rds"))
