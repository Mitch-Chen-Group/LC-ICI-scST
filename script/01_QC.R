## QC

## =====================================================
# libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(Matrix)

in_dir  <- "/path/to/in_folder"
out_dir <- "/path/to/out_folder"

slideIDs <- c("RUN1_NR1","RUN1_NR2","RUN1_R1","RUN1_R2","RUN2_NR3","RUN2_NR4","RUN2_R3")

add_global_ids <- function(meta, expr, slide_folder) {

  parts <- strsplit(slide_folder, "_")[[1]]
  run_tag   <- parts[1]                          
  slide_tag <- paste(parts[-1], collapse = "_")  
  
  meta$runID   <- sub("^RUN", "Run", run_tag)    
  meta$slideID <- slide_tag                      
  
  # global IDs
  meta$global_cellID <- paste0(meta$runID, "_", meta$cell)
  meta$global_fovID  <- paste0(meta$slideID, "_", meta$fov)
  
  # apply to rownames (must match expr rownames)
  rownames(meta) <- meta$global_cellID
  rownames(expr) <- meta$global_cellID
  
  # sanity
  stopifnot(anyDuplicated(rownames(meta)) == 0)
  stopifnot(identical(rownames(meta), rownames(expr)))
  
  list(meta = meta, expr = expr)
}

read_one_slide <- function(slideID, in_dir) {
  slide_dir <- file.path(in_dir, slideID)
  
  meta_path <- file.path(slide_dir, paste0("metadata_for_Seurat_", slideID, ".RDS"))
  expr_path <- file.path(slide_dir, paste0("exprMat_for_Seurat_", slideID, ".RDS"))
  
  meta <- readRDS(meta_path)
  expr <- readRDS(expr_path)
  
  stopifnot(identical(rownames(meta), rownames(expr)))
  
  add_global_ids(meta, expr, slideID)
}

slide_list <- lapply(slideIDs, read_one_slide, in_dir = in_dir)
names(slide_list) <- slideIDs

metadata <- bind_rows(lapply(slide_list, `[[`, "meta"))
exprMat  <- do.call(rbind, lapply(slide_list, \(x) as.matrix(x$expr)))

stopifnot(identical(rownames(metadata), rownames(exprMat)))

saveRDS(metadata, file = file.path(out_dir, "metadata_ALL_before_qc.rds"))
saveRDS(exprMat,  file = file.path(out_dir, "exprMat_ALL_before_qc.rds"))

## =====================================================
##  filter unselected FOVs
metadata <- metadata %>%
  filter(!is.na(region) == TRUE)
## filter FOV border cells (margin = mean cell diameter in px)
mean(metadata$Area.um2) # all = 97.19309
mean_diameter_um <- 2 * sqrt(mean(metadata$Area.um2, na.rm = TRUE) / pi)
mean_diameter_px <- mean_diameter_um / 0.12028 # um per pixel in CosMx
metadata = metadata %>%
  filter(
    CenterX_local_px > mean_diameter_px,
    CenterX_local_px < 4256 - mean_diameter_px,
    CenterY_local_px > mean_diameter_px,
    CenterY_local_px < 4256 - mean_diameter_px) # 500,846

# filter exprMat
exprMat <- exprMat[rownames(metadata), , drop = FALSE]
stopifnot(identical(rownames(metadata), rownames(exprMat)))

## =====================================================
## Quick visual checks
# Size (Area) 
ggplot(metadata %>% filter(PatientID =="P1"), aes(x=CenterX_global_px, y=CenterY_global_px, colour=log2(Area))) + 
  geom_point(size=0.005) + 
  theme_minimal() + scale_colour_gradientn(colours = c("grey50", "blue", "yellow", "red")) + 
  theme(legend.position="right") + ggtitle("Area (log2)")
# DAPI staining
ggplot(metadata %>% filter(PatientID =="P1"), aes(x=CenterX_global_px, y=CenterY_global_px, colour=log2(Mean.DAPI))) + 
  geom_point(size=0.005) + 
  theme_minimal() + scale_colour_gradientn(colours = c("grey50", "blue", "yellow", "red")) + 
  theme(legend.position="right") + ggtitle("Mean.DAPI (log2)")
# CD45 staining
ggplot(metadata %>% filter(PatientID =="P1"), aes(x=CenterX_global_px, y=CenterY_global_px, colour=log2(Mean.CD45))) + 
  geom_point(size=0.01) + 
  theme_minimal() + scale_colour_gradientn(colours = c("grey50", "blue", "yellow", "red")) + 
  theme(legend.position="right") + ggtitle("Mean.CD45 (log2)")
# CD68 staining
ggplot(metadata %>% filter(PatientID =="P1"), aes(x=CenterX_global_px, y=CenterY_global_px, colour=log2(Mean.CD68))) + 
  geom_point(size=0.01) + 
  theme_minimal() + scale_colour_gradientn(colours = c("grey50", "blue", "yellow", "red")) + 
  theme(legend.position="right") + ggtitle("Mean.CD68 (log2)")
# PanCK staining
ggplot(metadata %>% filter(PatientID =="P1"), aes(x=CenterX_global_px, y=CenterY_global_px, colour=log2(Mean.PanCK))) + 
  geom_point(size=0.01) + 
  theme_minimal() + scale_colour_gradientn(colours = c("grey50", "blue", "yellow", "red")) + 
  theme(legend.position="right") + ggtitle("Mean.PanCK (log2)")

## =====================================================
## Create Seurat
seu <- CreateSeuratObject(
  counts = t(exprMat),
  project = "CosMx",
  meta.data = metadata,
  min.cells = 3,
  min.features = 1
)
seu

slideID <- c("NR1","NR2","R1","R2","NR3","NR4","R3")
seu$slideID <- factor(seu$slideID, levels = slideID)

unique(seu$PatientID)
seu$PatientID = factor(seu$PatientID, levels = c("P1","P2","P3","P4","P5","P6", # non-responders
                                                 "P7","P8","P9","P10")) # responders

seu$response <- ifelse(seu$PatientID %in% c("P1","P2","P3","P4","P5","P6"), "Non-responders","Responders")
seu$response <- factor(seu$response,levels = c("Non-responders", "Responders"))

seu$response_lbl <- ifelse(seu$PatientID %in% c("P1","P2","P3","P4","P5","P6"), "PD","PR")
seu$response_lbl <- factor(seu$response_lbl,levels = c("PD", "PR"))

# optional 
saveRDS(seu, file = file.path(out_dir, "Seurat_ALL_before_qc.rds"))

## =====================================================
## Overview
seu[["percent.NegPrb"]] <- PercentageFeatureSet(seu, pattern = "^NegPrb")

# Number of FOVs per slide
seu@meta.data %>%
  dplyr::select(slideID, global_fovID) %>%
  dplyr::distinct() %>%
  dplyr::count(slideID, name = "num_FOVs")

# Number of FOVs per PatientID
seu@meta.data %>%
  dplyr::select(PatientID, fov) %>%
  dplyr::distinct() %>%
  dplyr::count(PatientID, name = "num_FOVs")

VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.NegPrb"), ncol = 3, pt.size = 0)

VlnPlot(seu, features = "nCount_RNA", pt.size = 0, group.by = "slideID") +
  stat_boxplot(geom = "boxplot", width = 0.05, fill = "white", colour = "black") +
  theme_classic()

# Summary stats by slide
seu@meta.data %>%
  group_by(slideID) %>%
  summarise(
    mean_nCount_RNA    = mean(nCount_RNA),
    mean_nFeature_RNA  = mean(nFeature_RNA),
    mean_percent_NegPrb= mean(percent.NegPrb),
    n_cells            = dplyr::n(),
    .groups = "drop"
  )

hist(seu$nFeature_RNA, breaks = 100, main = "nFeature_RNA")
abline(v = 10, col = "red", lwd = 3)

hist(seu$percent.NegPrb, breaks = 100, main = "percent.NegPrb")
abline(v = 5, col = "red", lwd = 3)

plot1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.NegPrb", raster = F)
plot2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = F)
plot1 + plot2

## =====================================================
## QC
## Apply filters
seu <- subset(seu, subset = nFeature_RNA >= 10 & percent.NegPrb <= 5)
seu <- subset(seu, subset = nCount_RNA >= 20)

## =====================================================
## Doublet removal (scDblFinder)
library(SingleCellExperiment)
library(scDblFinder)
library(BiocParallel)

set.seed(1234)

sce <- as.SingleCellExperiment(seu)
t <- system.time({
  sce <- scDblFinder(
    sce,
    samples = "PatientID",
    dbr = 0.05,
    dbr.sd = 1,
    BPPARAM = SerialParam(RNGseed = 1234)
  )
})

t

seu$scDblFinder_class <- sce$scDblFinder.class
seu$scDblFinder_score <- sce$scDblFinder.score

# visual before subsetting
FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = FALSE, pt.size = 0.05)
# subset
seu <- subset(seu, subset = scDblFinder_class == "singlet")
# visual after subsetting
FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = FALSE, pt.size = 0.05)

## =====================================================
## Background gene + NegPrb removal
NegPrb <- c("NegPrb1","NegPrb2","NegPrb3","NegPrb4","NegPrb5",
            "NegPrb6","NegPrb7","NegPrb8","NegPrb9","NegPrb10")

slide_col <- "slideID"
remove_if_below_in_more_than <- 4

X <- GetAssayData(seu, assay = "RNA", layer = "counts")  # features x cells

neg   <- NegPrb
genes <- setdiff(rownames(X), neg)

slide_id <- seu@meta.data[[slide_col]]
cells_by_slide <- split(colnames(X), slide_id)

below_count <- setNames(integer(length(genes)), genes)

for (s in names(cells_by_slide)) {
  cells_s <- cells_by_slide[[s]]
  
  bg_s <- mean(X[neg, cells_s, drop = FALSE])  # slide-specific NegPrb background
  gene_mean_s <- Matrix::rowMeans(X[genes, cells_s, drop = FALSE])
  
  below_count <- below_count + as.integer(gene_mean_s <= bg_s)
}

genes_to_remove <- names(below_count)[below_count > remove_if_below_in_more_than]
genes_to_keep   <- setdiff(genes, genes_to_remove)

head(sort(below_count[genes_to_remove], decreasing = TRUE), 30)

# sanity check
seu <- subset(seu, features = genes_to_keep)

saveRDS(seu, file = file.path(out_dir, "ALL_Seurat_after_qc.rds"))
