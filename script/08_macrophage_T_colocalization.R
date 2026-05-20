## =============================================================================
## Colocalization / mixing score workflow
## https://github.com/SGIlabes/CosMx-BRCA-/blob/main/TNBC/06.Colocalization_with_Mac.R

## =============================================================================
## Libraries
library(dplyr)
library(tidyr)
library(Seurat)
library(SingleCellExperiment)
library(SpatialExperiment)
library(SPIAT)
library(reshape2)
library(pheatmap)
library(ggplot2)
library(ggpubr)
library(stringr)
library(rstatix)

## =============================================================================
## Paths + input
base_dir <- "path/to/folder"
seu_path <- file.path(base_dir, "ALL_Seurat_annot_combined_mac_T.rds")

seu <- readRDS(seu_path)
md  <- seu@meta.data

## =============================================================================
## Basic checks
required_cols <- c(
  "CenterX_local_px", "CenterY_local_px", "global_fovID",
  "PatientID", "response", "response_lbl", "annot_combined"
)
stopifnot(all(required_cols %in% colnames(md)))

## Convert coords from px to um
md <- md %>%
  mutate(
    CenterX_um = CenterX_local_px * 0.12028,
    CenterY_um = CenterY_local_px * 0.12028
  )

## =============================================================================
## Build SpatialExperiment for SPIAT
sce <- Seurat::as.SingleCellExperiment(seu)

spatialCoords <- as.matrix(md[, c("CenterX_um", "CenterY_um")])
colnames(spatialCoords) <- c("Cell.X.Position", "Cell.Y.Position")

spe <- SpatialExperiment(
  assays        = assays(sce),
  rowData       = rowData(sce),
  colData       = colData(sce),
  metadata      = metadata(sce),
  reducedDims   = reducedDims(sce),
  altExps       = altExps(sce),
  sample_id     = "cosmx",
  spatialCoords = spatialCoords
)

## Add required SPIAT columns
spe$in_tissue       <- 1
spe$data_id         <- "cosmx"
spe$fov_id          <- md$global_fovID
spe$annot_combined  <- as.character(md$annot_combined)

stopifnot(!any(is.na(spe$fov_id)))
stopifnot(!any(is.na(spe$annot_combined)))

fov_ids <- unique(spe$fov_id)

## =============================================================================
## FOV-level metadata for joining
fov_info <- md %>%
  dplyr::select(global_fovID, PatientID, response, response_lbl) %>%
  dplyr::distinct() %>%
  dplyr::rename(fov_id = global_fovID)

## =============================================================================
## Define target groups
mac_targets <- c("MMP19_Mac", "SPP1_Mac", "MARCO_Mac", "CXCL9_Mac")
t_targets   <- c("Cytotoxic_CD8T", "Activated_CD4T", "Treg", "Ribosome_T")

radii <- c(10, 25, 50, 75, 100)

## =============================================================================
## Function: compute mixing scores
compute_mixing_scores <- function(spe, fov_ids, targets, fov_info, radii, feature_col = "annot_combined") {
  
  annot_values <- unique(as.character(colData(spe)[[feature_col]]))
  ref_levels <- setdiff(annot_values, c("Unassigned", targets))
  ref_levels <- ref_levels[!is.na(ref_levels)]
  
  cmixing <- tibble()
  
  for (id in fov_ids) {
    spe_tmp <- spe[, spe$fov_id == id]
    
    if (ncol(spe_tmp) == 0) next
    spe_tmp$sample_id <- as.character(id)
    
    mix_list <- lapply(radii, function(r) {
      out <- tryCatch(
        mixing_score_summary(
          spe_object         = spe_tmp,
          reference_celltype = ref_levels,
          target_celltype    = targets,
          radius             = r,
          feature_colname    = feature_col
        ),
        error = function(e) NULL
      )
      
      if (is.null(out)) return(NULL)
      
      out$radius <- r
      out$fov_id <- id
      out
    })
    
    mix_all <- bind_rows(mix_list)
    if (nrow(mix_all) > 0) {
      cmixing <- bind_rows(cmixing, mix_all)
    }
  }
  
  cmixing_clean <- cmixing %>%
    left_join(fov_info, by = "fov_id") %>%
    filter(Number_of_reference_cells > 0) %>%
    mutate(
      Normalised_mixing_score = replace(Normalised_mixing_score, is.na(Normalised_mixing_score), 0)
    )
  
  return(cmixing_clean)
}

## =============================================================================
## Run macrophage and T-cell analyses separately
cmix_mac <- compute_mixing_scores(
  spe      = spe,
  fov_ids  = fov_ids,
  targets  = mac_targets,
  fov_info = fov_info,
  radii    = radii
)

cmix_t <- compute_mixing_scores(
  spe      = spe,
  fov_ids  = fov_ids,
  targets  = t_targets,
  fov_info = fov_info,
  radii    = radii
)

## =============================================================================
## Save outputs
write.csv(cmix_mac, file.path(base_dir, "ALL_Mac_targets_colocalization_score.csv"), row.names = FALSE)
write.csv(cmix_t,   file.path(base_dir, "ALL_T_targets_colocalization_score.csv"),   row.names = FALSE)

## Optional: read back in
cmix_mac <- read.csv(file.path(base_dir, "ALL_Mac_targets_colocalization_score.csv"),stringsAsFactors = FALSE)
cmix_t <- read.csv(file.path(base_dir, "ALL_T_targets_colocalization_score.csv"),stringsAsFactors = FALSE)

## =============================================================================
## Function: average heatmap
plot_average_heatmap <- function(data,
                                 radius_use = 50,
                                 patient_use = NULL,
                                 response_use = NULL,
                                 plot_title = NULL) {
  
  cmix_sub <- data %>%
    filter(radius == radius_use)
  
  if (!is.null(patient_use)) {
    cmix_sub <- cmix_sub %>%
      filter(PatientID %in% patient_use)
  }
  
  if (!is.null(response_use)) {
    cmix_sub <- cmix_sub %>%
      filter(response_lbl %in% response_use)
  }
  
  if (nrow(cmix_sub) == 0) {
    stop("No rows found for the selected filters.")
  }
  
  data_avg <- cmix_sub %>%
    group_by(Target, Reference) %>%
    summarise(
      AverageColocalization = mean(Normalised_mixing_score, na.rm = TRUE),
      .groups = "drop"
    )
  
  data_matrix <- reshape2::dcast(
    data_avg,
    Target ~ Reference,
    value.var = "AverageColocalization"
  )
  
  rownames(data_matrix) <- data_matrix$Target
  data_matrix <- data_matrix %>%
    dplyr::select(-Target) %>%
    as.matrix()
  
  data_matrix[is.na(data_matrix)] <- 0
  data_matrix[data_matrix < 0] <- 0
  data_matrix[data_matrix > 1] <- 1
  
  if (is.null(plot_title)) {
    plot_title <- paste0("Average mixing score (", radius_use, " µm)")
  }
  
  pheatmap(
    data_matrix,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    main = plot_title,
    color = colorRampPalette(c("#08306b", "#1262b3", "white", "#cc3d3d", "#67000d"))(100)
  )
}

## =============================================================================
## Function: boxplot for either macrophage or T-cell targets
create_boxplot <- function(data,
                           reference,
                           target_levels,
                           y_label,
                           stat_label_y,
                           ylim = NULL,
                           radius_use = 50,
                           patient_use = NULL,
                           response_use = NULL,
                           x_text_angle = 30) {
  
  plot_data <- data %>%
    filter(
      radius == radius_use,
      Reference %in% reference,
      Target %in% target_levels
    )
  
  if (!is.null(patient_use)) {
    plot_data <- plot_data %>%
      filter(PatientID %in% patient_use)
  }
  
  if (!is.null(response_use)) {
    plot_data <- plot_data %>%
      filter(response_lbl %in% response_use)
  }
  
  plot_data$Target <- factor(plot_data$Target, levels = target_levels)
  
  p <- ggplot(plot_data, aes(x = Target, y = Normalised_mixing_score, fill = Target)) +
    geom_boxplot(outlier.shape = NA) +
    theme_minimal() +
    theme(
      axis.line = element_line(),
      axis.title = element_text(face = "bold"),
      legend.position = "none",
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(size = 10, angle = x_text_angle, hjust = 1),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = "",
      y = y_label,
      fill = ""
    ) +
    stat_compare_means(
      method = "kruskal.test",
      label.y = stat_label_y
    )
  
  if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }
  
  print(p)
}

## =============================================================================
# compare heatmap by response groups
plot_heatmap_by_response <- function(data,
                                     radius_use = 50,
                                     patient_use = NULL,
                                     response_levels = c("PD", "PR"),
                                     plot_title = NULL) {
  
  cmix_sub <- data %>%
    filter(radius == radius_use)
  
  if (!is.null(patient_use)) {
    cmix_sub <- cmix_sub %>%
      filter(PatientID %in% patient_use)
  }
  
  cmix_sub <- cmix_sub %>%
    filter(response_lbl %in% response_levels)
  
  if (nrow(cmix_sub) == 0) {
    stop("No rows found for the selected filters.")
  }
  
  data_avg <- cmix_sub %>%
    group_by(Target, response_lbl, Reference) %>%
    summarise(
      AverageColocalization = mean(Normalised_mixing_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Target_Response = paste0(Target, " (", response_lbl, ")")
    )
  
  data_matrix <- reshape2::dcast(
    data_avg,
    Target_Response ~ Reference,
    value.var = "AverageColocalization"
  )
  
  rownames(data_matrix) <- data_matrix$Target_Response
  data_matrix <- data_matrix %>%
    dplyr::select(-Target_Response) %>%
    as.matrix()
  
  data_matrix[is.na(data_matrix)] <- 0
  
  if (is.null(plot_title)) {
    plot_title <- paste0("Average mixing score by response (", radius_use, " µm)")
  }
  
  pheatmap(
    data_matrix,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    main = plot_title,
    color = colorRampPalette(c("#08306b", "#1262b3", "white", "#cc3d3d", "#67000d"))(100)
  )
}

## =============================================================================
# compare boxplots by response
plot_boxplot_by_response <- function(data,
                                     reference,
                                     target_levels,
                                     y_label,
                                     stat_label_y = NULL,
                                     ylim = NULL,
                                     radius_use = 50,
                                     patient_use = NULL,
                                     response_levels = c("PD", "PR"),
                                     x_text_angle = 30) {
  
  plot_data <- data %>%
    filter(radius == radius_use,
           Reference %in% reference,
           Target %in% target_levels,
           response_lbl %in% response_levels)
  
  if (!is.null(patient_use)) {
    plot_data <- plot_data %>%
      filter(PatientID %in% patient_use)
  }
  
  plot_data$Target <- factor(plot_data$Target,
                             levels = target_levels)
  
  p <- ggplot(plot_data,
              aes(x = Target,
                  y = Normalised_mixing_score,
                  fill = response_lbl)) +
    geom_boxplot(position = position_dodge(width = 0.8),
                 outlier.shape = NA) +
    stat_compare_means(aes(group = response_lbl),
                       method = "wilcox.test",
                       label = "p.format",
                       label.y = stat_label_y) +
    theme_minimal() +
    theme(axis.line = element_line(),
          axis.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(size = 10,
                                     angle = x_text_angle,
                                     hjust = 1),
          panel.grid.minor = element_blank()) +
    labs(x = "",
         y = y_label,
         fill = "") +
    scale_fill_manual(values = c("PD" = "#E76F6F",
                                 "PR" = "#7DA7D9"))
  
  if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }
  
  print(p)
}

## =============================================================================
## EXAMPLES
## Heatmap for T-cell targets
plot_average_heatmap(
  data = cmix_t,
  radius_use = 50,
  patient_use = c("P1", "P3", "P4", "P5", "P6"),
  response_use = "PD",
  plot_title = "Average mixing score (50 µm) - T targets - P1, P3, P4, P5, P6 - PD"
)

## Boxplot
create_boxplot(
  data = cmix_mac,
  reference = "Epi",
  target_levels = mac_targets,
  y_label = "Colocalization with Epi",
  patient_use = "P2",
  stat_label_y = 1.5,
  ylim = c(0, 1.5),
  radius_use = 50
)

## compare heatmap
plot_heatmap_by_response(
  data = cmix_t,
  radius_use = 50,
  patient_use = c("P2",
                  "P7"),
  response_levels = c("PD", "PR"),
  plot_title = "T-cell colocalization (within 50μm)"
)

plot_boxplot_by_response(
  data = cmix_t,
  reference = "Epi",
  target_levels = "Cytotoxic_CD8T",
  y_label = "Colocalization with Epi",
  stat_label_y = 1.4,
  ylim = c(0, 1.5),
  radius_use = 50,
  patient_use = c("P1", "P3", "P4", "P5", "P6", "P8", "P9", "P10")
)
