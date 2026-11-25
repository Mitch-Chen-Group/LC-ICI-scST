## =============================================================================
# script adapted from https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/spatial-plotting/
library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(grid)
library(tibble)

## =============================================================================
## Global inputs
base_dir <- "/path/to/folder"
seu_path <- file.path(base_dir, "ALL_Seurat_annot_combined_mac_T.rds")
seu <- readRDS(seu_path)

## =============================================================================
## One row per slide
slide_manifest <- tibble::tribble(
  ~slide_key, ~slide_name, ~raw_data_dir,
  "NR1", "RUN1_NR1", "/path/to/flatFiles_polygons_qced_data/RUN1_NR1",
  "NR2", "RUN1_NR2", "/path/to/flatFiles_polygons_qced_data/RUN1_NR2",
  "R1",  "RUN1_R1", "/path/to/flatFiles_polygons_qced_data/RUN1_R1",
  "R2",  "RUN1_R2", "/path/to/flatFiles_polygons_qced_data/RUN1_R2",
  "NR3", "RUN2_NR3", "/path/to/flatFiles_polygons_qced_data/RUN2_NR3",
  "NR4", "RUN2_NR4", "/path/to/flatFiles_polygons_qced_data/RUN2_NR4",
  "R3",  "RUN2_R3", "/path/to/flatFiles_polygons_qced_data/RUN2_R3"
)
slide_manifest$slide_key <- trimws(as.character(slide_manifest$slide_key))
slide_manifest$slide_name <- trimws(as.character(slide_manifest$slide_name))
slide_manifest$raw_data_dir <- trimws(as.character(slide_manifest$raw_data_dir))

## =============================================================================
## Palettes
dark_cols <- c("Epi"        = "#FF6B6B",
               "Fib"        = "#4CAF50",
               "Mono/Mac"   = "#4FC3F7",
               "Endo"       = "#D4A017",
               "Neu"        = "#81D4FA",
               "T"          = "#BA68C8",
               "B/Plasma"   = "#F48FB1",
               "Unassigned" = "grey40")

light_cols <- c("Epi"        = "#B2182B",
                "Endo"       = "#2166AC",
                "Fib"        = "#FD8D3C",
                "Mono/Mac"   = "#1B9E77",
                "Neu"        = "#6BAED6",
                "T"          = "#984EA3",
                "B/Plasma"   = "#8C8C00",
                "Unassigned" = "grey70")

# dark_cols_combined <- c("Epi"            = "#FF6B6B",
#                         "Endo"           = "#D4A017",
#                         "Fib"            = "#4CAF50",
#                         "Neu"            = "#81D4FA",
#                         "B/Plasma"       = "#F48FB1",
#                         "Cytotoxic_CD8T" = "#9C27B0",
#                         "Activated_CD4T" = "#BA68C8",
#                         "Treg"           = "#6A1B9A",
#                         "Cycling_like_T" = "#CE93D8",
#                         "MMP19_Mac"      = "#29B6F6",
#                         "SPP1_Mac"       = "#FB8C00",
#                         "MARCO_Mac"      = "#00ACC1",
#                         "CXCL9_Mac"      = "#26A69A",
#                         "Unassigned"     = "grey55")

light_cols_combined <- c("Epi"            = "#D73027", 
                         "Endo"           = "#2166AC",
                         "Fib"            = "#4DAF4A",
                         "Neu"            = "#6BAED6",
                         "B/Plasma"       = "#CC79A7",
                         "Cytotoxic_CD8T" = "#542788",
                         "Activated_CD4T" = "#984EA3",
                         "Treg"           = "#542788",
                         "Cycling_like_T" = "#C2A5CF",
                         "MMP19_Mac"      = "yellow",
                         "SPP1_Mac"       = "yellow",
                         "MARCO_Mac"      = "yellow",
                         "CXCL9_Mac"      = "#00FFFF",
                         "Unassigned"     = "grey90")

dark_cols_combined <- c("Epi"            = "#D73027",
                        "Endo"           = "#2166AC",
                        "Fib"            = "#4DAF4A",
                        "Neu"            = "#6BAED6",
                        "B/Plasma"       = "#CC79A7",
                        "Cytotoxic_CD8T" = "#00E5FF",
                        "Activated_CD4T" = "#984EA3",
                        "Treg"           = "#542788",
                        "Cycling_like_T" = "#C2A5CF",
                        "MMP19_Mac"      = "yellow",
                        "SPP1_Mac"       = "yellow",
                        "MARCO_Mac"      = "yellow",
                        "CXCL9_Mac"      = "#00FFFF",
                        "Unassigned"     = "grey90")

## =============================================================================
## Scale bar helper
make_scale_bar <- function(x_vals, y_vals, microns_per_pixel = 0.12028) {
  x_range <- range(x_vals, na.rm = TRUE)
  y_range <- range(y_vals, na.rm = TRUE)
  
  x_length_um <- diff(x_range) * microns_per_pixel
  target <- x_length_um / 4
  
  order <- 10^floor(log10(target))
  mantissa <- target / order
  nice_mantissa <- if (mantissa < 1.5) {
    1
  } else if (mantissa < 3.5) {
    2
  } else if (mantissa < 7.5) {
    5
  } else {
    10
  }
  
  scale_length_um <- nice_mantissa * order
  scale_length_px <- scale_length_um / microns_per_pixel
  
  scale_label <- if (scale_length_um >= 1000) {
    paste0(scale_length_um / 1000, " mm")
  } else {
    paste0(scale_length_um, " \u00B5m")
  }
  
  x_start <- x_range[2] - scale_length_px * 1.1
  x_end   <- x_range[2] - scale_length_px * 0.1
  y_pos   <- y_range[1] + diff(y_range) * 0.04
  
  list(
    bg = annotation_custom(
      grob = rectGrob(gp = gpar(fill = "white", alpha = 0.80, col = NA)),
      xmin = x_start - scale_length_px * 0.08,
      xmax = x_end + scale_length_px * 0.08,
      ymin = y_pos - scale_length_px * 0.08,
      ymax = y_pos + scale_length_px * 0.35
    ),
    rect = annotation_custom(
      grob = rectGrob(gp = gpar(fill = "black", col = NA)),
      xmin = x_start,
      xmax = x_end,
      ymin = y_pos,
      ymax = y_pos + scale_length_px * 0.05
    ),
    label = annotation_custom(
      grob = textGrob(scale_label, gp = gpar(col = "black"), just = "center"),
      xmin = (x_start + x_end) / 2,
      xmax = (x_start + x_end) / 2,
      ymin = y_pos + scale_length_px * 0.14,
      ymax = y_pos + scale_length_px * 0.14
    )
  )
}

## =============================================================================
## Cache environment so raw files are loaded only once per slide
.cosmx_cache <- new.env(parent = emptyenv())

load_cosmx_slide <- function(slide_key, slide_manifest) {
  slide_key <- trimws(as.character(slide_key))
  
  idx <- which(trimws(as.character(slide_manifest$slide_key)) == slide_key)
  
  if (length(idx) != 1) {
    stop(
      "slide_key not found uniquely in slide_manifest: '", slide_key, "'. ",
      "Manifest keys are: ", paste(slide_manifest$slide_key, collapse = ", ")
    )
  }
  
  if (exists(slide_key, envir = .cosmx_cache, inherits = FALSE)) {
    return(get(slide_key, envir = .cosmx_cache, inherits = FALSE))
  }
  
  slide_name   <- slide_manifest$slide_name[idx]
  raw_data_dir <- slide_manifest$raw_data_dir[idx]
  
  polygons_file <- file.path(raw_data_dir, paste0(slide_name, "-polygons.csv.gz"))

  if (!file.exists(polygons_file)) stop("Polygon file not found: ", polygons_file)

  obj <- list(
    cell_polygons = fread(polygons_file)
  )
  
  assign(slide_key, obj, envir = .cosmx_cache)
  obj
}

## =============================================================================
## Helper to find which slide a patient/FOV belongs to
find_patient_fov_slide <- function(seu, patient_id, fov_id, slide_col = "slideID") {
  meta <- seu@meta.data
  
  req_cols <- c("PatientID", "fov", slide_col)
  missing_cols <- setdiff(req_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop("Missing metadata columns in Seurat object: ", paste(missing_cols, collapse = ", "))
  }
  
  hit <- meta %>%
    dplyr::filter(.data$PatientID == patient_id, .data$fov == fov_id) %>%
    dplyr::distinct(.data[[slide_col]])
  
  if (nrow(hit) == 0) {
    stop("No cells found for patient ", patient_id, " and FOV ", fov_id, ".")
  }
  
  if (nrow(hit) > 1) {
    stop("This patient/FOV combination maps to more than one slide. Check metadata.")
  }
  
  trimws(as.character(hit[[slide_col]][[1]]))
}

## =============================================================================
## Main plotting function
plot_cosmx_fov <- function(
    seu,
    patient_id,
    fov_id,
    slide_manifest,
    slide_col = "slideID",
    annot_col = "major_annot",
    exclude_unassigned = TRUE,
    theme_style = c("dark", "light"),
    border_color_dark = "grey80",
    border_color_light = "black",
    border_size = 0.2,
    microns_per_pixel = 0.12028,
    highlight_celltypes = NULL,
    non_highlight_color_dark = "grey20",
    non_highlight_color_light = "grey85",
    subset_highlight = FALSE
) {
  theme_style <- match.arg(theme_style)
  
  meta <- seu@meta.data
  
  req_cols <- c("PatientID", "fov", "cell", annot_col, slide_col)
  missing_cols <- setdiff(req_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop("Missing metadata columns in Seurat object: ", paste(missing_cols, collapse = ", "))
  }
  
  slide_key <- find_patient_fov_slide(
    seu = seu,
    patient_id = patient_id,
    fov_id = fov_id,
    slide_col = slide_col
  )
  
  raw_obj <- load_cosmx_slide(slide_key, slide_manifest)
  cell_polygons <- raw_obj$cell_polygons
  
  meta_sub <- meta %>%
    dplyr::filter(.data$PatientID == patient_id, .data$fov == fov_id)
  
  if (exclude_unassigned) {
    meta_sub <- meta_sub %>%
      dplyr::filter(.data[[annot_col]] != "Unassigned")
  }
  
  if (nrow(meta_sub) == 0) {
    stop("No Seurat cells left after filtering for patient/FOV.")
  }
  
  cells_keep <- as.character(meta_sub$cell)
  
  poly_df <- cell_polygons[
    as.character(cell_polygons[["cell"]]) %in% cells_keep &
      cell_polygons[["fov"]] %in% fov_id,
  ] %>%
    as.data.frame()
  
  if (nrow(poly_df) == 0) {
    stop("No polygon rows found for patient ", patient_id, ", FOV ", fov_id, ", slide ", slide_key, ".")
  }
  
  ## SAFE annotation lookup
  annot_lookup <- setNames(as.character(meta_sub[[annot_col]]), as.character(meta_sub$cell))
  poly_df$cell_type <- unname(annot_lookup[as.character(poly_df$cell)])
  
  ## Anything unmatched becomes Unassigned
  poly_df$cell_type[is.na(poly_df$cell_type)] <- "Unassigned"
  
  ## base palette
  pal <- if (annot_col == "major_annot") {
    if (theme_style == "dark") dark_cols else light_cols
  } else if (annot_col == "annot_combined") {
    if (theme_style == "dark") dark_cols_combined else light_cols_combined
  } else {
    stop("annot_col must be either 'major_annot' or 'annot_combined'.")
  }
  
  ## optional highlighting
  if (!is.null(highlight_celltypes)) {
    highlight_celltypes <- as.character(highlight_celltypes)
    
    missing_highlight <- setdiff(highlight_celltypes, unique(poly_df$cell_type))
    if (length(missing_highlight) > 0) {
      message(
        "These highlighted cell types are not present in this FOV: ",
        paste(missing_highlight, collapse = ", ")
      )
    }
    
    if (subset_highlight) {
      poly_df <- poly_df %>%
        dplyr::filter(.data$cell_type %in% highlight_celltypes)
      
      if (nrow(poly_df) == 0) {
        stop("After subsetting, no highlighted cell types remain in this FOV.")
      }
      
      pal_use <- pal[names(pal) %in% highlight_celltypes]
      poly_df$plot_type <- poly_df$cell_type
    } else {
      non_highlight_color <- if (theme_style == "dark") {
        non_highlight_color_dark
      } else {
        non_highlight_color_light
      }
      
      poly_df$plot_type <- dplyr::case_when(
        poly_df$cell_type %in% highlight_celltypes ~ poly_df$cell_type,
        TRUE ~ "Other"
      )
      
      pal_use <- c(
        pal[names(pal) %in% highlight_celltypes],
        "Other" = non_highlight_color
      )
    }
  } else {
    poly_df$plot_type <- poly_df$cell_type
    pal_use <- pal
  }
  
  poly_df$plot_type <- factor(poly_df$plot_type, levels = names(pal_use))
  
  scale_bar <- make_scale_bar(
    x_vals = poly_df$x_global_px,
    y_vals = poly_df$y_global_px,
    microns_per_pixel = microns_per_pixel
  )
  
  p <- ggplot(
    poly_df,
    aes(x = x_global_px, y = y_global_px, group = cell, fill = plot_type)
  ) +
    geom_polygon(
      color = if (theme_style == "dark") border_color_dark else border_color_light,
      linewidth = border_size
    ) +
    scale_fill_manual(values = pal_use, drop = TRUE, na.translate = FALSE) +
    coord_fixed() +
    theme_void() +
    labs(
      fill = "Cell type",
      title = paste0("Patient ", patient_id, " | FOV ", fov_id, " | Slide ", slide_key)
    ) +
    scale_bar$bg + scale_bar$rect + scale_bar$label
  
  if (theme_style == "dark") {
    p <- p +
      theme(
        plot.background   = element_rect(fill = "#111111", colour = NA),
        panel.background  = element_rect(fill = "#111111", colour = NA),
        legend.background = element_rect(fill = "#111111", colour = NA),
        legend.key        = element_rect(fill = "#111111", colour = NA),
        legend.text       = element_text(color = "white"),
        legend.title      = element_text(color = "white"),
        plot.title        = element_text(color = "white", hjust = 0.5)
      )
  } else {
    p <- p +
      theme(
        plot.background   = element_rect(fill = "white", colour = NA),
        panel.background  = element_rect(fill = "white", colour = NA),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key        = element_rect(fill = "white", colour = NA),
        legend.text       = element_text(color = "black"),
        legend.title      = element_text(color = "black"),
        plot.title        = element_text(color = "black", hjust = 0.5)
      )
  }
  
  return(p)
}

## =============================================================================
## Example plots

p1 <- plot_cosmx_fov(
  seu = seu,
  patient_id = "P1",
  fov_id = 173,
  slide_manifest = slide_manifest,
  slide_col = "slideID",
  annot_col = "annot_combined",
  theme_style = "dark",
  exclude_unassigned = FALSE,
  highlight_celltypes = c("Epi", "Fib", "Cytotoxic_CD8T")
)
p1

p3 <- plot_cosmx_fov(
  seu = seu,
  patient_id = "P3",
  fov_id = 151,
  slide_manifest = slide_manifest,
  slide_col = "slideID",
  annot_col = "annot_combined",
  theme_style = "dark",
  exclude_unassigned = FALSE,
  highlight_celltypes = c("Epi", "Fib", "Cytotoxic_CD8T", "SPP1_Mac")
)
p3