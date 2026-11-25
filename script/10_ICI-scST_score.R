## =====================================================
## Libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(rstatix)
library(ggpubr)

## =====================================================
## Paths
base_dir <- "/path/to/folder"
seu <- readRDS(file.path(base_dir, "ALL_Seurat_annot_combined_mac_T.rds"))
cmix_mac <- read.csv(file.path(base_dir, "ALL_Mac_targets_colocalization_score.csv"),
                     stringsAsFactors = FALSE)

## =====================================================
## Define cell types
non_cxcl9_mac <- c("MMP19_Mac", "SPP1_Mac", "MARCO_Mac")
cd8_label <- "Cytotoxic_CD8T"

## =====================================================
## PD-L1 table
pdl1_df <- tibble::tribble(
  ~PatientID, ~PDL1_percent,
  "P1",  5,
  "P2",  100,
  "P3",  0,
  "P4",  0,
  "P5",  1,
  "P6",  5,
  "P7",  15,
  "P8",  90,
  "P9",  65,
  "P10", 90
) %>%
  mutate(PDL1_fraction = PDL1_percent / 100)

## =====================================================
## Extract metadata
meta_df <- seu@meta.data %>%
  as.data.frame() %>%
  dplyr::select(PatientID, response_lbl, annot_combined) %>%
  filter(!is.na(PatientID), !is.na(response_lbl), !is.na(annot_combined))

## =====================================================
## Total cells per patient
total_cells <- meta_df %>%
  count(PatientID, name = "total_cells")

## =====================================================
## CD8 counts + fraction
cd8_df <- meta_df %>%
  filter(annot_combined == cd8_label) %>%
  count(PatientID, response_lbl, name = "cd8_count") %>%
  right_join(total_cells, by = "PatientID") %>%
  mutate(cd8_count = ifelse(is.na(cd8_count), 0, cd8_count),
         cd8_fraction = cd8_count / total_cells)

## =====================================================
## Aggregate macrophage–epithelial colocalization
coloc_patient <- cmix_mac %>%
  filter(radius == 50,
         Reference == "Epi",
         Target %in% non_cxcl9_mac) %>%
  group_by(PatientID) %>%
  summarise(mac_epi_coloc = mean(Normalised_mixing_score, na.rm = TRUE),
            .groups = "drop")

## =====================================================
## Merge all data
plot_df <- cd8_df %>%
  left_join(coloc_patient, by = "PatientID") %>%
  left_join(pdl1_df, by = "PatientID") %>%
  mutate(
    cd8_percent = cd8_fraction * 100,
    mac_epi_coloc = ifelse(mac_epi_coloc <= 0, NA, mac_epi_coloc),
    ## Final metrics
    metric_cd8_pdl1 = cd8_percent * PDL1_percent,
    metric_cd8_pdl1_coloc = (cd8_percent * PDL1_percent) / mac_epi_coloc
  ) %>%
  filter(!is.na(metric_cd8_pdl1), !is.na(metric_cd8_pdl1_coloc))
print(plot_df)

## =========================================
#STATISTICS 
stat_cd8_pdl1 <- plot_df %>%
  wilcox_test(metric_cd8_pdl1 ~ response_lbl) %>%
  add_significance("p")

stat_cd8_pdl1_coloc <- plot_df %>%
  wilcox_test(metric_cd8_pdl1_coloc ~ response_lbl) %>%
  add_significance("p")

stat_coloc <- plot_df %>%
  wilcox_test(mac_epi_coloc ~ response_lbl) %>%
  add_significance("p")

stat_cd8_pdl1$y.position <-
  max(plot_df$metric_cd8_pdl1,
      na.rm = TRUE) * 1.1

stat_cd8_pdl1_coloc$y.position <-
  max(plot_df$metric_cd8_pdl1_coloc,
      na.rm = TRUE) * 1.1

stat_coloc$y.position <-
  max(plot_df$mac_epi_coloc,
      na.rm = TRUE) * 1.1

## ===============================================
# plots
## Plot 1
p1 <- ggplot(plot_df,
             aes(x = response_lbl, y = metric_cd8_pdl1, fill = response_lbl)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.7, coef=Inf) +
  geom_jitter(width = 0.12, size = 3) +
  ggrepel::geom_text_repel(aes(label = PatientID), size = 3) +
  stat_pvalue_manual(
    stat_cd8_pdl1,
    label = "p",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position"
  ) +
  labs(
    x = NULL,
    y = "PD-L1 X CD8T",
    fill = "Response"
  ) +
  theme_classic(base_size = 18)
print(p1)

## Plot 2: ICI-scST score
p2 <- ggplot(plot_df,
             aes(x = response_lbl, y = metric_cd8_pdl1_coloc, fill = response_lbl)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.7, coef=Inf) +
  geom_jitter(width = 0.12, size = 3) +
  ggrepel::geom_text_repel(aes(label = PatientID), size = 3) +
  stat_pvalue_manual(
    stat_cd8_pdl1_coloc,
    label = "p",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position"
  ) +
  labs(
    x = NULL,
    y = "ICI-scST Score",
    fill = "Response"
  ) +
  theme_classic(base_size = 18)
print(p2)

## Plot 3: coloc score only
p3 <- ggplot(plot_df,
             aes(x = response_lbl, y = mac_epi_coloc, fill = response_lbl)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.7, coef=Inf) +
  geom_jitter(width = 0.12, size = 3) +
  ggrepel::geom_text_repel(aes(label = PatientID), size = 3) +
  stat_pvalue_manual(
    stat_coloc,
    label = "p",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position"
  ) +
  labs(
    x = NULL,
    y = "macrophage–epithelial colocalization (NMS)",
    fill = "Response"
  ) +
  theme_classic(base_size = 18)
print(p3)

# ===============================================================================
# ROC curves (smoothed, density method)
library(pROC)

## Make sure the response is binary
plot_df$response_binary <- ifelse(plot_df$response_lbl == "PR", 1, 0)

## Helper: build empirical ROC (for threshold), smooth it (for AUC/curve), and plot
plot_smooth_roc <- function(predictor, label, col_line = "blue",
                            top_mar = 3.5,          # <- room for 2-line title
                            title_line = 1.5,       # <- gap between title block and box
                            y_pad = 0.04,          # <- small y-axis padding (top & bottom)
                            x_pad = 0.04,
                            text_x_buffer = -0.08,  # <- gap between AUC text and right edge
                            text_y = 0.03,         # <- vertical position of AUC/threshold text
                            text_cex = 1.2) {
  
  ## Empirical ROC — used to get the Youden threshold on the original predictor scale
  roc_emp <- roc(
    response  = plot_df$response_binary,
    predictor = predictor,
    levels    = c(0, 1),
    direction = "<",
    quiet     = TRUE
  )
  
  ## Smoothed ROC — used for the curve itself and its AUC/CI
  roc_sm <- smooth(roc_emp, method = "density")
  
  auc_val <- auc(roc_sm)
  ci_val  <- ci.auc(roc_sm, method = "bootstrap", boot.n = 2000)
  
  best <- coords(
    roc_emp,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity")
  )
  
  old_par <- par(mar = c(5.1, 4.1, top_mar, 2.1))
  on.exit(par(old_par))
  
  plot(
    roc_sm,
    col  = col_line,
    lwd  = 3,
    main = "",
    xlim = c(1 + x_pad, 0 - x_pad),
    ylim = c(0 - y_pad, 1 + y_pad),   # small, controlled padding instead of none/default
    xaxs = "i",
    yaxs = "i"
  )
  
  ## Draw each title line separately, stacked upward, so they never dip into the box
  main_lines <- strsplit(paste0("ROC Curve: ", label, "\nin Predicting ICI Outcome"), "\n")[[1]]
  n <- length(main_lines)
  for (i in seq_along(main_lines)) {
    mtext(main_lines[i], side = 3, line = title_line + (n - i), cex = 1.2, font = 2)
  }
  
  label_text <- sprintf(
    "AUC = %.3f (95%% CI: %.3f-%.3f)\nOptimal threshold = %.3f",
    auc_val, ci_val[1], ci_val[3], best["threshold"]
  )
  
  ## Right-justify at a fixed small offset from x = 0 (the right edge, since axis is reversed)
  ## adj = c(1, 0) anchors the text's RIGHT edge at text_x_buffer, so it always fits
  ## regardless of how long the numbers turn out to be
  text(
    text_x_buffer, text_y,
    labels = label_text,
    cex = text_cex,
    adj = c(1, 0)
  )
  
  invisible(list(roc_emp = roc_emp, roc_sm = roc_sm,
                 auc = auc_val, ci = ci_val, best = best))
}

## ROC 1: ICI-scST Score
res_coloc <- plot_smooth_roc(
  predictor = plot_df$metric_cd8_pdl1_coloc,
  label     = "ICI-scST Score"
)

## ROC 2: PDL1 X CD8T
res_cd8pdl1 <- plot_smooth_roc(
  predictor = plot_df$metric_cd8_pdl1,
  label     = "PD-L1 x CD8T"
)
