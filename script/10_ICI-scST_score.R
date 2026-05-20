## =====================================================
## Libraries
library(Seurat)
library(dplyr)
library(tidyr)
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
    mac_epi_coloc = ifelse(mac_epi_coloc <= 0, NA, mac_epi_coloc),
    ## Final metrics
    metric_cd8_pdl1 = cd8_fraction * PDL1_fraction,
    metric_cd8_pdl1_coloc = (cd8_fraction * PDL1_fraction) / mac_epi_coloc
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

stat_cd8_pdl1$y.position <-
  max(plot_df$metric_cd8_pdl1,
      na.rm = TRUE) * 1.1

stat_cd8_pdl1_coloc$y.position <-
  max(plot_df$metric_cd8_pdl1_coloc,
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
    y = "CD8 X PD-L1",
    fill = "Response"
  ) +
  theme_classic(base_size = 12)
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
    y = "(PD-L1 × CD8) / macrophage–epithelial colocalization",
    fill = "Response"
  ) +
  theme_classic(base_size = 12)
print(p2)

