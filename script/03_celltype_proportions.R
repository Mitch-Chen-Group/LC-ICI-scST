## =====================================================
## Libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(rstatix)

## =====================================================
## Paths + input
dir <- "/path/to/folder"
seu <- readRDS(file.path(dir, "ALL_Seurat_major_annot.rds"))

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

## =====================================================
# cell number / proportion per cluster
group_col <- "response"
stopifnot(group_col %in% colnames(seu@meta.data))

df2 <- seu@meta.data %>%
  as.data.frame() %>%
  count(major_annot, !!rlang::sym(group_col), name = "cells") %>%
  group_by(major_annot) %>%
  mutate(prop = cells / sum(cells)) %>%
  ungroup() %>%
  group_by(major_annot) %>%
  summarize(total = sum(cells), .groups = "drop") %>%
  right_join(seu@meta.data %>% as.data.frame() %>%
               count(major_annot, !!rlang::sym(group_col), name = "cells"),
             by = c("major_annot")) %>%
  group_by(major_annot) %>%
  mutate(prop = cells / sum(cells)) %>%
  ungroup()

# order clusters by total size
order_df <- df2 %>% group_by(major_annot) %>% summarize(total = sum(cells), .groups="drop") %>%
  arrange(desc(total))
df2$major_annot <- factor(df2$major_annot, levels = order_df$major_annot)

# stacked counts
p_stack_counts <- ggplot(df2, aes(x = major_annot, y = cells, fill = .data[[group_col]])) +
  geom_col() +
  coord_flip() +
  labs(x = "Cluster", y = "Cell count", fill = group_col,
       title = "Cell numbers per cluster (stacked)")

# stacked proportions (each bar sums to 1)
p_stack_props <- ggplot(df2, aes(x = major_annot, y = prop, fill = .data[[group_col]])) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Cluster", y = "Proportion within cluster", fill = group_col,
       title = "Proportion per cluster (stacked)")

p_stack_counts
p_stack_props

## =====================================================
# cell number / proportion per patient
# check required columns
unique(seu$PatientID)
stopifnot(all(c("PatientID", "major_annot") %in% colnames(seu@meta.data)))

df <- seu@meta.data %>%
  as.data.frame() %>%
  count(PatientID, major_annot, name = "cells") %>%
  group_by(PatientID) %>%
  mutate(prop = cells / sum(cells)) %>%
  ungroup()

# stacked barplot: proportion of each cell type per patient
p <- ggplot(df, aes(x = PatientID, y = prop, fill = major_annot)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = light_cluster_colors) +
  labs(x = "Patient", y = "Proportion", fill = "Cell type",
       title = "Cell type proportions per patient") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p

# stacked barplot: proportion of each cell type per patient (without 'Unassigned')
df_no_unassigned <- df %>%
  filter(major_annot != "Unassigned")

p_no_unassigned <- ggplot(df_no_unassigned,
                          aes(x = PatientID, y = prop, fill = major_annot)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = light_cluster_colors) +
  labs(x = "Patient", y = "Proportion", fill = "Cell type",
       title = "Cell type proportions per patient (excluding Unassigned)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_no_unassigned

## =====================================================
# cell number / proportion per response
# check required columns
stopifnot(all(c("response_lbl", "major_annot") %in% colnames(seu@meta.data)))

df <- seu@meta.data %>%
  as.data.frame() %>%
  count(response_lbl, major_annot, name = "cells") %>%
  group_by(response_lbl) %>%
  mutate(prop = cells / sum(cells)) %>%
  ungroup()

# stacked barplot: proportion of each cell type per response group
p <- ggplot(df, aes(x = response_lbl, y = prop, fill = major_annot)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = light_cluster_colors) +
  labs(x = "Response group", y = "Proportion", fill = "Cell type",
       title = "Cell type proportions per response group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p

# stacked barplot: proportion of each cell type per response group (without 'Unassigned')
df_no_unassigned <- df %>%
  filter(major_annot != "Unassigned")

p_no_unassigned <- ggplot(df_no_unassigned,
                          aes(x = response_lbl, y = prop, fill = major_annot)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = light_cluster_colors) +
  labs(x = "Response group", y = "Proportion", fill = "Cell type",
       title = "Cell type proportions per response group (excluding Unassigned)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_no_unassigned

## =====================================================
## Boxplot of proportion per sample, split by response (PD/PR)
stopifnot(all(c("PatientID", "major_annot", "response_lbl") %in% colnames(seu@meta.data)))

# proportion per sample (PatientID) for each major_annot
df_prop <- seu@meta.data %>%
  as.data.frame() %>%
  filter(major_annot != "Unassigned") %>%   # remove Unassigned
  count(PatientID, response_lbl, major_annot, name = "cells") %>%
  group_by(PatientID) %>%
  mutate(prop = cells / sum(cells)) %>%
  ungroup()
y_pos <- df_prop %>%
  group_by(major_annot) %>%
  summarize(y = max(prop) * 1.05)
stats_df <- df_prop %>%
  group_by(major_annot) %>%
  wilcox_test(prop ~ response_lbl) %>%     # PD vs PR
  adjust_pvalue(method = "BH") %>%         # FDR across cell types
  add_significance("p.adj") %>%
  ungroup()
stats_df <- stats_df %>%
  left_join(y_pos, by = "major_annot")

# plot: boxplot + points
p_box <- ggplot(df_prop, aes(x = major_annot, y = prop, fill = response_lbl)) +
  geom_boxplot(position = position_dodge(width = 0.8),
    width = 0.65,
    outlier.shape = NA,
    coef = Inf,
    alpha = 0.85) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
             size = 2, shape = 21, stroke = 0.8) +
  labs(x = NULL, y = "Proportion per sample", fill = "Group") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank()) +
  stat_pvalue_manual(stats_df,
                     label = "p.adj.signif",
                     y.position = "y",
                     tip.length = 0.01,
                     hide.ns = TRUE)
p_box

## =====================================================
## Table: cell number + percentage per PatientID (major_annot level)

stopifnot(all(c("PatientID", "major_annot") %in% colnames(seu@meta.data)))

tab_long <- seu@meta.data %>%
  as.data.frame() %>%
  count(PatientID, major_annot, name = "n_cells") %>%
  group_by(PatientID) %>%
  mutate(
    pct = 100 * n_cells / sum(n_cells),
    pct = round(pct, 2)
  ) %>%
  ungroup() %>%
  arrange(PatientID, desc(n_cells))

tab_wide <- tab_long %>%
  tidyr::pivot_wider(
    names_from = major_annot,
    values_from = c(n_cells, pct),
    values_fill = 0,
    names_glue = "{major_annot}_{.value}"
  )

tab_long <- as.data.frame(tab_long)
tab_wide <- as.data.frame(tab_wide)
write.csv(tab_wide, file.path(dir, "major_annot_count_percentage_PatientID.csv"))

## exclude Unassigned
tab_long_no_unassigned <- tab_long %>% filter(major_annot != "Unassigned")

tab_wide_no_unassigned <- tab_long_no_unassigned %>%
  group_by(PatientID) %>%
  mutate(
    pct = 100 * n_cells / sum(n_cells),
    pct = round(pct, 2)
  ) %>%
  ungroup() %>%
  tidyr::pivot_wider(
    names_from = major_annot,
    values_from = c(n_cells, pct),
    values_fill = 0,
    names_glue = "{major_annot}_{.value}"
  )

tab_wide_no_unassigned <- as.data.frame(tab_wide_no_unassigned)
write.csv(tab_wide_no_unassigned, file.path(dir, "major_annot_count_percentage_xUnassigned_PatientID.csv"))

