## ===================================================== # all 10 samples
## Libraries
library(Seurat)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(ggplot2)
library(scales)
library(readxl)
library(survival)
library(survminer)
library(ggrepel)

## =====================================================
## Paths + input
dir <- "/path/to/folder"
mac_only <- readRDS(file.path(dir, "ALL_Seurat_Mac_Recluster_2.rds"))

## =====================================================
## metadata
meta <- mac_only@meta.data %>%
  as_tibble(rownames = "cell_id") %>%
  transmute(cell_id,
            PatientID = PatientID,
            response  = response,
            response_lbl = response_lbl,
            mac_subtype = mac_subtype)

## =====================================================
## Composition
## By response (stacked bar, within macrophages)
mac_subtype_prop_by_response <- meta %>%
  dplyr::count(response_lbl, mac_subtype, name = "n") %>%
  dplyr::group_by(response_lbl) %>%
  dplyr::mutate(prop_within_mac = n / sum(n)) %>%
  dplyr::ungroup()

p_comp_response <- ggplot(mac_subtype_prop_by_response,
                          aes(x = response_lbl, y = prop_within_mac, fill = mac_subtype)) +
  geom_col() +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "Response", y = "Macrophage subtype proportion",
       fill = "Mac subtype") +
  theme_classic()

print(p_comp_response)

## By patient (stacked bar, within macrophages)
subtype_prop_by_patient_plot <- meta %>%
  dplyr::count(PatientID, response_lbl, mac_subtype, name = "n") %>%
  dplyr::group_by(PatientID, response_lbl) %>%
  dplyr::mutate(prop_within_patient = n / sum(n)) %>%
  dplyr::ungroup()

p_comp_patient <- ggplot(subtype_prop_by_patient_plot,
                         aes(x = PatientID, y = prop_within_patient, fill = mac_subtype)) +
  geom_col() +
  facet_wrap(~ response_lbl, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "Patient", y = "Macrophage subtype proportion") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_comp_patient)

## Boxplots: per-patient subtype proportion within macrophages, grouped by PD/PR
subtype_prop_by_patient <- meta %>%
  dplyr::group_by(PatientID, response_lbl, mac_subtype) %>%
  dplyr::summarise(n = n(), .groups = "drop_last") %>%
  dplyr::mutate(total_mac = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(prop = n / total_mac)

p_box <- ggplot(subtype_prop_by_patient,
                aes(x = mac_subtype, y = prop, fill = response_lbl)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               width = 0.7, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = response_lbl),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 2, alpha = 0.7) +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "Macrophage subtype", y = "Proportion within macrophages", fill = "Group") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_box)

## Wilcoxon per subtype (table)
stats_results <- subtype_prop_by_patient %>%
  dplyr::group_by(mac_subtype) %>%
  dplyr::summarise(
    p.value = wilcox.test(prop ~ response_lbl)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(p.adj = p.adjust(p.value, method = "BH"))

print(stats_results)

## =====================================================
# composition without P7
meta_noP7 <- meta %>%
  dplyr::filter(PatientID != "P7")

mac_subtype_prop_by_response_xP7 <- meta_noP7 %>%
  dplyr::count(response_lbl, mac_subtype, name = "n") %>%
  dplyr::group_by(response_lbl) %>%
  dplyr::mutate(prop_within_mac = n / sum(n)) %>%
  dplyr::ungroup()

p_comp_response_xP7 <- ggplot(mac_subtype_prop_by_response_xP7,
                              aes(x = response_lbl, y = prop_within_mac, fill = mac_subtype)) +
  geom_col() +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "Response", y = "Macrophage subtype proportion",
       fill = "Mac subtype") +
  theme_classic()

print(p_comp_response_xP7)

subtype_prop_by_patient_plot_xP7 <- meta_noP7 %>%
  dplyr::count(PatientID, response_lbl, mac_subtype, name = "n") %>%
  dplyr::group_by(PatientID, response_lbl) %>%
  dplyr::mutate(prop_within_patient = n / sum(n)) %>%
  dplyr::ungroup()

p_comp_patient_xP7 <- ggplot(subtype_prop_by_patient_plot_xP7,
                             aes(x = PatientID, y = prop_within_patient, fill = mac_subtype)) +
  geom_col() +
  facet_wrap(~ response_lbl, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "Patient", y = "Macrophage subtype proportion") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_comp_patient_xP7)

subtype_prop_by_patient_xP7 <- meta_noP7 %>%
  dplyr::group_by(PatientID, response_lbl, mac_subtype) %>%
  dplyr::summarise(n = n(), .groups = "drop_last") %>%
  dplyr::mutate(total_mac = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(prop = n / total_mac)

p_box_xP7 <- ggplot(subtype_prop_by_patient_xP7,
                    aes(x = mac_subtype, y = prop, fill = response_lbl)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               width = 0.7, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = response_lbl),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 2, alpha = 0.7) +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "Macrophage subtype", y = "Proportion within macrophages", fill = "Group") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_box_xP7)

subtype_prop_by_patient_xP7 <- meta_noP7 %>%
  dplyr::group_by(PatientID, response_lbl, mac_subtype) %>%
  dplyr::summarise(n = n(), .groups = "drop_last") %>%
  dplyr::mutate(total_mac = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(prop = n / total_mac)

p_box_xP7 <- ggplot(subtype_prop_by_patient_xP7,
                    aes(x = mac_subtype, y = prop, fill = response_lbl)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               width = 0.7, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = response_lbl),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 2, alpha = 0.7) +
  geom_text(aes(label = PatientID),
            position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.95),
            size = 2.1, vjust = -0.7, color = "black", show.legend = FALSE) +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "Macrophage subtype", y = "Proportion within macrophages", fill = "Group") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_box_xP7)

stats_results_xP7 <- subtype_prop_by_patient_xP7 %>%
  dplyr::group_by(mac_subtype) %>%
  dplyr::summarise(
    p.value = wilcox.test(prop ~ response_lbl)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(p.adj = p.adjust(p.value, method = "BH"))

print(stats_results_xP7)


## =====================================================
## Ratio
props_wide <- subtype_prop_by_patient %>%
  # filter(mac_subtype %in% c("CXCL9_Mac","MMP19_Mac","SPP1_Mac","MARCO_Mac")) %>%
  select(PatientID, response_lbl, mac_subtype, prop) %>%
  pivot_wider(names_from = mac_subtype, values_from = prop, values_fill = 0)

mac_ratios <- props_wide %>%
  mutate(
    CXCL9_over_MMP19 = `CXCL9_Mac` / `MMP19_Mac`,
    CXCL9_over_SPP1 = `CXCL9_Mac` / `SPP1_Mac`,
    CXCL9_over_MARCO = `CXCL9_Mac` / `MARCO_Mac`
  ) %>%
  select(PatientID, response_lbl, CXCL9_over_MMP19, CXCL9_over_SPP1, CXCL9_over_MARCO)


mac_ratios_long <- mac_ratios %>%
  pivot_longer(cols = starts_with("CXCL9_over_"),
               names_to = "ratio_type", values_to = "ratio") %>%
  mutate(
    ratio_type = recode(ratio_type,
                        CXCL9_over_MMP19 = "CXCL9_Mac / MMP19_Mac",
                        CXCL9_over_SPP1 = "CXCL9_Mac / SPP1_Mac",
                        CXCL9_over_MARCO = "CXCL9_Mac / MARCO_Mac"),
    ratio_type = factor(ratio_type, levels = c("CXCL9_Mac / MMP19_Mac", "CXCL9_Mac / SPP1_Mac", "CXCL9_Mac / MARCO_Mac"))
  )

# per-facet Wilcoxon
ratio_stats <- mac_ratios_long %>%
  dplyr::group_by(ratio_type) %>%
  dplyr::summarise(p = wilcox.test(ratio ~ response_lbl)$p.value, .groups = "drop") %>%
  dplyr::mutate(label = paste0("p = ", signif(p, 3)))

ann_df <- mac_ratios_long %>%
  dplyr::group_by(ratio_type) %>%
  dplyr::summarise(y = max(ratio, na.rm = TRUE) * 1.05, .groups = "drop") %>%
  dplyr::left_join(ratio_stats, by = "ratio_type") %>%
  dplyr::mutate(x = 1.5)

p_ratio <- ggplot(mac_ratios_long, aes(x = response_lbl, y = ratio, fill = response_lbl)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.75, coef=Inf) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  facet_wrap(~ ratio_type, scales = "free_y") +
  geom_text(data = ann_df, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, vjust = 0, size = 3.6) +
  labs(x = "Response", y = "Ratio", title = "Macrophage ratios per patient") +
  theme_classic() +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6))

print(p_ratio)

# with patient label
p_ratio <- ggplot(mac_ratios_long, aes(x = response_lbl, y = ratio, fill = response_lbl)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.75, coef = Inf) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  geom_text_repel(
    aes(label = PatientID),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    show.legend = FALSE) +
  facet_wrap(~ ratio_type, scales = "free_y") +
  geom_text(
    data = ann_df, aes(x = x, y = y, label = label),
    inherit.aes = FALSE, vjust = 0, size = 3.6) +
  labs(x = "Response", y = "Ratio", title = "Macrophage ratios per patient") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6))

print(p_ratio)

## =====================================================
## Ratio (FLIPPED: CXCL9_Mac is denominator)
props_wide <- subtype_prop_by_patient %>%
  # filter(mac_subtype %in% c("CXCL9_Mac","MMP19_Mac","SPP1_Mac","MARCO_Mac")) %>%
  select(PatientID, response_lbl, mac_subtype, prop) %>%
  pivot_wider(names_from = mac_subtype, values_from = prop, values_fill = 0)

mac_ratios <- props_wide %>%
  mutate(
    MMP19_over_CXCL9 = `MMP19_Mac` / `CXCL9_Mac`,
    SPP1_over_CXCL9  = `SPP1_Mac` / `CXCL9_Mac`,
    MARCO_over_CXCL9 = `MARCO_Mac` / `CXCL9_Mac`
  ) %>%
  select(PatientID, response_lbl, MMP19_over_CXCL9, SPP1_over_CXCL9, MARCO_over_CXCL9)

mac_ratios_long <- mac_ratios %>%
  pivot_longer(cols = c(MMP19_over_CXCL9, SPP1_over_CXCL9, MARCO_over_CXCL9),
               names_to = "ratio_type", values_to = "ratio") %>%
  mutate(
    ratio_type = recode(ratio_type,
                        MMP19_over_CXCL9  = "MMP19_Mac / CXCL9_Mac",
                        SPP1_over_CXCL9  = "SPP1_Mac / CXCL9_Mac",
                        MARCO_over_CXCL9 = "MARCO_Mac / CXCL9_Mac"),
    ratio_type = factor(ratio_type,
                        levels = c("MMP19_Mac / CXCL9_Mac", 
                                   "SPP1_Mac / CXCL9_Mac",
                                   "MARCO_Mac / CXCL9_Mac"))
  )

# per-facet Wilcoxon
ratio_stats <- mac_ratios_long %>%
  dplyr::group_by(ratio_type) %>%
  dplyr::summarise(p = wilcox.test(ratio ~ response_lbl)$p.value, .groups = "drop") %>%
  dplyr::mutate(label = paste0("p = ", signif(p, 3)))

ann_df <- mac_ratios_long %>%
  dplyr::group_by(ratio_type) %>%
  dplyr::summarise(y = max(ratio, na.rm = TRUE) * 1.05, .groups = "drop") %>%
  dplyr::left_join(ratio_stats, by = "ratio_type") %>%
  dplyr::mutate(x = 1.5)

p_ratio <- ggplot(mac_ratios_long, aes(x = response_lbl, y = ratio, fill = response_lbl)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.75, coef = Inf) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  facet_wrap(~ ratio_type, scales = "free_y") +
  geom_text(data = ann_df, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, vjust = 0, size = 3.6) +
  labs(x = "Response", y = "Ratio", title = "Macrophage ratios per patient (CXCL9 as denominator)") +
  theme_classic() +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6))

print(p_ratio)

# with patient label
p_ratio <- ggplot(mac_ratios_long,
                  aes(x = response_lbl, y = ratio, fill = response_lbl)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.75, coef = Inf) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  geom_text_repel(
    aes(label = PatientID),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    show.legend = FALSE) +
  facet_wrap(~ ratio_type, scales = "free_y") +
  geom_text(
    data = ann_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    vjust = 0,
    size = 3.6) +
  labs(
    x = "Response",
    y = "Ratio",
    title = "Macrophage ratios per patient (CXCL9 as denominator)") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6))

print(p_ratio)

## =====================================================
## Non-CXCL9 / CXCL9 ratio (including P7)

props_wide <- subtype_prop_by_patient %>%
  select(PatientID, response_lbl, mac_subtype, prop) %>%
  pivot_wider(names_from = mac_subtype, values_from = prop, values_fill = 0)

mac_ratios_extra <- props_wide %>%
  mutate(
    nonCXCL9 = MMP19_Mac + SPP1_Mac + MARCO_Mac,
    ratio = nonCXCL9 / CXCL9_Mac
  ) %>%
  select(PatientID, response_lbl, ratio)

ratio_stats_extra <- mac_ratios_extra %>%
  summarise(p = wilcox.test(ratio ~ response_lbl)$p.value) %>%
  mutate(label = paste0("p = ", signif(p, 3)))

ann_df_extra <- mac_ratios_extra %>%
  summarise(y = max(ratio, na.rm = TRUE) * 1.05) %>%
  bind_cols(ratio_stats_extra) %>%
  mutate(x = 1.5)

p_ratio_extra <- ggplot(mac_ratios_extra,
                        aes(x = response_lbl, y = ratio, fill = response_lbl)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.75, coef = Inf) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  geom_text_repel(
    aes(label = PatientID),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    show.legend = FALSE
  ) +
  geom_text(
    data = ann_df_extra,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    vjust = 0,
    size = 3.6
  ) +
  labs(
    x = "Response",
    y = "Non-CXCL9 / CXCL9 ratio",
    title = "Macrophage balance per patient"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

print(p_ratio_extra)

## =====================================================
## All ratios together in one faceted plot
props_wide <- subtype_prop_by_patient %>%
  select(PatientID, response_lbl, mac_subtype, prop) %>%
  pivot_wider(names_from = mac_subtype, values_from = prop, values_fill = 0)

all_ratios_long <- props_wide %>%
  mutate(
    MMP19_over_CXCL9   = MMP19_Mac / CXCL9_Mac,
    SPP1_over_CXCL9    = SPP1_Mac  / CXCL9_Mac,
    MARCO_over_CXCL9   = MARCO_Mac / CXCL9_Mac,
    nonCXCL9_over_CXCL9 = (MMP19_Mac + SPP1_Mac + MARCO_Mac) / CXCL9_Mac
  ) %>%
  select(
    PatientID, response_lbl,
    MMP19_over_CXCL9, SPP1_over_CXCL9,
    MARCO_over_CXCL9, nonCXCL9_over_CXCL9
  ) %>%
  pivot_longer(
    cols = c(
      MMP19_over_CXCL9,
      SPP1_over_CXCL9,
      MARCO_over_CXCL9,
      nonCXCL9_over_CXCL9
    ),
    names_to = "ratio_type",
    values_to = "ratio"
  ) %>%
  mutate(
    ratio_type = recode(
      ratio_type,
      MMP19_over_CXCL9    = "MMP19_Mac / CXCL9_Mac",
      SPP1_over_CXCL9     = "SPP1_Mac / CXCL9_Mac",
      MARCO_over_CXCL9    = "MARCO_Mac / CXCL9_Mac",
      nonCXCL9_over_CXCL9 = "NonCXCL9_Mac / CXCL9_Mac"
    ),
    ratio_type = factor(
      ratio_type,
      levels = c(
        "MMP19_Mac / CXCL9_Mac",
        "SPP1_Mac / CXCL9_Mac",
        "MARCO_Mac / CXCL9_Mac",
        "NonCXCL9_Mac / CXCL9_Mac"
      )
    )
  )

## per-facet Wilcoxon
ratio_stats_all <- all_ratios_long %>%
  dplyr::group_by(ratio_type) %>%
  dplyr::summarise(
    p = wilcox.test(ratio ~ response_lbl)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(label = paste0("p = ", signif(p, 3)))

## annotation y-position per facet
ann_df_all <- all_ratios_long %>%
  dplyr::group_by(ratio_type) %>%
  dplyr::summarise(
    y = max(ratio, na.rm = TRUE) * 1.05,
    .groups = "drop"
  ) %>%
  dplyr::left_join(ratio_stats_all, by = "ratio_type") %>%
  dplyr::mutate(x = 1.5)

## plot
p_ratio_all <- ggplot(
  all_ratios_long,
  aes(x = response_lbl, y = ratio, fill = response_lbl)
) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.75, coef = Inf) +
  geom_jitter(width = 0.12, size = 1.5, alpha = 0.8) +
  geom_text_repel(
    aes(label = PatientID),
    size = 2.5,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    show.legend = FALSE
  ) +
  facet_wrap(~ ratio_type, scales = "free_y", nrow=1) +
  geom_text(
    data = ann_df_all,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    vjust = 0,
    size = 3.6
  ) +
  labs(
    x = "Response",
    y = "Ratio",
    title = "Macrophage ratios per patient"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  ) 

print(p_ratio_all)

## =====================================================
## Exclude P7 first
subtype_prop_by_patient_xP7 <- subtype_prop_by_patient %>%
  dplyr::filter(PatientID != "P7")

## =====================================================
## Ratio (p7 excluded)
props_wide <- subtype_prop_by_patient_xP7 %>%
  # filter(mac_subtype %in% c("CXCL9_Mac","MMP19_Mac","SPP1_Mac","MARCO_Mac")) %>%
  select(PatientID, response_lbl, mac_subtype, prop) %>%
  pivot_wider(names_from = mac_subtype, values_from = prop, values_fill = 0)

mac_ratios <- props_wide %>%
  mutate(
    CXCL9_over_MMP19 = `CXCL9_Mac` / `MMP19_Mac`,
    CXCL9_over_SPP1  = `CXCL9_Mac` / `SPP1_Mac`,
    CXCL9_over_MARCO = `CXCL9_Mac` / `MARCO_Mac`
  ) %>%
  select(PatientID, response_lbl, CXCL9_over_MMP19, CXCL9_over_SPP1, CXCL9_over_MARCO)

mac_ratios_long <- mac_ratios %>%
  pivot_longer(cols = starts_with("CXCL9_over_"),
               names_to = "ratio_type", values_to = "ratio") %>%
  mutate(
    ratio_type = recode(ratio_type,
                        CXCL9_over_MMP19 = "CXCL9_Mac / MMP19_Mac",
                        CXCL9_over_SPP1  = "CXCL9_Mac / SPP1_Mac",
                        CXCL9_over_MARCO = "CXCL9_Mac / MARCO_Mac"),
    ratio_type = factor(ratio_type,
                        levels = c("CXCL9_Mac / MMP19_Mac",
                                   "CXCL9_Mac / SPP1_Mac",
                                   "CXCL9_Mac / MARCO_Mac"))
  )

# per-facet Wilcoxon
ratio_stats <- mac_ratios_long %>%
  dplyr::group_by(ratio_type) %>%
  dplyr::summarise(
    p = wilcox.test(ratio ~ response_lbl)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(label = paste0("p = ", signif(p, 3)))

ann_df <- mac_ratios_long %>%
  dplyr::group_by(ratio_type) %>%
  dplyr::summarise(
    y = max(ratio, na.rm = TRUE) * 1.05,
    .groups = "drop"
  ) %>%
  dplyr::left_join(ratio_stats, by = "ratio_type") %>%
  dplyr::mutate(x = 1.5)

p_ratio <- ggplot(mac_ratios_long, aes(x = response_lbl, y = ratio, fill = response_lbl)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.75, coef = Inf) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  facet_wrap(~ ratio_type, scales = "free_y") +
  geom_text(data = ann_df, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, vjust = 0, size = 3.6) +
  labs(x = "Response", y = "Ratio", title = "Macrophage ratios per patient (P7 excluded)") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

print(p_ratio)

# with patient label
p_ratio <- ggplot(mac_ratios_long, aes(x = response_lbl, y = ratio, fill = response_lbl)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.75, coef = Inf) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  geom_text_repel(
    aes(label = PatientID),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    show.legend = FALSE
  ) +
  facet_wrap(~ ratio_type, scales = "free_y") +
  geom_text(
    data = ann_df, aes(x = x, y = y, label = label),
    inherit.aes = FALSE, vjust = 0, size = 3.6
  ) +
  labs(x = "Response", y = "Ratio", title = "Macrophage ratios per patient (P7 excluded)") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

print(p_ratio)


## =====================================================
## Ratio (FLIPPED: CXCL9_Mac is denominator) (p7 excluded)
props_wide <- subtype_prop_by_patient_xP7 %>%
  # filter(mac_subtype %in% c("CXCL9_Mac","MMP19_Mac","SPP1_Mac","MARCO_Mac")) %>%
  select(PatientID, response_lbl, mac_subtype, prop) %>%
  pivot_wider(names_from = mac_subtype, values_from = prop, values_fill = 0)

mac_ratios <- props_wide %>%
  mutate(
    MMP19_over_CXCL9 = `MMP19_Mac` / `CXCL9_Mac`,
    SPP1_over_CXCL9  = `SPP1_Mac` / `CXCL9_Mac`,
    MARCO_over_CXCL9 = `MARCO_Mac` / `CXCL9_Mac`
  ) %>%
  select(PatientID, response_lbl, MMP19_over_CXCL9, SPP1_over_CXCL9, MARCO_over_CXCL9)

mac_ratios_long <- mac_ratios %>%
  pivot_longer(cols = c(MMP19_over_CXCL9, SPP1_over_CXCL9, MARCO_over_CXCL9),
               names_to = "ratio_type", values_to = "ratio") %>%
  mutate(
    ratio_type = recode(ratio_type,
                        MMP19_over_CXCL9 = "MMP19_Mac / CXCL9_Mac",
                        SPP1_over_CXCL9  = "SPP1_Mac / CXCL9_Mac",
                        MARCO_over_CXCL9 = "MARCO_Mac / CXCL9_Mac"),
    ratio_type = factor(ratio_type,
                        levels = c("MMP19_Mac / CXCL9_Mac",
                                   "SPP1_Mac / CXCL9_Mac",
                                   "MARCO_Mac / CXCL9_Mac"))
  )

# per-facet Wilcoxon
ratio_stats <- mac_ratios_long %>%
  dplyr::group_by(ratio_type) %>%
  dplyr::summarise(
    p = wilcox.test(ratio ~ response_lbl)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(label = paste0("p = ", signif(p, 3)))

ann_df <- mac_ratios_long %>%
  dplyr::group_by(ratio_type) %>%
  dplyr::summarise(
    y = max(ratio, na.rm = TRUE) * 1.05,
    .groups = "drop"
  ) %>%
  dplyr::left_join(ratio_stats, by = "ratio_type") %>%
  dplyr::mutate(x = 1.5)

p_ratio <- ggplot(mac_ratios_long, aes(x = response_lbl, y = ratio, fill = response_lbl)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.75, coef = Inf) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  facet_wrap(~ ratio_type, scales = "free_y") +
  geom_text(data = ann_df, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, vjust = 0, size = 3.6) +
  labs(x = "Response", y = "Ratio", title = "Macrophage ratios per patient (CXCL9 as denominator, P7 excluded)") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

print(p_ratio)

# with patient label
p_ratio <- ggplot(mac_ratios_long,
                  aes(x = response_lbl, y = ratio, fill = response_lbl)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.75, coef = Inf) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  geom_text_repel(
    aes(label = PatientID),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    show.legend = FALSE
  ) +
  facet_wrap(~ ratio_type, scales = "free_y") +
  geom_text(
    data = ann_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    vjust = 0,
    size = 3.6
  ) +
  labs(
    x = "Response",
    y = "Ratio",
    title = "Macrophage ratios per patient (CXCL9 as denominator, P7 excluded)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

print(p_ratio)

## =====================================================
## Non-CXCL9 / CXCL9 ratio (P7 excluded)
props_wide <- subtype_prop_by_patient_xP7 %>%
  select(PatientID, response_lbl, mac_subtype, prop) %>%
  pivot_wider(names_from = mac_subtype, values_from = prop, values_fill = 0)

mac_ratios_extra <- props_wide %>%
  mutate(
    nonCXCL9 = MMP19_Mac + SPP1_Mac + MARCO_Mac,
    ratio = nonCXCL9 / CXCL9_Mac
  ) %>%
  select(PatientID, response_lbl, ratio)

## Wilcoxon test
ratio_stats_extra <- mac_ratios_extra %>%
  summarise(p = wilcox.test(ratio ~ response_lbl)$p.value) %>%
  mutate(label = paste0("p = ", signif(p, 3)))

## Annotation position
ann_df_extra <- mac_ratios_extra %>%
  summarise(y = max(ratio, na.rm = TRUE) * 1.05) %>%
  bind_cols(ratio_stats_extra) %>%
  mutate(x = 1.5)

p_ratio_extra <- ggplot(mac_ratios_extra,
                        aes(x = response_lbl, y = ratio, fill = response_lbl)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.75, coef = Inf) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  geom_text_repel(
    aes(label = PatientID),
    size = 2.5,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.5,
    show.legend = FALSE
  ) +
  geom_text(
    data = ann_df_extra,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    vjust = 0,
    size = 3.6
  ) +
  labs(
    x = "Response",
    y = "Non-CXCL9 / CXCL9 ratio",
    title = "Macrophage balance per patient (P7 excluded)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

print(p_ratio_extra)

## =====================================================
## All ratios together in one faceted plot (P7 excluded)
props_wide <- subtype_prop_by_patient_xP7 %>%
  select(PatientID, response_lbl, mac_subtype, prop) %>%
  pivot_wider(names_from = mac_subtype, values_from = prop, values_fill = 0)

all_ratios_long <- props_wide %>%
  mutate(
    MMP19_over_CXCL9    = MMP19_Mac / CXCL9_Mac,
    SPP1_over_CXCL9     = SPP1_Mac / CXCL9_Mac,
    MARCO_over_CXCL9    = MARCO_Mac / CXCL9_Mac,
    nonCXCL9_over_CXCL9 = (MMP19_Mac + SPP1_Mac + MARCO_Mac) / CXCL9_Mac
  ) %>%
  select(
    PatientID, response_lbl,
    MMP19_over_CXCL9, SPP1_over_CXCL9,
    MARCO_over_CXCL9, nonCXCL9_over_CXCL9
  ) %>%
  pivot_longer(
    cols = c(
      MMP19_over_CXCL9,
      SPP1_over_CXCL9,
      MARCO_over_CXCL9,
      nonCXCL9_over_CXCL9
    ),
    names_to = "ratio_type",
    values_to = "ratio"
  ) %>%
  mutate(
    ratio_type = recode(
      ratio_type,
      MMP19_over_CXCL9    = "MMP19_Mac / CXCL9_Mac",
      SPP1_over_CXCL9     = "SPP1_Mac / CXCL9_Mac",
      MARCO_over_CXCL9    = "MARCO_Mac / CXCL9_Mac",
      nonCXCL9_over_CXCL9 = "NonCXCL9_Mac / CXCL9_Mac"
    ),
    ratio_type = factor(
      ratio_type,
      levels = c(
        "MMP19_Mac / CXCL9_Mac",
        "SPP1_Mac / CXCL9_Mac",
        "MARCO_Mac / CXCL9_Mac",
        "NonCXCL9_Mac / CXCL9_Mac"
      )
    )
  )

## per-facet Wilcoxon
ratio_stats_all <- all_ratios_long %>%
  dplyr::group_by(ratio_type) %>%
  dplyr::summarise(
    p = wilcox.test(ratio ~ response_lbl)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(label = paste0("p = ", signif(p, 3)))

## annotation y-position per facet
ann_df_all <- all_ratios_long %>%
  dplyr::group_by(ratio_type) %>%
  dplyr::summarise(
    y = max(ratio, na.rm = TRUE) * 1.05,
    .groups = "drop"
  ) %>%
  dplyr::left_join(ratio_stats_all, by = "ratio_type") %>%
  dplyr::mutate(x = 1.5)

## plot
p_ratio_all_xP7 <- ggplot(
  all_ratios_long,
  aes(x = response_lbl, y = ratio, fill = response_lbl)
) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.75, coef = Inf) +
  geom_jitter(width = 0.12, size = 1.5, alpha = 0.8) +
  geom_text_repel(
    aes(label = PatientID),
    size = 2.5,
    max.overlaps = Inf,
    box.padding = 0.1,
    point.padding = 0.1,
    show.legend = FALSE
  ) +
  facet_wrap(~ ratio_type, scales = "free_y", nrow = 1) +
  geom_text(
    data = ann_df_all,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    vjust = 0,
    size = 3.6
  ) +
  labs(
    x = "Response",
    y = "Ratio",
    title = "Macrophage ratios per patient (P7 excluded)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

print(p_ratio_all_xP7)


## =====================================================
## Put macrophage subtypes back into the original Seurat
## (overwrite only the reclustered macrophage cells)
seu <- readRDS(file.path(dir, "ALL_Seurat_major_annot.rds"))

annot_combined <- setNames(as.character(seu$major_annot), colnames(seu))
mac_named <- setNames(as.character(mac_only$mac_subtype), colnames(mac_only))
annot_combined[names(mac_named)] <- mac_named

# Force NA / blanks to "Unassigned" --> these are 'Unassigned' cells from Recluster_1
annot_combined[annot_combined == "Mono/Mac"] <- "Unassigned"

lvl <- c(
  "Epi","Endo","Fib","Neu","T","B/Plasma",
  "MMP19_Mac","SPP1_Mac","MARCO_Mac","CXCL9_Mac",
  "Unassigned"
)
seu$annot_combined <- factor(annot_combined, levels = lvl)
unique(seu$annot_combined) 

# save
saveRDS(seu, file = file.path(dir, "ALL_Seurat_annot_combined.rds"))
