## ===================================================== # all 10 samples
## Libraries
library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)
library(scales)

## =====================================================
## Paths + input
dir <- "/path/to/folder"
T_only <- readRDS(file.path(dir, "ALL_Seurat_T_Recluster_2.rds"))

## =====================================================
## metadata
meta <- T_only@meta.data %>%
  as_tibble(rownames = "cell_id") %>%
  transmute(cell_id,
            PatientID = PatientID,
            response  = response,
            response_lbl = response_lbl,
            T_subtype = T_subtype)

palette.pbmc <- c("#b30024",
                  "#08519c",
                  "#f28ac6",
                  "#20b2aa")

## =====================================================
## Composition
## By response (stacked bar, within T cells)
T_subtype_prop_by_response <- meta %>%
  dplyr::count(response_lbl, T_subtype, name = "n") %>%
  dplyr::group_by(response_lbl) %>%
  dplyr::mutate(prop_within_T = n / sum(n)) %>%
  dplyr::ungroup()

p_comp_response <- ggplot(T_subtype_prop_by_response,
                          aes(x = response_lbl, y = prop_within_T, fill = T_subtype)) +
  geom_col() +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = palette.pbmc) +
  labs(x = "Response", y = "T subtype proportion",
       fill = "T subtype") +
  theme_classic()

print(p_comp_response)

## By patient (stacked bar, within T cells)
subtype_prop_by_patient_plot <- meta %>%
  dplyr::count(PatientID, response_lbl, T_subtype, name = "n") %>%
  dplyr::group_by(PatientID, response_lbl) %>%
  dplyr::mutate(prop_within_patient = n / sum(n)) %>%
  dplyr::ungroup()

p_comp_patient <- ggplot(subtype_prop_by_patient_plot,
                         aes(x = PatientID, y = prop_within_patient, fill = T_subtype)) +
  geom_col() +
  facet_grid(~ response_lbl, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = palette.pbmc) +
  labs(x = "Patient", y = "T subtype proportion") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_comp_patient)

## Boxplots: per-patient subtype proportion within T cells, grouped by PD/PR
subtype_prop_by_patient <- meta %>%
  dplyr::group_by(PatientID, response_lbl, T_subtype) %>%
  dplyr::summarise(n = n(), .groups = "drop_last") %>%
  dplyr::mutate(total_T = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(prop = n / total_T)

# no patient label
p_box <- ggplot(subtype_prop_by_patient,
                aes(x = T_subtype, y = prop, fill = response_lbl)) +
  geom_boxplot(position = position_dodge(width = 0.8), coef = Inf,
               width = 0.7, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = response_lbl),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 2, alpha = 0.7) +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "T subtype", y = "Proportion within T", fill = "Group") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_box)

# with patient label
p_box <- ggplot(subtype_prop_by_patient,
                aes(x = T_subtype, y = prop, fill = response_lbl)) +
  geom_boxplot(position = position_dodge(width = 0.8), coef = Inf,
               width = 0.7, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = response_lbl),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 2, alpha = 0.7) +
  geom_text(aes(label = PatientID),
            position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.95),
            size = 2.1, vjust = -0.7, color = "black", show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "T subtype", y = "Proportion within T", fill = "Group") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_box)

## Wilcoxon per subtype (table)
stats_results <- subtype_prop_by_patient %>%
  dplyr::group_by(T_subtype) %>%
  dplyr::summarise(
    p.value = wilcox.test(prop ~ response_lbl)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(p.adj = p.adjust(p.value, method = "BH"))

print(stats_results)

median_prop_by_response <- subtype_prop_by_patient %>%
  dplyr::group_by(response_lbl, T_subtype) %>%
  dplyr::summarise(
    median_prop = median(prop, na.rm = TRUE),
    .groups = "drop"
  )

print(median_prop_by_response)

## =====================================================
## Put T subtypes into mac-reclustered seurat
## (overwrite only the reclustered T cells)
seu <- readRDS(file.path(dir, "ALL_Seurat_annot_combined.rds"))

annot_combined <- setNames(as.character(seu$annot_combined), colnames(seu))
T_named <- setNames(as.character(T_only$T_subtype), colnames(T_only))
annot_combined[names(T_named)] <- T_named

# remaining T cells = 'Unassigned' in Recluster_1 -> "Unassigned"
annot_combined[annot_combined == "T"] <- "Unassigned"

lvl <- c(
  "Epi","Endo","Fib","Neu","B/Plasma",
  "Cytotoxic_CD8T", "Activated_CD4T", "Treg", "Cycling_like_T",
  "MMP19_Mac","SPP1_Mac","MARCO_Mac","CXCL9_Mac",
  "Unassigned"
)
seu$annot_combined <- factor(annot_combined, levels = lvl)
unique(seu$annot_combined) 

# save
saveRDS(seu, file = file.path(dir, "ALL_Seurat_annot_combined_mac_T.rds"))
