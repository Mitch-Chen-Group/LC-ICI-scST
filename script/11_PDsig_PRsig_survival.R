## =====================================================
## Libraries
library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(Matrix)
library(readxl)
library(survival)
library(survminer)

set.seed(1234)

## =====================================================
## Paths + load data
base_dir <- "/path/to/folder"
setwd(base_dir)
seu <- readRDS(file.path(base_dir, "ALL_Seurat_annot_combined_mac_T.rds"))
mac_only <- readRDS(file.path(base_dir, "ALL_Seurat_Mac_Recluster_2.rds"))

## =====================================================
## DE: PR vs PD
## Positive log2FC_resp = up in Responders
Idents(mac_only) <- mac_only$response_lbl

de_resp <- FindMarkers(
  mac_only,
  ident.1 = "PR",
  ident.2 = "PD",
  logfc.threshold = 0,
  min.pct = 0,
  test.use = "wilcox",
  verbose = FALSE
) %>%
  rownames_to_column("gene") %>%
  transmute(
    gene,
    log2FC_resp = avg_log2FC,
    pval_resp   = p_val,
    padj_resp   = p_val_adj
  )

## =====================================================
## DE:CXCL9_Mac vs all other macrophages
## Positive log2FC_mac = up in CXCL9_Mac
mac_only$mac_vs_rest <- ifelse(
  mac_only$mac_subtype == "CXCL9_Mac",
  "CXCL9_Mac",
  "nonCXCL9_Mac"
)

Idents(mac_only) <- mac_only$mac_vs_rest

de_mac <- FindMarkers(
  mac_only,
  ident.1 = "CXCL9_Mac",
  ident.2 = "nonCXCL9_Mac",
  logfc.threshold = 0,
  min.pct = 0,
  test.use = "wilcox",
  verbose = FALSE
) %>%
  rownames_to_column("gene") %>%
  transmute(
    gene,
    log2FC_mac = avg_log2FC,
    pval_mac   = p_val,
    padj_mac   = p_val_adj
  )

## =====================================================
## Merge DE tables
merged <- full_join(de_resp, de_mac, by = "gene") %>%
  filter(padj_resp < 0.05 | padj_mac < 0.05) %>%
  mutate(
    sig_both = padj_resp < 0.05 & padj_mac < 0.05,
    cat = case_when(
      sig_both & log2FC_resp >  0.5 & log2FC_mac >  0.5 ~ "R_up & CXCL9_Mac_up",
      sig_both & log2FC_resp < -0.5 & log2FC_mac < -0.5 ~ "NR_up & nonCXCL9_Mac_up",
      TRUE ~ "other"
    )
  )

write.csv(merged, "NRvsR_and_CXCL9_VS_nonCXCL9_mac_DE_merged.csv", row.names = F) # inspect DEs

## =====================================================
## Quadrant plot (log2FC_mac vs log2FC_resp)
mac_gene_tbl <- read.csv("significant_genes_macrophage_associated.csv",stringsAsFactors = FALSE) # macrophage-specifc/associated genes only

## keep only macrophage-associated genes
genes_to_label <- mac_gene_tbl %>%
  dplyr::filter(macrophage_associated == "Yes") %>%
  dplyr::pull(gene)

length(genes_to_label)
head(genes_to_label)

lab_df <- merged %>% filter(gene %in% genes_to_label)
cols <- c("R_up & CXCL9_Mac_up" = "#00BFC4",
          "NR_up & nonCXCL9_Mac_up"    = "#F8766D",
          "other"                  = "grey90")

p_quad <- ggplot(merged, aes(x = log2FC_mac, y = log2FC_resp)) +
  geom_point(aes(color = cat), alpha = 0.7, size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  scale_color_manual(values = cols) +
  geom_text_repel(data = lab_df,
                  aes(label = gene),
                  size = 3.0,
                  fontface = "bold",
                  color = "black",
                  max.overlaps = Inf,
                  min.segment.length = 0) +
  labs(x = "avg log2FC (CXCL9_Mac vs nonCXCL9_Mac; + = CXCL9_Mac up)",
       y = "avg log2FC (Responders vs Non-responders; + = Responders up)",
       color = NULL) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top")
print(p_quad)

## =====================================================
## Weighted signatures + PFS (our cohort)
## Score is computed per cell and averaged per patient
expr <- GetAssayData(seu, layer = "data")  # genes x cells
patients <- seu$PatientID

pb_list <- lapply(unique(patients), function(pid) {
  cells <- which(patients == pid)
  Matrix::rowMeans(expr[, cells, drop = FALSE])
})
names(pb_list) <- unique(patients)
pb_expr <- do.call(cbind, pb_list)

gene_weights_NR <- merged %>%
  filter(cat == "NR_up & nonCXCL9_Mac_up", gene %in% genes_to_label) %>%
  mutate(weight = abs(log2FC_resp * log2FC_mac) * -log10(padj_resp * padj_mac + 1e-300)) %>%
  select(gene, weight)

gene_weights_R <- merged %>%
  filter(cat == "R_up & CXCL9_Mac_up", gene %in% genes_to_label) %>%
  mutate(weight = abs(log2FC_resp * log2FC_mac) * -log10(padj_resp * padj_mac + 1e-300)) %>%
  select(gene, weight)

# calculate scores
compute_pb_score <- function(pb_mat, gene_weights) {
  common_genes <- intersect(rownames(pb_mat), gene_weights$gene)
  if(length(common_genes) == 0){
    stop("No overlap between pb genes and gene_weights genes.")
  }
  w <- gene_weights$weight[match(common_genes, gene_weights$gene)]
  scores <- colSums(pb_mat[common_genes, , drop = FALSE] * w)
  tibble(
    PatientID = colnames(pb_mat),
    score = scores
  )
}

score_NR <- compute_pb_score(pb_expr, gene_weights_NR) %>%
  rename(NRup_nonCXCL9up_score = score) # PD_sig

score_R <- compute_pb_score(pb_expr, gene_weights_R) %>%
  rename(Rup_CXCL9up_score = score) # PR_sig

## =====================================================
# our cohort
# PFS
dir <- "/path/to/folder"
pfs <- read_excel(file.path(dir,"Patient_PFS.xlsx"), sheet = 1)
colnames(pfs) <- c("PatientID", "Progression", "PFS_months")
pfs$PatientID <- as.character(pfs$PatientID)

surv_df <- pfs %>%
  left_join(score_NR, by = "PatientID") %>%
  left_join(score_R,  by = "PatientID")

## KM: NR signature (median split)
surv_df$group_NR <- ifelse(surv_df$NRup_nonCXCL9up_score >= median(surv_df$NRup_nonCXCL9up_score, na.rm = TRUE),
                           "High", "Low")
surv_df$group_NR <- factor(surv_df$group_NR, levels = c("Low", "High"))
fit_NR <- survfit(Surv(PFS_months, Progression) ~ group_NR, data = surv_df)
ggsurvplot(
  fit_NR,
  data = surv_df,
  pval = T,                   # log-rank p-value
  conf.int = F,
  break.time.by = 6,          # x-axis breaks every 6 months
  risk.table.y.text.col = T,
  risk.table.y.text = T,
  legend.labs = c("Low", "High"),
  risk.table = T,
  palette = "jco",
  xlab = "Time in Months",
  ylab = "PFS Probability",
  surv.median.line = "v",
  cumevents = F,
  surv.scale = "percent"
)

## KM: R signature (median split)
surv_df$group_R <- ifelse(
  surv_df$Rup_CXCL9up_score >= median(surv_df$Rup_CXCL9up_score, na.rm = TRUE),
  "High", "Low"
)
surv_df$group_R <- factor(surv_df$group_R, levels = c("Low", "High"))
fit_R <- survfit(Surv(PFS_months, Progression) ~ group_R, data = surv_df)
ggsurvplot(
  fit_R,
  data = surv_df,
  pval = T,                   # log-rank p-value
  conf.int = F,
  break.time.by = 6,          # x-axis breaks every 6 months
  risk.table.y.text.col = T,
  risk.table.y.text = T,
  legend.labs = c("Low", "High"),
  risk.table = T,
  palette = "jco",
  xlab = "Time in Months",
  ylab = "PFS Probability",
  surv.median.line = "v",
  cumevents = F,
  surv.scale = "percent"
)

## Cox
cox_fit_NR <- coxph(Surv(PFS_months, Progression) ~ NRup_nonCXCL9up_score, data = surv_df)
print(summary(cox_fit_NR))
cox_fit_R <- coxph(Surv(PFS_months, Progression) ~ Rup_CXCL9up_score, data = surv_df)
print(summary(cox_fit_R))
## Save cohort scores
write.csv(surv_df, file.path(base_dir,"our_cohort_CXCL9_vs_nonCXCL9_signature_scores_with_PFS.csv"), row.names = FALSE)

## =====================================================
# SU2C-MARK COHORT
## =====================================================
## =====================================================
## Bulk RNA-seq validation (SU2C-MARK)
dir <- "/path/to/folder"
## input matrix & metadata
expr_mat <- readRDS(file.path(dir, "SU2C_MARK_log2tpm_gene_mat.rds"))
rna_metadata <- readRDS(file.path(dir, "SU2C_MARK_metadata.rds"))

## =====================================================
## Prepare expression + metadata
bulk_expr <- expr_mat
bulk_meta <- rna_metadata

common_samples <- intersect(colnames(bulk_expr), bulk_meta$SampleID)
bulk_expr <- bulk_expr[, common_samples]

bulk_meta <- bulk_meta %>%
  filter(SampleID %in% common_samples) %>%
  slice(match(common_samples, SampleID))
bulk_meta$Harmonized_Confirmed_BOR <- factor(bulk_meta$Harmonized_Confirmed_BOR,
                                             levels = c("PD", "SD", "PR", "CR"))

bulk_meta <- bulk_meta %>%
  mutate(response_lbl = case_when(Harmonized_Confirmed_BOR %in% c("PR", "CR") ~ "Responder",
                                  TRUE ~ "Non-responder"))

## =====================================================
## Stats helpers
run_all_tests <- function(df, score_col){
  cat("\n=============================\n")
  cat("Score:", score_col, "\n")
  cat("=============================\n")
  
  print(table(df$Harmonized_Confirmed_BOR))
  
  print(kruskal.test(as.formula(paste(score_col, "~ Harmonized_Confirmed_BOR")), data = df))
  
  print(pairwise.wilcox.test(
    df[[score_col]],
    df$Harmonized_Confirmed_BOR,
    p.adjust.method = "BH"
  ))
}

run_binary_test <- function(df, score_col){
  cat("\n=============================\n")
  cat("Binary comparison:", score_col, "\n")
  cat("=============================\n")
  
  print(table(df$response_lbl))
  
  print(wilcox.test(as.formula(paste(score_col, "~ response_lbl")), data = df))
}

## =====================================================
## calculate WEIGHTED score
compute_bulk_score_weighted <- function(expr_mat, gene_weights, score_name){
  common_genes <- intersect(rownames(expr_mat), gene_weights$gene)
  w <- gene_weights$weight[match(common_genes, gene_weights$gene)]
  w[!is.finite(w)] <- 0
  scores <- colSums(expr_mat[common_genes, , drop = FALSE] * w) / sum(abs(w))
  tibble(SampleID = colnames(expr_mat),
         !!score_name := scores)
}

score_NR_w <- compute_bulk_score_weighted(bulk_expr, gene_weights_NR, "NRup_nonCXCL9up_score_w")
score_R_w  <- compute_bulk_score_weighted(bulk_expr, gene_weights_R,  "Rup_CXCL9up_score_w")
df_w <- bulk_meta %>%
  inner_join(inner_join(score_NR_w, score_R_w, by = "SampleID"), by = "SampleID")

## =====================================================
## IO (Pembro + Atezo)
df_w_io <- df_w %>% filter(grepl("Pembrolizumab|Atezolizumab", Agent_PD1))

print(ggplot(df_w_io, aes(response_lbl, NRup_nonCXCL9up_score_w, fill=response_lbl)) +
        geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.15)+theme_bw())
print(ggplot(df_w_io, aes(response_lbl, Rup_CXCL9up_score_w, fill=response_lbl)) +
        geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.15)+theme_bw())

run_binary_test(df_w_io, "NRup_nonCXCL9up_score_w")
run_binary_test(df_w_io, "Rup_CXCL9up_score_w")

print(ggplot(df_w_io, aes(Harmonized_Confirmed_BOR, NRup_nonCXCL9up_score_w, fill=Harmonized_Confirmed_BOR)) +
        geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.15)+theme_bw())
print(ggplot(df_w_io, aes(Harmonized_Confirmed_BOR, Rup_CXCL9up_score_w, fill=Harmonized_Confirmed_BOR)) +
        geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.15)+theme_bw())

run_all_tests(df_w_io, "NRup_nonCXCL9up_score_w")
run_all_tests(df_w_io, "Rup_CXCL9up_score_w")

## =====================================================
## PD-to-PR sig
## =====================================================
## Shared color scheme:
## Non-responder / PD / SD -> red-family shades
## Responder / PR / CR     -> teal-family shades
binary_colors <- c(
  "Non-responder" = "#F8766D",   # base red-family (matches PD/SD family below)
  "Responder"     = "#00BFC4"    # base teal-family (matches PR/CR family below)
)

bor_colors <- c(
  "PD" = "#E07A7A",   # darker red  (more severe end of non-responder family)
  "SD" = "#F0A8A8",   # lighter red (milder end of non-responder family)
  "PR" = "#7FC7C6",   # lighter teal (milder end of responder family)
  "CR" = "#4FA9A8"    # darker teal (deeper end of responder family)
)

plot_NR_R_axis <- function(meta_df){
  df <- meta_df %>%
    mutate(ratio_NR_R = NRup_nonCXCL9up_score_w - Rup_CXCL9up_score_w)
  df$response_lbl <- factor(df$response_lbl,
                            levels = c("Non-responder","Responder"))
  df$Harmonized_Confirmed_BOR <- factor(df$Harmonized_Confirmed_BOR,
                                        levels = c("PD","SD","PR","CR"))
  
  ## Binary
  print(ggplot(df,aes(response_lbl,ratio_NR_R,fill = response_lbl)) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(width = 0.15) +
          scale_fill_manual(values = binary_colors) +
          theme_bw() +
          labs(y = "PD-to-PR score",x = NULL))
  print(wilcox.test(ratio_NR_R ~ response_lbl,
                    data = df))
  
  ## iRECIST
  print(ggplot(df,aes(Harmonized_Confirmed_BOR,ratio_NR_R,fill = Harmonized_Confirmed_BOR)) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(width = 0.15) +
          scale_fill_manual(values = bor_colors) +
          theme_bw() +
          labs(y = "PD-to-PR score",x = NULL))
  
  ## global test
  print(kruskal.test(ratio_NR_R ~ Harmonized_Confirmed_BOR,data = df))
  
  ## pairwise comparisons
  pw <- pairwise.wilcox.test(x = df$ratio_NR_R,
                             g = df$Harmonized_Confirmed_BOR,
                             p.adjust.method = "BH")
  print(pw)
}

plot_NR_R_axis(df_w_io)
df_w_io <- df_w_io %>%
  mutate(ratio_NR_R = NRup_nonCXCL9up_score_w -Rup_CXCL9up_score_w)

## =====================================================
## survival analysis
##Prepare dataframe
surv_all <- df_w_io %>%
  select(SampleID,
         Harmonized_PFS_Event, Harmonized_PFS_Days,
         Harmonized_OS_Event, Harmonized_OS_Days,
         NRup_nonCXCL9up_score_w, Rup_CXCL9up_score_w, ratio_NR_R) %>%
  mutate(PFS_event = as.numeric(Harmonized_PFS_Event),
         PFS_time  = as.numeric(Harmonized_PFS_Days) / 30.44,
         OS_event  = as.numeric(Harmonized_OS_Event),
         OS_time   = as.numeric(Harmonized_OS_Days) / 30.44)
df <- surv_all %>%
  filter(!is.na(PFS_time),
         !is.na(PFS_event),
         !is.na(OS_time),
         !is.na(OS_event),
         !is.na(ratio_NR_R))

## =====================================================
## optimal cutoff
minprop_val <- 0.2
## PFS ===
cut_pfs <- surv_cutpoint(df,time = "PFS_time",event = "PFS_event",
                         variables = "ratio_NR_R",minprop = minprop_val)
cut_val_pfs <- cut_pfs$cutpoint$cutpoint
df$ratio_NR_R_opt_pfs <- factor(
  ifelse(df$ratio_NR_R > cut_val_pfs, "High", "Low"),
  levels = c("Low","High"))
## KM
fit_pfs_opt <- survfit(Surv(PFS_time, PFS_event) ~ ratio_NR_R_opt_pfs, data = df)
ggsurvplot(
  fit_pfs_opt,
  data = df,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("blue","red"),
  title = "PFS (optimal cutoff)")
## Cox
print(summary(coxph(Surv(PFS_time, PFS_event) ~ ratio_NR_R_opt_pfs, data = df)))

## OS ===
cut_os <- surv_cutpoint(df,
                        time = "OS_time",
                        event = "OS_event",
                        variables = "ratio_NR_R",
                        minprop = minprop_val)
cut_val_os <- cut_os$cutpoint$cutpoint
df$ratio_NR_R_opt_os <- factor(ifelse(df$ratio_NR_R > cut_val_os, "High", "Low"),
                               levels = c("Low","High"))
## KM
fit_os_opt <- survfit(Surv(OS_time, OS_event) ~ ratio_NR_R_opt_os, data = df)
ggsurvplot(
  fit_os_opt,
  data = df,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("blue","red"),
  title = "OS (optimal cutoff)"
)
## Cox
print(summary(coxph(Surv(OS_time, OS_event) ~ ratio_NR_R_opt_os, data = df)))

## =====================================================
## Administrative censoring at 20 months

cutoff <- 20

## PFS
df$PFS_time20  <- pmin(df$PFS_time, cutoff)
df$PFS_event20 <- ifelse(df$PFS_time <= cutoff & df$PFS_event == 1, 1, 0)

## OS
df$OS_time20  <- pmin(df$OS_time, cutoff)
df$OS_event20 <- ifelse(df$OS_time <= cutoff & df$OS_event == 1, 1, 0)


## =====================================================
## Optimal cutoff
minprop_val <- 0.15

## PFS ===
cut_pfs <- surv_cutpoint(
  df,
  time = "PFS_time20",
  event = "PFS_event20",
  variables = "ratio_NR_R",
  minprop = minprop_val
)

cut_val_pfs <- cut_pfs$cutpoint$cutpoint

df$ratio_NR_R_opt_pfs <- factor(
  ifelse(df$ratio_NR_R > cut_val_pfs, "High", "Low"),
  levels = c("Low", "High")
)

## KM
fit_pfs_opt <- survfit(
  Surv(PFS_time20, PFS_event20) ~ ratio_NR_R_opt_pfs,
  data = df
)

ggsurvplot(
  fit_pfs_opt,
  data = df,
  pval = F,
  risk.table = TRUE,
  palette = c("#1F4E99", "#FF1F1F"),
  legend.title = "",
  legend = "none",
  legend.labs = c("Low", "High"),
  xlab = "Months",
  title = ""
)

## Cox
print(summary(
  coxph(Surv(PFS_time20, PFS_event20) ~ ratio_NR_R_opt_pfs, data = df)
))

print(summary(
  coxph(Surv(PFS_time20, PFS_event20) ~ ratio_NR_R, data = df)
))



## =====================================================
## OS ===
cut_os <- surv_cutpoint(
  df,
  time = "OS_time20",
  event = "OS_event20",
  variables = "ratio_NR_R",
  minprop = minprop_val
)

cut_val_os <- cut_os$cutpoint$cutpoint

df$ratio_NR_R_opt_os <- factor(
  ifelse(df$ratio_NR_R > cut_val_os, "High", "Low"),
  levels = c("Low", "High")
)

## KM
fit_os_opt <- survfit(
  Surv(OS_time20, OS_event20) ~ ratio_NR_R_opt_os,
  data = df
)

ggsurvplot(
  fit_os_opt,
  data = df,
  pval = F,
  risk.table = TRUE,
  palette = c("#1F4E99", "#FF1F1F"),
  legend.title = "",
  legend = "none",
  legend.labs = c("Low", "High"),
  xlab = "Months",
  title = ""
)

## Cox
print(summary(
  coxph(Surv(OS_time20, OS_event20) ~ ratio_NR_R_opt_os, data = df)
))

print(summary(
  coxph(Surv(OS_time20, OS_event20) ~ ratio_NR_R, data = df)
))
