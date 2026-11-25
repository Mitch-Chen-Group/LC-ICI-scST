## =====================================================
## Epithelial SPP1 downstream score vs distance to the nearest SPP1+ macrophage, PD vs PR
## =====================================================
## Libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(FNN)
library(lme4)
library(lmerTest)
library(broom.mixed)

set.seed(1)

## =====================================================
## Paths + input
base_dir <- "/path/to/folder"
seu <- readRDS(file.path(base_dir, "ALL_Seurat_annot_combined_mac_T.rds"))

## =====================================================
## SPP1 downstream score (15 genes, AddModuleScore on log-normalised data, x100)
SPP1_downstream <- c(
  "AKT1", "STAT3", "PTK2", "SRC", "FN1", "MYC", "CCND1", "GSK3B", "MTOR",
  "SPP1", "CD44", "ITGAV", "ITGB5", "ITGA5", "ITGB1"
)
SPP1_downstream <- intersect(SPP1_downstream, rownames(seu))

seu$SPP1_rc <- FetchData(seu, vars = "SPP1")[, 1]   # SPP1 as normalised in script 02 (used for the FOV figure)

median_nFeature <- median(seu$nFeature_RNA)
seu <- NormalizeData(seu, normalization.method = "LogNormalize",
                     scale.factor = median_nFeature, verbose = FALSE)

nbin_use <- 10
ctrl_use <- max(5, floor(nrow(seu) / nbin_use / 2))
seu <- AddModuleScore(
  object   = seu,
  features = list(SPP1_downstream),
  name     = "SPP1_DownstreamScore",
  nbin     = nbin_use,
  ctrl     = ctrl_use,
  seed     = 1
)
seu$SPP1_DownstreamScore <- seu$SPP1_DownstreamScore1 * 100
seu$SPP1_expr <- FetchData(seu, vars = "SPP1")[, 1]

px_to_um <- 0.12028
meta_df <- seu@meta.data %>%
  mutate(cell_id = rownames(seu@meta.data)) %>%
  transmute(cell_id, PatientID, response_lbl, major_annot, slideID, fov,
            x = CenterX_global_px * px_to_um,
            y = CenterY_global_px * px_to_um,
            SPP1_DownstreamScore, SPP1_expr, SPP1_rc)
rm(seu); invisible(gc())

## =====================================================
## Distance from each epithelial cell to the nearest SPP1+ macrophage (same slide and FOV)
## SPP1+ macrophage = Mono/Mac cell with SPP1 detected (log-normalised SPP1 > 0)
epi_df <- meta_df %>%
  group_by(slideID, fov) %>%
  group_modify(~{
    src <- .x %>% filter(major_annot == "Mono/Mac", SPP1_expr > 0)
    epi <- .x %>% filter(major_annot == "Epi")
    if (nrow(src) < 3 | nrow(epi) < 3) return(tibble())
    nn <- get.knnx(data = as.matrix(src[, c("x", "y")]),
                   query = as.matrix(epi[, c("x", "y")]), k = 1)
    epi %>% mutate(dist_to_nearest_SPP1posMac = nn$nn.dist[, 1])
  }) %>%
  ungroup()

## =====================================================
## Score vs distance plots (score shifted for display only; statistics use the unshifted score)
SCORE_SHIFT <- ceiling(abs(min(epi_df$SPP1_DownstreamScore, na.rm = TRUE)))
epi_df <- epi_df %>% mutate(SPP1_DownstreamScore_plot = SPP1_DownstreamScore + SCORE_SHIFT)

plot_distance_decay <- function(cell_df, title_txt) {
  ggplot(cell_df, aes(x = dist_to_nearest_SPP1posMac, y = SPP1_DownstreamScore_plot,
                      colour = response_lbl, fill = response_lbl)) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 1.3,
                method.args = list(control = loess.control(surface = "direct"))) +
    scale_colour_manual(values = c("PD" = "#D55E00", "PR" = "#0072B2")) +
    scale_fill_manual(values = c("PD" = "#D55E00", "PR" = "#0072B2")) +
    labs(title = title_txt,
         x = "Distance to nearest SPP1+ Mac (µm)",
         y = "Epi SPP1 downstream score",
         colour = "Response", fill = "Response") +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
}

p_nocap <- plot_distance_decay(epi_df, "SPP1 downstream score vs distance to nearest SPP1+ Mac")
print(p_nocap)

cap_um <- 200   # distance cap for the capped plot
p_cap <- plot_distance_decay(epi_df %>% filter(dist_to_nearest_SPP1posMac <= cap_um),
                             paste0("SPP1 downstream score vs distance (capped ", cap_um, "µm)")) +
  coord_cartesian(xlim = c(0, cap_um))
print(p_cap)

## =====================================================
## Statistics by distance cutoff
## slope for PD alone, PR alone, and PD vs PR (interaction), each with
## random-intercept (RI) and random-slope (RS) mixed models (patient as random effect)
cutoffs <- c(100, 200, 300, 400, Inf)

get_slope_row <- function(model, term_name) {
  tidy(model, effects = "fixed", conf.int = TRUE) %>% filter(term == term_name)
}

results <- list()
for (cutoff in cutoffs) {
  d <- if (is.finite(cutoff)) epi_df %>% filter(dist_to_nearest_SPP1posMac <= cutoff) else epi_df
  cutoff_lbl <- if (is.finite(cutoff)) paste0(cutoff, "um") else "no_cutoff"
  pd_d <- d %>% filter(response_lbl == "PD")
  pr_d <- d %>% filter(response_lbl == "PR")

  m_pd_ri <- lmer(SPP1_DownstreamScore ~ dist_to_nearest_SPP1posMac + (1 | PatientID), data = pd_d)
  m_pd_rs <- lmer(SPP1_DownstreamScore ~ dist_to_nearest_SPP1posMac + (1 + dist_to_nearest_SPP1posMac | PatientID),
                  data = pd_d, control = lmerControl(optimizer = "bobyqa"))
  m_pr_ri <- lmer(SPP1_DownstreamScore ~ dist_to_nearest_SPP1posMac + (1 | PatientID), data = pr_d)
  m_pr_rs <- lmer(SPP1_DownstreamScore ~ dist_to_nearest_SPP1posMac + (1 + dist_to_nearest_SPP1posMac | PatientID),
                  data = pr_d, control = lmerControl(optimizer = "bobyqa"))
  m_int_ri <- lmer(SPP1_DownstreamScore ~ dist_to_nearest_SPP1posMac * response_lbl + (1 | PatientID), data = d)
  m_int_rs <- lmer(SPP1_DownstreamScore ~ dist_to_nearest_SPP1posMac * response_lbl +
                     (1 + dist_to_nearest_SPP1posMac | PatientID),
                   data = d, control = lmerControl(optimizer = "bobyqa"))

  pd_ri <- get_slope_row(m_pd_ri, "dist_to_nearest_SPP1posMac")
  pd_rs <- get_slope_row(m_pd_rs, "dist_to_nearest_SPP1posMac")
  pr_ri <- get_slope_row(m_pr_ri, "dist_to_nearest_SPP1posMac")
  pr_rs <- get_slope_row(m_pr_rs, "dist_to_nearest_SPP1posMac")
  int_ri <- get_slope_row(m_int_ri, "dist_to_nearest_SPP1posMac:response_lblPR")
  int_rs <- get_slope_row(m_int_rs, "dist_to_nearest_SPP1posMac:response_lblPR")

  results[[cutoff_lbl]] <- data.frame(
    cutoff = cutoff_lbl,
    n_total = nrow(d), n_PD = nrow(pd_d), n_PR = nrow(pr_d),
    PD_slope_RI = pd_ri$estimate, PD_p_RI = pd_ri$p.value,
    PD_slope_RS = pd_rs$estimate, PD_p_RS = pd_rs$p.value,
    PR_slope_RI = pr_ri$estimate, PR_p_RI = pr_ri$p.value,
    PR_slope_RS = pr_rs$estimate, PR_p_RS = pr_rs$p.value,
    # interaction = PR slope minus PD slope; sign flipped so positive = PD declines more than PR
    PDvsPR_diff_RI = -int_ri$estimate, PDvsPR_p_RI = int_ri$p.value,
    PDvsPR_diff_RS = -int_rs$estimate, PDvsPR_p_RS = int_rs$p.value
  )
}
cutoff_tbl <- bind_rows(results)
print(cutoff_tbl)
write.csv(cutoff_tbl, file.path(base_dir, "SPP1_AMS15_SPP1posMac_allcutoffs_stats.csv"), row.names = FALSE)

## =====================================================
## Example FOVs: epithelial score (diverging, shifted as in the line plots) and SPP1 in macrophages (green)
gp <- quantile(epi_df$SPP1_DownstreamScore, c(0.01, 0.99), na.rm = TRUE)
epi_limits <- SCORE_SHIFT + c(-max(abs(gp)), max(abs(gp)))
mac_limits <- c(0, 2)   # log1p(SPP1)

plot_fov_spatial <- function(patient_id, fov_id, title_txt) {
  epi_d <- epi_df %>% filter(PatientID == patient_id, fov == fov_id)
  other_all <- meta_df %>% filter(PatientID == patient_id, fov == fov_id)
  mac_d <- other_all %>% filter(major_annot == "Mono/Mac")
  bg_d <- other_all %>% filter(!major_annot %in% c("Mono/Mac", "Epi"))

  ggplot() +
    geom_point(data = bg_d, aes(x, y, shape = "Other cells"),
               colour = "grey85", size = 0.8, alpha = 0.5) +
    geom_point(data = epi_d, aes(x, y, colour = SPP1_DownstreamScore_plot, shape = "Epithelial"),
               size = 1.8, alpha = 0.9) +
    scale_colour_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = SCORE_SHIFT,
      limits = epi_limits, oob = scales::squish, na.value = "grey90",
      name = paste0("Epithelial\nSPP1 downstream\nscore (+", SCORE_SHIFT, " shift)\n(",
                    SCORE_SHIFT, " = signature-null)")
    ) +
    ggnewscale::new_scale_colour() +
    geom_point(data = mac_d, aes(x, y, colour = log1p(SPP1_rc), shape = "Macrophage"),
               size = 2.6, alpha = 0.95) +
    scale_colour_gradient(low = "#C7E9C0", high = "#00441B",
                          limits = mac_limits, oob = scales::squish,
                          name = "Macrophage\nlog1p(SPP1)") +
    scale_shape_manual(name = "Cell type",
                       values = c("Other cells" = 16, "Epithelial" = 16, "Macrophage" = 17)) +
    coord_fixed() +
    labs(title = title_txt, x = expression(paste("x (", mu, "m)")), y = expression(paste("y (", mu, "m)"))) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold"), legend.box = "vertical")
}

p_pd <- plot_fov_spatial("P6", 44, "P6 | FOV 44 (PD)")
p_pr <- plot_fov_spatial("P9", 26, "P9 | FOV 26 (PR)")
print(p_pd)
print(p_pr)
