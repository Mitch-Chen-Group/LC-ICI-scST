# LC-ICI-scST

End-to-end analysis pipeline for CosMx spatial single-cell transcriptomics readout for precision immune checkpoint inhibition (ICI) in non-small cell lung cancer.

For the manuscript "Spatial single-cell analysis reveals macrophage microenvironment associated with immunotherapy resistance in non-small cell lung cancer". The manuscript is in press at NPJ Precision Oncology.

## Cohort

10 patients with NSCLC treated with ICI, profiled with the CosMx 1,000-plex RNA panel on 7 slides.

| Slide    | Patients | Response |
|----------|----------|----------|
| RUN1_NR1 | P1, P2   | PD       |
| RUN1_NR2 | P3, P4   | PD       |
| RUN2_NR3 | P5       | PD       |
| RUN2_NR4 | P6       | PD       |
| RUN1_R1  | P7, P8   | PR       |
| RUN1_R2  | P9       | PR       |
| RUN2_R3  | P10      | PR       |

## Repository structure

```         
data/
  slide_level_pre-qc_data/        per slide: expression matrix + cell metadata after FOV selection, before QC (input to script 01)
  patient_level_qced_data/        per patient: the 328,283 cells that passed QC, with final cell-type labels
  flatFiles_polygons_qced_data/   per slide: cell outlines of the QC-passed cells (input to script 09)
others/
  Patient_PFS.xlsx                progression-free survival of the 10 patients (script 11)
  SU2C_MARK_log2tpm_gene_mat.rds  SU2C-MARK bulk RNA-seq, log2(TPM + 1), validation cohort (script 11)
  SU2C_MARK_metadata.rds          SU2C-MARK clinical metadata (script 11)
script/                           analysis scripts 01-14, run in order
```

Large data files are stored with Git LFS. Install Git LFS before cloning, or run `git lfs pull` after cloning.

### Patient-level data (`data/patient_level_qced_data/P*/`)

- `P*_exprs.rds`: genes x cells sparse matrix (962 genes). Values are relative-count normalised (each cell scaled to the median number of detected genes per cell, 60), not raw counts. Raw counts can be recovered as `value * nCount_RNA / 60`.
- `P*_metadata.rds`: one row per cell. Key columns:
  - `PatientID`, `slideID`, `fov`, `global_fovID`, `response_lbl` (PD/PR)
  - `CenterX_global_px`, `CenterY_global_px`, `CenterX_local_px`, `CenterY_local_px` (1 px = 0.12028 µm). Global coordinates are per slide, so they are only comparable within a slide.
  - `scDblFinder_class`, `scDblFinder_score`: doublet calls used in QC
  - `major_annot`: 8 major cell types (Epi, Endo, Fib, Mono/Mac, Neu, T, B/Plasma, Unassigned)
  - `annot_combined`: final labels used in the manuscript, including the macrophage and T-cell subtypes

## Scripts

Each script sets its file locations at the top (`/path/to/folder`). Edit these before running.

| Script | What it does | Main input | Main output |
|------------------|------------------|------------------|------------------|
| `01_QC.R` | combine slides, FOV-border filter, cell QC, doublet removal, background-gene removal | `data/slide_level_pre-qc_data` | `ALL_Seurat_after_qc.rds` |
| `02_normalization_banksy_clustering_annotation.R` | normalisation, BANKSY, Harmony, clustering, major cell types | `ALL_Seurat_after_qc.rds` | `ALL_Seurat_major_annot.rds` |
| `03_celltype_proportions.R` | cell-type proportions per patient and response group | `ALL_Seurat_major_annot.rds` | plots, tables |
| `04_macrophage_recluster.R` | macrophage re-clustering into 4 subtypes | `ALL_Seurat_major_annot.rds` | `ALL_Seurat_Mac_Recluster_2.rds` |
| `05_T_recluster.R` | T-cell re-clustering into 4 subtypes | `ALL_Seurat_major_annot.rds` | `ALL_Seurat_T_Recluster_2.rds` |
| `06_macrophage_composition_ratio.R` | adds macrophage subtypes to the main object; macrophage composition and ratios (CXCL9_Mac as denominator) | scripts 02, 04 | `ALL_Seurat_annot_combined.rds` |
| `07_T_composition.R` | adds T-cell subtypes; T-cell composition | scripts 05, 06 | `ALL_Seurat_annot_combined_mac_T.rds` (main annotated object) |
| `08_macrophage_T_colocalization.R` | macrophage and T-cell colocalisation (SPIAT mixing score, 50 µm) | main object | `ALL_Mac_targets_colocalization_score.csv`, `ALL_T_targets_colocalization_score.csv` |
| `09_plot_fov_polygons.R` | cell-outline plots of selected FOVs | main object, `data/flatFiles_polygons_qced_data` | plots |
| `10_ICI-scST_score.R` | ICI-scST score ((PD-L1 x CD8) / macrophage-epithelial colocalisation), PD vs PR, ROC | main object, script 08 output | plots |
| `11_PDsig_PRsig_survival.R` | macrophage PD and PR signatures; survival in this cohort and in SU2C-MARK | main object, script 04 output, `others/` | signature scores, survival plots |
| `12_GO_analysis.R` | GO enrichment of the signature genes | script 11 output | plots |
| `13_CellChat.R` | spatial CellChat, PR vs PD (P7 excluded) | script 06 output | CellChat objects, plots |
| `14_SPP1_downstream_distance.R` | epithelial SPP1 downstream score vs distance to the nearest SPP1+ macrophage, PD vs PR (mixed models) | main object | plots, `SPP1_AMS15_SPP1posMac_allcutoffs_stats.csv` |

Scripts 11 and 12 also need `significant_genes_macrophage_associated.csv` (curated list of macrophage-associated genes), which is available in the Supplementary Materials of the manuscript.

## Software

Analyses were run in R. Versions of the main packages used to check these scripts:

| Package           | Version            | Package         | Version    |
|-------------------|--------------------|-----------------|------------|
| R                 | 4.6.1              | CellChat        | 2.2.0.9001 |
| Seurat            | 5.6.0              | SPIAT           | 1.14.0     |
| SeuratObject      | 5.4.0              | clusterProfiler | 4.20.3     |
| SeuratWrappers    | GitHub             | org.Hs.eg.db    | 3.23.1     |
| Banksy            | 1.8.1              | survival        | 3.8.6      |
| harmony           | 1.2.4              | survminer       | 0.5.2      |
| scDblFinder       | 1.26.9             | lme4            | 2.0.6      |
| SpatialExperiment | 1.22.0             | lmerTest        | 3.2.1      |
| pROC              | 1.19.1             | ggplot2         | 4.0.3      |