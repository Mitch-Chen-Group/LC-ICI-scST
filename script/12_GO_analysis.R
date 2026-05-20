## =====================================================
## Libraries
library(dplyr)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

base_dir <- "/path/to/folder"
## =====================================================
## Load macrophage-associated gene list
mac_gene_tbl <- read.csv("significant_genes_macrophage_associated.csv",
                         stringsAsFactors = FALSE)
mac_genes <- mac_gene_tbl %>%
  filter(macrophage_associated == "Yes") %>%
  pull(gene) %>%
  unique()

## =====================================================
## CLEAN gene names + collapse probes
merged <- read.csv(file.path(base_dir,
                             "NRvsR_and_CXCL9_VS_nonCXCL9_mac_DE_merged.csv"),stringsAsFactors = FALSE)
merged_clean <- merged %>%
  mutate(gene_clean = str_replace(gene, "\\.\\d+$", "")) %>%
  group_by(gene_clean) %>%
  slice_max(order_by = abs(log2FC_resp) + abs(log2FC_mac), n = 1) %>%
  ungroup() %>%
  filter(gene_clean %in% mac_genes)

## =====================================================
## Extract gene sets
genes_R_CXCL9 <- merged_clean %>%
  filter(cat == "R_up & CXCL9_Mac_up") %>%
  pull(gene_clean) %>%
  unique()
genes_NR_nonCXCL9 <- merged_clean %>%
  filter(cat == "NR_up & nonCXCL9_Mac_up") %>%
  pull(gene_clean) %>%
  unique()

## =====================================================
## Convert to ENTREZ IDs
conv_R <- bitr(genes_R_CXCL9,
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = org.Hs.eg.db) %>% 
  distinct(ENTREZID, .keep_all = TRUE)

conv_NR <- bitr(genes_NR_nonCXCL9,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db) %>% 
  distinct(ENTREZID, .keep_all = TRUE)

## =====================================================
## GO enrichment
ego_R <- enrichGO(gene = conv_R$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP", pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,readable = TRUE)

ego_NR <- enrichGO(gene = conv_NR$ENTREZID,
                   OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,readable = TRUE)

## =====================================================
## GO plots (clean + consistent)
wrap_fun <- function(x) stringr::str_wrap(x, width = 45)
p_go_R <- dotplot(ego_R, showCategory = 15) +
  scale_y_discrete(labels = wrap_fun) +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  ggtitle("Response-associated Macrophage Signature") +
  theme_bw(base_size = 10) +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 9),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))

p_go_NR <- dotplot(ego_NR, showCategory = 15) +
  scale_y_discrete(labels = wrap_fun) +
  ggtitle("Resistance-associated Macrophage Signature") +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  theme_bw(base_size = 10) +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 9),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))
print(p_go_R)
print(p_go_NR)
