#####
# Extended Data Figure 1D
# 
# 23/7/2026
####
library(tidyverse)
library(ggpubr)
library(SummarizedExperiment)
library(limma)
library(FragPipeAnalystR)
library(ggrepel) 

load("FragPipe Output/Imputed_SE.RData")
FPRime <- RData
FPRime <- test_limma(FPRime, type = "all")
FPRime <- add_rejections(FPRime, alpha = 0.01, lfc = 1)


growth_map <- c(
  "PrimaryCD4"         = "Suspension",
  "Jurkat"             = "Suspension",
  "PrimaryBreastEpis"  = "Adherent",
  "MCF7"               = "Adherent",
  "231"                = "Adherent",
  "MCF10A"             = "Adherent",
  "BB7"                = "PDX",
  "BB3RC31"            = "PDX",
  "HBC34"              = "PDX",
  "KMBC2"              = "Adherent",
  "PrimaryNHU"         = "Adherent"
)

meta_df <- as.data.frame(colData(FPRime))
meta_df$Prefix <- sub("_.*", "", meta_df$label)
meta_df$Prefix <- sub(" $", "", meta_df$Prefix) 
meta_df$Growth <- growth_map[meta_df$Prefix]

colnames(FPRime) <- meta_df$label

assay_mat <- assay(FPRime)
colnames(assay_mat) <- meta_df$label

non_pdx_samples <- meta_df$Growth != "PDX" & !is.na(meta_df$Growth)
assay_mat <- assay_mat[, non_pdx_samples]
meta_df    <- meta_df[non_pdx_samples, ]

growth_factor <- factor(meta_df$Growth)
design_to_protect <- model.matrix(~ 0 + growth_factor)
cell_line_batch <- factor(meta_df$Prefix)

residual_matrix <- removeBatchEffect(
  assay_mat, 
  batch = cell_line_batch, 
  design = design_to_protect
)

residual_df <- as.data.frame(residual_matrix)
residual_df$protein_id <- rownames(assay_mat) 

row_metadata <- as.data.frame(rowData(FPRime))
row_metadata$protein_id <- rownames(FPRime)

gene_symbol_map <- row_metadata %>% 
  dplyr::select(protein_id, name)

padj_cols <- grep("GR_vs_.*_IgG_p\\.adj|IgG_vs_.*_GR_p\\.adj", colnames(row_metadata), value = TRUE)

contrast_pvals_long <- row_metadata %>%
  dplyr::select(protein_id, all_of(padj_cols)) %>%
  pivot_longer(
    cols = -protein_id,
    names_to = "contrast_name",
    values_to = "p_adj"
  ) %>%
  mutate(Prefix = sub("_.*", "", contrast_name)) %>%
  mutate(Prefix = trimws(sub(" $", "", Prefix))) %>%
  mutate(Growth = growth_map[Prefix]) %>%
  filter(!is.na(Growth))

long_adjusted <- residual_df %>%
  pivot_longer(
    cols = -protein_id,
    names_to = "sample_name",
    values_to = "Adjusted_Intensity"
  ) %>%
  filter(grepl("GR", sample_name)) %>%
  mutate(
    Prefix = trimws(sub("_.*", "", sample_name)),
    Growth = growth_map[Prefix]
  ) %>%
  filter(!is.na(Growth))

# 1-to-many warning
matched_df <- long_adjusted %>%
  left_join(
    contrast_pvals_long %>%
      dplyr::select(protein_id, Prefix, p_adj),
    by = c("protein_id", "Prefix")
  )

plot_ready <- matched_df %>%
  group_by(protein_id, Growth) %>%
  summarise(
    Mean_Adjusted = mean(Adjusted_Intensity, na.rm = TRUE),
    Is_Significant = sum(p_adj < 0.01, na.rm = TRUE) > 30,
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Growth,
    values_from = c(Mean_Adjusted, Is_Significant),
    names_glue = "{Growth}_{.value}"
  ) %>%
  drop_na(
    Suspension_Mean_Adjusted,
    Adherent_Mean_Adjusted
  ) %>%
  mutate(
    Protein_Class = if_else(
      Suspension_Is_Significant & Adherent_Is_Significant,
      "Conserved GR binder",
      "Other"
    )
  ) %>%
  left_join(gene_symbol_map, by = "protein_id")

# ==========================================
# 7. Label proteins of interest
# ==========================================

interesting <- c(
  "NR3C1",
  "NCOA1",
  "NCOA2",
  "NCOA3",
  "MED1",
  "EP300",
  "CREBBP",
  "SMARCA4"
  )

labels_df <- plot_ready %>%
  filter(name %in% interesting)

# ==========================================
# 8. Plot
# ==========================================

ggscatter(
  plot_ready,
  x = "Suspension_Mean_Adjusted",
  y = "Adherent_Mean_Adjusted",
  color = "Protein_Class",
  palette = c("#4C78A8","grey75"),
  alpha = 0.55,
  add = "reg.line",
  conf.int = TRUE,
  add.params = list(
    color = "firebrick",
    fill = "grey90"
  ),
  cor.coef = TRUE,
  cor.method = "pearson",
  xlab = "Mean cell-line-adjusted DIANNeR intensity (Suspension)",
  ylab = "Mean cell-line-adjusted DIANNeR intensity (Adherent)"
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = 2,
    colour = "grey50"
  ) +
  geom_text_repel(
    data = labels_df,
    aes(label = name),
    size = 3.5,
    fontface = "bold",
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "black",
    segment.size = 0.4,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
   theme(legend.position = "bottom")+labs(color = NULL)


