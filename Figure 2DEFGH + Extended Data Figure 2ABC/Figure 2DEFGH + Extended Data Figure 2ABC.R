#####
# Figure 2DEFGH + Extended Data Figure 2ABC 
# Analysis of key contributing proteins to PC1-3
# 24/7/2026
####

#Average protein intensities (unsued)
library(tidyverse)
intensity<-read.delim("FragPipe Input/E733_DIA_Report_MinImp_pool_removed.tsv")


df_long <- intensity %>%
  pivot_longer(
    cols = 7:ncol(.),  # intensity columns
    names_to = "Sample",
    values_to = "Intensity"
  ) %>%
  mutate(
    # Extract IP type (IgG or GR)
    IP = case_when(
      str_detect(Sample, "_GR_") | str_detect(Sample, "_GR_IP") ~ "GR",
      str_detect(Sample, "_IgG_") | str_detect(Sample, "_IgG_IP") ~ "IgG",
      TRUE ~ "Other"
    ),
    # Extract cell line before the IP tag
    CellLine = Sample %>%
      str_remove("_GR_.*|_IgG_.*|_GR_IP_.*|_IgG_IP_.*") %>%
      str_replace_all("\\.+", "")  # remove any stray dots from names
  ) %>%
  filter(IP != "Other")  # keep only IgG and GR samples

df_long$log10Intensity <- log10(df_long$Intensity)

df_long <- df_long %>%mutate(
  CellLine = str_replace(CellLine, "^X231", "MDA-MB-231")
)



p<-ggplot(df_long, aes(x = interaction(CellLine, IP), y = log10Intensity, fill = IP)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.3) +
  scale_fill_manual(values = c("IgG" = "#1f77b4", "GR" = "#ff7f0e")) +
  labs(
    x = "Cell Line and IP Type",
    y = expression(log[10]*"(Protein Intensity)"),
    title = "Protein Intensities by Cell Line and IP Type"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=4))

#Unused
#ggsave(p, file="Rebuttal Figure X - Protein Intesities By Cell line.png",dpi=300, width=170, height=120, unit="mm")

df_long_pdx <- df_long %>%
  mutate(
    # Assign CellLine: if underscore exists, extract before it, else keep as is
    CellLine = if_else(str_detect(CellLine, "_"), str_remove(CellLine, "_.*$"), CellLine),
    # Assign PDX label for those 3 cell lines
    CellLine = if_else(CellLine %in% c("BB3RC31", "BB7", "HBC34"), "PDX", CellLine),
    # Extract replicate number from sample name where possible
    Replicate = str_extract(Sample, "_\\d+$") %>% str_remove("_")
  ) %>%
  # Now assign replicate numbers explicitly for those PDX samples based on original cell line
  mutate(
    Replicate = case_when(
      CellLine == "PDX" & str_detect(Sample, "BB3RC31") ~ "1",
      CellLine == "PDX" & str_detect(Sample, "BB7") ~ "2",
      CellLine == "PDX" & str_detect(Sample, "HBC34") ~ "3",
      TRUE ~ Replicate
    )
  )

df_summary <- df_long_pdx %>%
  group_by(CellLine, IP, Replicate) %>%
  summarise(
    AvgIntensity = mean(Intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(log10AvgIntensity = log10(AvgIntensity + 1e-6))

p2 <- ggplot(df_summary, aes(x = CellLine, y = log10AvgIntensity, color = IP)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.6, outlier.shape = NA) +  # boxplot without outliers
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 2, alpha = 0.8) +  # overlay points
  scale_color_manual(values = c("IgG" = "#1f77b4", "GR" = "#ff7f0e")) +
  labs(
    x = "Cell Line",
    y = expression(log[10]*"(Average Protein Intensity)"),
    title = "Average Protein Intensity per Replicate by Cell Line and IP"
  ) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2
##Figure 2DEFGH
library(RColorBrewer)
library(pheatmap)

DEResults<-read.csv("FragPipe Output/DE_results.csv")

contrasts<-colnames(DEResults)

conditions<-unique(read.csv('FragPipe Input/E733_All_FP-A_experimental_annotation_pool_removed.tsv', sep="\t")$condition)
conditionsGRonly<-conditions[grepl('_GR$',conditions)]
tissues<-gsub('.{3}$', '', conditionsGRonly)
tissues[tissues=='231']<-"X231" #Safe naming
contrastsLFC<-paste0(tissues,"_IgG_vs_",tissues,"_GR_log2.fold.change")
contrastsPval<-paste0(tissues,"_IgG_vs_",tissues,"_GR_p.val")

#Check contrasts exist
contrastsLFC %in%   contrasts
contrastsPval %in%   contrasts

#LFC is negative for NR3C1 as the contasts are the wrong way around
#DEResults[DEResults$Gene.Name == "NR3C1",]
#View(DEResults[1108,contrastsLFC])
DEResultsFcGR<-  -DEResults[, contrastsLFC] 
#Remove any protein with LFC < 1 as it is IgG Specific
DEResultsFcGR[DEResultsFcGR < 0]<-0



##PCA
library(ggplot2)
library(ggpubr)
pca_input <- t(DEResultsFcGR)        
colnames(pca_input)<-DEResults$Gene.Name


rownames(pca_input)<-sapply(strsplit(rownames(pca_input), "_"), "[[", 1)
pca_res <- prcomp(pca_input)

tissue_map <- c(
  "PrimaryCD4" = "CD4+",
  "PrimaryBreastEpis" = "Breast Epithelial",
  "PrimaryNHU" = "Urothelial",
  "MCF7" = "MCF7",
  "X231" = "MDA-MB-231",
  "Jurkat" = "Jurkat",
  "KMBC2" = "KMBC2",
  "MCF10A" = "MCF10A",
  "BB7" = "PDX",
  "BB3RC31" = "PDX",
  "HBC34" = "PDX"
)


origin_map <- c(
  "PrimaryCD4" = "Blood",
  "Jurkat" = "Blood",
  
  "PrimaryBreastEpis" = "Breast",
  "MCF7" = "Breast",
  "X231" = "Breast",
  "MCF10A" = "Breast",
  "BB7" = "Breast",
  "BB3RC31" = "Breast",
  "HBC34" = "Breast",
  
  "KMBC2" = "Ureter",
  "PrimaryNHU" = "Ureter"
)

pca_df <- data.frame(
  Sample = tissue_map[rownames(pca_res$x)],
  Origin= origin_map[rownames(pca_res$x)],
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  PC3 = pca_res$x[, 3],
  PC4 = pca_res$x[, 4]
)

pca12<-ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample, color=Origin)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.5) +
  labs(title = "Principal Components 1 & 2",
       x = paste0("PC1 (", round(100 * summary(pca_res)$importance[2, 1], 1), "% variance)"),
       y = paste0("PC2 (", round(100 * summary(pca_res)$importance[2, 2], 1), "% variance)")) +
  theme_pubr()

pca23<-ggplot(pca_df, aes(x = PC2, y = PC3, label = Sample, color=Origin)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.5) +
  labs(title = "Principal Components 2 & 3",
       x = paste0("PC2 (", round(100 * summary(pca_res)$importance[2, 2], 1), "% variance)"),
       y = paste0("PC3 (", round(100 * summary(pca_res)$importance[2, 3], 1), "% variance)")) +
  theme_pubr()


pca123<-ggarrange(nrow=2,pca12,pca23)


loadings <- pca_res$rotation
top_PC1 <- sort(abs(loadings[, 1]), decreasing = TRUE)[1:25]
top_PC1_proteins <- names(top_PC1)

top_PC2 <- sort(abs(loadings[, 2]), decreasing = TRUE)[1:25]
top_PC2_proteins <- names(top_PC2)

top_PC3 <- sort(abs(loadings[, 3]), decreasing = TRUE)[1:25]
top_PC3_proteins <- names(top_PC3)
top_PC3_proteins[top_PC3_proteins=="NR3C1;NR3C2"]<-"NR3C2" #Mannually asigned as NR3C1 is the bait

PCproteins<-c(top_PC1_proteins,top_PC2_proteins,top_PC3_proteins)


#PCA Heatmap
DEResultsFcGRPvalFiltered<-DEResultsFcGR[apply(DEResultsFcGR > 2, 1, any) &
                                           apply(DEResults[, contrastsPval] < 0.01, 1 , any)
                                         ,]

DEResultsPC<-DEResultsFcGRPvalFiltered[DEResults$Gene.Name[as.numeric(rownames(DEResultsFcGRPvalFiltered))] %in% PCproteins,]

breaks <- c(0, 2, 4, 8, max(DEResultsPC) + 0.1) # Example breaks
my_colors_custom <- colorRampPalette(c("white", "lightblue", "steelblue", "darkblue"))(length(breaks) - 1)


colnames(DEResultsPC)<-
  tissue_map[sapply(strsplit(colnames(DEResultsPC), "_"), "[[", 1)]

p<-grid::grid.grabExpr(pheatmap(
  mat = as.matrix(DEResultsPC),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Heatmap of Top PC1-3 Contributors",
  color = my_colors_custom,
  breaks = breaks, # Uncomment to use custom breaks
  labels_row=DEResults$Gene.Name[as.numeric(rownames(DEResultsPC))],
  #  scale = "row" # Uncomment to apply row scaling (Z-score)
  fontsize_row = 5)
 )



library(ggsci)
#Bar plot
plot_gene_expression <- function(gene_expression_matrix, gene_names) {
  gene_names <- as.character(gene_names)
  gene_labels_all <- DEResults$Gene.Name[as.numeric(rownames(gene_expression_matrix))]
  gene_expression_subset <- gene_expression_matrix[gene_labels_all %in% gene_names, , drop = FALSE]
  gene_labels_found <- gene_labels_all[gene_labels_all %in% gene_names]
  
  gene_expression_subset <- as.matrix(gene_expression_subset)
  storage.mode(gene_expression_subset) <- "numeric"
  rownames(gene_expression_subset) <- gene_labels_found
  df <- data.frame(
    Sample = rep(colnames(gene_expression_subset), each = nrow(gene_expression_subset)),
    Expression = as.numeric(gene_expression_subset),
    Gene = rep(rownames(gene_expression_subset), times = ncol(gene_expression_subset))
  )
  
  ordered_levels <- c("CD4+", "Jurkat", "MCF7", "Breast Epithelial", 
                      "Urothelial", "KMBC2", "MDA-MB-231", "MCF10A")
  df$Sample <- factor(df$Sample , levels = ordered_levels)
  
  # Plot
  p <- ggplot(df, aes(x = Sample, y = Expression, fill = Gene)) +
    geom_bar(stat = "identity", position = "dodge",color='black') +
    labs(
     # title = paste("GR interaction with", paste(gene_names, collapse = ", ")),
      x = "",
      y = "Log2FC over IgG Control"
    ) + ylim(0, 10) +
    theme_pubr() +
    scale_fill_jco() +
    geom_hline(yintercept = 2, linetype = "dashed")
    
  return(p)
}

# Plot ESR1 - Added in AR and PGR for completeness, but not part of PC1-3 contributors in heatmap
DEResultsBarPlot<-DEResultsFcGRPvalFiltered
colnames(DEResultsBarPlot)<-sapply(strsplit(colnames(DEResultsFcGRPvalFiltered), "_"), "[[", 1)
colnames(DEResultsBarPlot)<-tissue_map[colnames(DEResultsBarPlot)]
esr1_plot <- plot_gene_expression(DEResultsBarPlot, gene_names = c("ESR1","PGR","AR"))

# Plot PPARG
#pparg_plot <- plot_gene_expression(DEResultsPC, gene_names = "PPARG")

# Plot E2F2 and E2F3 together (your example had them joined with ;)
e2f_plot <- plot_gene_expression(DEResultsPC, gene_names = c("E2F2;E2F3"))

# Plot GRHL1 and GRHL2 together
grhl_plot <- plot_gene_expression(DEResultsPC, gene_names = c("GRHL1", "GRHL2","ZEB1","ZEB2"))

#Total interactions numbers
fc_pass <- DEResultsFcGR > 2
pval_pass <- DEResults[, contrastsPval] < 0.01
col_sums <- colSums(fc_pass & pval_pass)
names(col_sums)<-
  tissue_map[sapply(strsplit(names(col_sums), "_"), "[[", 1)]
counts_df <- data.frame(
  Model = names(col_sums),
  Count = as.numeric(col_sums)
)

ordered_levels <- c("CD4+", "Jurkat", "MCF7", "Breast Epithelial", 
                    "Urothelial", "KMBC2", "MDA-MB-231", "MCF10A")
counts_df$Model <- factor(counts_df$Model, levels = ordered_levels)

interctors_plot <-
  ggplot(counts_df, aes(x = Model, y = Count)) +
  geom_bar(stat = "identity", position = "dodge", color='black') +
  labs(
    #title = "Specific GR Interactions (LFC > 2, adj-P value < 0.01)",
    x = "",
    y = "Number of Proteins"
  ) +
  theme_pubr() +
  scale_fill_jco() 

# Arrange all interaction plots
combined_plot <- ggarrange(interctors_plot,e2f_plot, esr1_plot, grhl_plot, nrow = 4)

### combine with heatmap
Figure2DEFGH<-ggarrange(p,combined_plot, ncol=2)

#ggsave("Figure4.svg",Figure4,units="mm",dpi=300,width=180, height=215)
Figure2DEFGH


#Full Heatmap for supplimentary
DEResultsFcGRNamed<-DEResultsFcGRPvalFiltered
colnames(DEResultsFcGRNamed)<-
  tissue_map[sapply(strsplit(colnames(DEResultsFcGRPvalFiltered), "_"), "[[", 1)]





fullHM<-grid::grid.grabExpr(pheatmap(
  mat = as.matrix(DEResultsFcGRNamed),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  main = "LFC > 2 and adjusted P-value < 0.01 for ≥ 1 sample",
  color = my_colors_custom,
  breaks = breaks, 
  fontsize_row = 5)
)

#Extended Data Figure 2ABC
ExtendedFigure2 <- ggarrange(ncol=2,fullHM,pca123)
ExtendedFigure2
#ggsave("ExtendedFigure2.svg",ExtendedFigure2,units="mm",dpi=300,width=180, height=215)

