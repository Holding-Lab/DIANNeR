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
combined_plot <- ggarrange(interctors_plot, esr1_plot,  e2f_plot, grhl_plot, nrow = 4)

### combine with heatmap
Figure4<-ggarrange(p,combined_plot, ncol=2)

ggsave("Figure4.svg",Figure4,units="mm",dpi=300,width=180, height=215)



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

FigureS7 <- ggarrange(ncol=2,fullHM,pca123)
ggsave("FigureS7.svg",FigureS7,units="mm",dpi=300,width=180, height=215)
