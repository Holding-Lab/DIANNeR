library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(tidyr)
library(biomaRt)
library(DGEobj.utils)
library(pheatmap)


#Load CellxGene data
cgData<-read.delim("CellxGene/census_v4_summary_incremental.tsv", sep="\t")

#Filter for the cells types only in the tissue of interest for reduce figure complexity
cgDataFiltered<- cgData[(
                    (cgData$tissue == "breast" & 
                       cgData$cell_type_group == "luminal_mammary_epithelial_cell"
                     ) |
                    (cgData$tissue == "blood"  &
                       cgData$cell_type_group == "cd4_tcell"
                     ) |
                    (cgData$tissue == "ureter"  & 
                       cgData$cell_type_group == "urothelial_cell"
                     )
                 ),]

#Get gene symbols
genes <- mapIds(org.Hs.eg.db,
                keys = cgDataFiltered$gene_ensembl,
                column = "SYMBOL",
                keytype = "ENSEMBL",
                multiVals = "first")

cgDataFiltered$gene <- genes
cgDataFiltered<-cgDataFiltered[!cgDataFiltered$gene %in% c("GAPDH", "ACTB","RPL27", "RPLP1"),]
cgDataFiltered$gene <- factor(cgDataFiltered$gene, levels = unique(cgDataFiltered$gene))


#### Plots single cell data

cppt_df <- cgDataFiltered %>%
  dplyr::select(c('tissue','gene','mean_logCPTT_allcells_dataset_avg')) %>%
  tidyr::pivot_wider(names_from =  gene, values_from = mean_logCPTT_allcells_dataset_avg )

cppt_mat <- t(as.matrix(cppt_df[,-1]))
colnames(cppt_mat)<- t(cppt_df[1:3,1])

cppt_z_by_sample <- apply(cppt_mat, 2, scale)
rownames(cppt_z_by_sample )<-rownames(cppt_mat)
txbreaks <- seq(-1, 1, length.out = 101)

scColors <- colorRampPalette(c("#4575B4", "#91BFDB", "white", "#FEE090", "#D73027"))(length(txbreaks) - 1)



#Combined both on same plot.

tissue_map_tx <- c(
  "CD4." = "CD4+",
  "Breast.Epithelial" = "Breast Epithelial",
  "Urothelial" = "Urothelial",
  "MCF7" = "MCF7",
  "MDA.MB.231" = "MDA-MB-231",
  "Jurkat" = "Jurkat",
  "KMBC2" = "KMBC2",
  "MCF10A" = "MCF10A",
  "BB7" = "PDX",
  "BB3RC31" = "PDX",
  "HBC34" = "PDX"
)

combined_plot<-read.csv(file="TableS3/TableS3.csv", row.names = 1)

colnames(combined_plot)<-tissue_map_tx[colnames(combined_plot)]
txHeatmap  <- grid::grid.grabExpr(
  pheatmap(
    combined_plot,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    color = scColors,
    breaks = txbreaks,
    fontsize_row = 5
  )
)


###Load proteomics

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

DEResults<-read.csv("FragPipe Output/DE_results.csv")

contrasts<-colnames(DEResults)
conditions<-unique(read.csv('FragPipe Input/E733_All_FP-A_experimental_annotation_pool_removed.tsv', sep="\t")$condition)
conditionsGRonly<-conditions[grepl('_GR$',conditions)]
tissues<-gsub('.{3}$', '', conditionsGRonly)
tissues[tissues=='231']<-"X231" #Safe naming
contrastsLFC<-paste0(tissues,"_IgG_vs_",tissues,"_GR_log2.fold.change")
contrastsPval<-paste0(tissues,"_IgG_vs_",tissues,"_GR_p.val")

contrastsLFC %in%   contrasts
contrastsPval %in%   contrasts

DEResultsFcGR<-  -DEResults[, contrastsLFC] 
#Remove any protein with LFC < 1 as it is IgG Specific
DEResultsFcGR[DEResultsFcGR < 0]<-0

pca_input <- t(DEResultsFcGR)        
colnames(pca_input)<-DEResults$Gene.Name
rownames(pca_input)<-sapply(strsplit(rownames(pca_input), "_"), "[[", 1)

#levels(cgDataFiltered$gene)[!levels(cgDataFiltered$gene) %in% colnames(pca_input) ]
#levels(cgDataFiltered$gene)[levels(cgDataFiltered$gene) %in% colnames(pca_input) ]




DEResultsMat<-t(pca_input)[rownames(t(pca_input)) %in% genes,!colnames(t(pca_input)) %in% "MCF10A"]

breaks <- c(-1, -0.75, -0.5, -0.25, 0.25,0.5, 0.75, 1)
my_colors_custom <- colorRampPalette(c( "gold1", "white",  "darkblue"))(length(breaks) - 1)


DEResultsMat_z_by_sample <- apply( as.matrix(DEResultsMat), 2, scale)
colnames(DEResultsMat_z_by_sample)<-tissue_map[colnames(DEResultsMat_z_by_sample)]
rownames(DEResultsMat_z_by_sample)<-rownames(DEResultsMat)


#Correlation plots

results <- data.frame(Protein = character(),
                      Correlation = numeric(),
                      P_value = numeric(),
                      stringsAsFactors = FALSE)

common_proteins <- intersect(rownames(DEResultsMat_z_by_sample), rownames(combined_plot))
df1 <- DEResultsMat_z_by_sample[common_proteins, ]
df2 <- combined_plot[common_proteins, ]

for (p in common_proteins) {
  x <- df1[p, colnames(df1)] 
  y <- df2[p,colnames(df1) ]

  test <- cor.test(x, as.numeric(y), method = "pearson")
  results <- rbind(results,
                   data.frame(Protein = p,
                              Correlation = unname(test$estimate),
                              P_value = test$p.value))

}

results <- results[order(-results$Correlation), ]


### Correlation results plots

#Distributions of correlations
distribution<-
  ggplot(results, aes(x = Correlation)) +
  geom_density(fill = "steelblue", alpha = 0.7, adjust = 0.6, outline.type = "full",bounds = c(-1,1)) +
  geom_vline(xintercept = 0, linetype = 2) +
  coord_cartesian(xlim = c(-1, 1)) +
  labs(   x = "Pearson r", y = "Density") +
  theme_minimal()


#Lollipop
to_plot <- results %>%
  arrange(desc(Correlation)) %>%
  mutate(Protein = factor(Protein, levels = Protein[order(Correlation)]),
         sig = P_value < 0.05)

lollipop<-ggplot(to_plot, aes(x = Protein, y = Correlation, color = sig)) +
  geom_segment(aes(xend = Protein, y = 0, yend = Correlation)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "deeppink")) +
  coord_flip() +
  labs(
       x = "", y = "", color = "Significant") +
  theme_minimal()+
  theme(text = element_text(size = 6.25))



# choose proteins
proteins <- c("GRHL1","GRHL2","ZEB1","ZEB2","PGR","ESR1")

# build long data frame (one row per cell type per protein)
df_list <- lapply(proteins, function(p) {
  x <- as.numeric(DEResultsMat_z_by_sample[p, , drop = TRUE])
  y <- as.numeric(combined_plot[p, colnames(DEResultsMat_z_by_sample), drop = TRUE])
  data.frame(Protein = p, CellType = colnames(DEResultsMat_z_by_sample), x = x, y = y)
})
dfp <- do.call(rbind, df_list)

# per-protein stats + positions to place the labels nicely
stats <- dfp %>%
  group_by(Protein) %>%
  summarise(
    r = cor(x, y, use = "pairwise.complete.obs", method = "pearson"),
    p = cor.test(x, y, method = "pearson")$p.value,
    x_lab = min(x, na.rm = TRUE),
    y_lab = max(y, na.rm = TRUE),
    .groups = "drop"
  )

# plot (facet by protein)


dfp$CellType <- factor(dfp$CellType, levels = c(
  "CD4+", "Jurkat",               # blood / immune
  "MCF7", "MDA-MB-231", "Breast Epithelial",  # breast
  "KMBC2", "Urothelial"           # biliary / bladder
))

# plot (facet by protein)
scatter<-ggplot(dfp, aes(x = x, y = y, colour = CellType, label = CellType)) +
  geom_point(size = 0.5, alpha = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  coord_equal() +
  facet_wrap(~ Protein, ncol = 3, scales = "fixed")  +
  scale_color_manual(
    values = c(
      "CD4+"              = "#3B4CC0",  # deep blue
      "Jurkat"            = "#8B1A9A",  # strong purple
      "MCF7"              = "#E41A1C",  # bright red
      "MDA-MB-231"        = "#FB8072",  # coral pink
      "Breast Epithelial" = "#FDB462",  # peach/orange
      "KMBC2"             = "#FFD92F",  # vivid yellow
      "Urothelial"        = "#B3DE69"   # green-yellow
    ))+
  labs( x = "DIANNeR z-Score", y = "Expression (RNA) z-Score", colour = "Cell type") +
  theme_minimal(base_size = 10)

library(ggpubr)
#FigureBCE<-ggarrange(nrow=3,distribution,lollipop,scatter)
FigureBC<-ggarrange(nrow=2,lollipop,scatter, heights= c(2, 1))
plot<-ggarrange(ncol=2,txHeatmap,FigureBC)


#Figure 8 Combined
plot

ggsave("Figure8_merged.svg",plot,units="mm",dpi=300,width=180, height=215)


## Figure S10: AR 

# choose proteins
proteins <- c("AR")

# build long data frame (one row per cell type per protein)
df_list <- lapply(proteins, function(p) {
  x <- as.numeric(DEResultsMat_z_by_sample[p, , drop = TRUE])
  y <- as.numeric(combined_plot[p, colnames(DEResultsMat_z_by_sample), drop = TRUE])
  data.frame(Protein = p, CellType = colnames(DEResultsMat_z_by_sample), x = x, y = y)
})
dfp <- do.call(rbind, df_list)

# per-protein stats + positions to place the labels nicely
stats <- dfp %>%
  group_by(Protein) %>%
  summarise(
    r = cor(x, y, use = "pairwise.complete.obs", method = "pearson"),
    p = cor.test(x, y, method = "pearson")$p.value,
    x_lab = min(x, na.rm = TRUE),
    y_lab = max(y, na.rm = TRUE),
    .groups = "drop"
  )

# plot (facet by protein)


dfp$CellType <- factor(dfp$CellType, levels = c(
  "CD4+", "Jurkat",               # blood / immune
  "MCF7", "MDA-MB-231", "Breast Epithelial",  # breast
  "KMBC2", "Urothelial"           # biliary / bladder
))

# plot (facet by protein)
scatter_supp<-ggplot(dfp, aes(x = x, y = y, colour = CellType, label = CellType)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  coord_equal() +
  facet_wrap(~ Protein, ncol = 3, scales = "fixed")  +
  scale_color_manual(
    values = c(
      "CD4+"              = "#3B4CC0",  # deep blue
      "Jurkat"            = "#8B1A9A",  # strong purple
      "MCF7"              = "#E41A1C",  # bright red
      "MDA-MB-231"        = "#FB8072",  # coral pink
      "Breast Epithelial" = "#FDB462",  # peach/orange
      "KMBC2"             = "#FFD92F",  # vivid yellow
      "Urothelial"        = "#B3DE69"   # green-yellow
    ))+
  labs( x = "DIANNeR z-Score", y = "Expression (RNA) z-Score", colour = "Cell type") +
  theme_minimal(base_size = 10)

scatter_supp

ggsave("FigureS10.svg",scatter_supp,units="mm",dpi=300,width=150, height=150)

mean(dfp[c(2,3,5,6,7),]$y)
#[1] 1.214344 # mean epithelial expression
mean(dfp[c(1,4),]$y)
#[1] 0.5684459


mean(dfp[c(2,3,5,6,7),]$x)
#[1] 1.214344 # mean epithelial expression
mean(dfp[c(1,4),]$x)
#[1] 0.5684459  # mean lymphoid expression

DEResults



#Figure S9


# choose proteins
proteins <- c("TEAD4","WWTR1")

# build long data frame (one row per cell type per protein)
df_list <- lapply(proteins, function(p) {
  x <- as.numeric(DEResultsMat_z_by_sample[p, , drop = TRUE])
  y <- as.numeric(combined_plot[p, colnames(DEResultsMat_z_by_sample), drop = TRUE])
  data.frame(Protein = p, CellType = colnames(DEResultsMat_z_by_sample), x = x, y = y)
})
dfp <- do.call(rbind, df_list)

# per-protein stats + positions to place the labels nicely
stats <- dfp %>%
  group_by(Protein) %>%
  summarise(
    r = cor(x, y, use = "pairwise.complete.obs", method = "pearson"),
    p = cor.test(x, y, method = "pearson")$p.value,
    x_lab = min(x, na.rm = TRUE),
    y_lab = max(y, na.rm = TRUE),
    .groups = "drop"
  )

# plot (facet by protein)


dfp$CellType <- factor(dfp$CellType, levels = c(
  "CD4+", "Jurkat",               # blood / immune
  "MCF7", "MDA-MB-231", "Breast Epithelial",  # breast
  "KMBC2", "Urothelial"           # biliary / bladder
))

# plot (facet by protein)
scatter_supp<-ggplot(dfp, aes(x = x, y = y, colour = CellType, label = CellType)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  coord_equal() +
  facet_wrap(~ Protein, ncol = 3, scales = "fixed")  +
  scale_color_manual(
    values = c(
      "CD4+"              = "#3B4CC0",  # deep blue
      "Jurkat"            = "#8B1A9A",  # strong purple
      "MCF7"              = "#E41A1C",  # bright red
      "MDA-MB-231"        = "#FB8072",  # coral pink
      "Breast Epithelial" = "#FDB462",  # peach/orange
      "KMBC2"             = "#FFD92F",  # vivid yellow
      "Urothelial"        = "#B3DE69"   # green-yellow
    ))+
  labs( x = "DIANNeR z-Score", y = "Expression (RNA) z-Score", colour = "Cell type") +
  theme_minimal(base_size = 10)

scatter_supp

ggsave("FigureS9.svg",scatter_supp,units="mm",dpi=300,width=150, height=150)



#Figure S11:Linage for lympoid


# choose proteins
proteins <- c("FOXP3", "IKZF3",  "SATB1")

# build long data frame (one row per cell type per protein)
df_list <- lapply(proteins, function(p) {
  x <- as.numeric(DEResultsMat_z_by_sample[p, , drop = TRUE])
  y <- as.numeric(combined_plot[p, colnames(DEResultsMat_z_by_sample), drop = TRUE])
  data.frame(Protein = p, CellType = colnames(DEResultsMat_z_by_sample), x = x, y = y)
})
dfp <- do.call(rbind, df_list)

# per-protein stats + positions to place the labels nicely
stats <- dfp %>%
  group_by(Protein) %>%
  summarise(
    r = cor(x, y, use = "pairwise.complete.obs", method = "pearson"),
    p = cor.test(x, y, method = "pearson")$p.value,
    x_lab = min(x, na.rm = TRUE),
    y_lab = max(y, na.rm = TRUE),
    .groups = "drop"
  )

# plot (facet by protein)


dfp$CellType <- factor(dfp$CellType, levels = c(
  "CD4+", "Jurkat",               # blood / immune
  "MCF7", "MDA-MB-231", "Breast Epithelial",  # breast
  "KMBC2", "Urothelial"           # biliary / bladder
))

# plot (facet by protein)
scatter_supp<-ggplot(dfp, aes(x = x, y = y, colour = CellType, label = CellType)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  coord_equal() +
  facet_wrap(~ Protein, ncol = 3, scales = "fixed")  +
  scale_color_manual(
    values = c(
      "CD4+"              = "#3B4CC0",  # deep blue
      "Jurkat"            = "#8B1A9A",  # strong purple
      "MCF7"              = "#E41A1C",  # bright red
      "MDA-MB-231"        = "#FB8072",  # coral pink
      "Breast Epithelial" = "#FDB462",  # peach/orange
      "KMBC2"             = "#FFD92F",  # vivid yellow
      "Urothelial"        = "#B3DE69"   # green-yellow
    ))+
  labs( x = "DIANNeR z-Score", y = "Expression (RNA) z-Score", colour = "Cell type") +
  theme_minimal(base_size = 10)

scatter_supp

ggsave("FigureS11.svg",scatter_supp,units="mm",dpi=300,width=150, height=150)

