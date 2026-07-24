#####
# Figure 3EFG + Extended Data Figure 3AB
# Transcript correlation plots
# 24/7/2026
####

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



# Load RNA-seq for cell lines
txData<-read.delim("GRSeq/counts/featureCounts.txt", sep="\t", skip=1)
txData$Geneid_clean <- sub("\\..*", "", txData$Geneid)

#Get gene symbols for RNA-seq
TxGenes <- mapIds(org.Hs.eg.db,
                keys = txData$Geneid_clean,
                column = "SYMBOL",
                keytype = "ENSEMBL",
                multiVals = "first")

txData$symbol <- TxGenes


sampleTable<-read.csv(file="GRSeq/sampletable.csv")

sampleTable
colnames(txData)

txDataTidy <- txData
sampleCols <- colnames(txDataTidy)[7:(ncol(txDataTidy)-2)]

# Clean up the BAM column names to match
cleanNames <-  gsub("\\.", "-", 
                   gsub("_Aligned\\.sortedByCoord\\.out\\.bam", "", 
                   gsub("\\.\\.BAMs.", "", sampleCols, fixed = FALSE)
                  ))
  
# Assign these cleaned names back to txDataTidy
colnames(txDataTidy)[7:(ncol(txDataTidy)-2)] <- cleanNames

# Verify alignment with your sample table
all(cleanNames %in% sampleTable$X)
# should return TRUE if everything matches

samples<-sampleTable[sampleTable$condition == "DMSO" & sampleTable$cellline %in% c("Jurkat","M231","KMBC2","MCF7") ,]
         
#Get only untreated cells 
txDataDMSO<-txDataTidy[,colnames(txDataTidy) %in% c("Geneid_clean", "symbol",samples$X)]


#Get transcript length
#mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#tx_lengths <- getBM(
#  attributes = c("ensembl_gene_id", "ensembl_transcript_id", "transcript_length"),
#  mart = mart
#)

#Get average Transcript length
#gene_lengths <- aggregate(
#  transcript_length ~ ensembl_gene_id,
#  data = tx_lengths,
#  FUN = mean
#)
#saveRDS(gene_lengths, "tx_lengths.rds")
gene_lengths<-readRDS("tx_lengths.rds")

txDataDMSOfiltered <- txDataDMSO[txDataDMSO$Geneid_clean %in% gene_lengths$ensembl_gene_id,]

rownames(gene_lengths) <- gene_lengths$ensembl_gene_id
lengths<-gene_lengths[txDataDMSOfiltered$Geneid_clean,]$transcript_length

TPMsDMSO<-cbind(
            convertCounts(counts=as.matrix(txDataDMSOfiltered[,1:(ncol(txDataDMSOfiltered)-2)]),
              geneLength=lengths,
              unit ="TPM",
              log = FALSE,
              normalize = "none"),
            txDataDMSOfiltered[,(ncol(txDataDMSOfiltered)-1):ncol(txDataDMSOfiltered)]
            )

TPMsDMSOtoPlot<-TPMsDMSO[TPMsDMSO$symbol %in% levels(cgDataFiltered$gene),]  
head(TPMsDMSOtoPlot)
#### Plot simlar to figure 4

num_cols <- setdiff(colnames(TPMsDMSOtoPlot), c("Geneid_clean", "symbol"))
tpm_mat  <- as.matrix(TPMsDMSOtoPlot[, num_cols, drop = FALSE])
rownames(tpm_mat) <- make.unique(TPMsDMSOtoPlot$symbol)

#log transform
base_names <- sub("-.*", "", colnames(tpm_mat))

tpm_df <- as.data.frame(tpm_mat)
tpm_avg <- sapply(unique(base_names), function(cell) {
  cols <- which(base_names == cell)
  rowMeans(tpm_df[, cols, drop = FALSE])
})

tpm_avg <- as.matrix(tpm_avg)
tpm_log <- log2(tpm_avg + 1)
tpm_z_by_sample <- apply(tpm_log, 2, scale)
rownames(tpm_z_by_sample)<-rownames(tpm_log)
breaks <- c(-1, -0.5, -0.25 ,0, 0.25, 0.5, 1)

my_colors_custom <- colorRampPalette(c("purple4", "mediumpurple2", "thistle","white", "lightgreen", "seagreen3", "darkgreen"))(length(breaks) - 1)


#Plot expression heatmap
# txHeatmap <- grid::grid.grabExpr(
#   pheatmap::pheatmap(
#     mat = tpm_z_by_sample,
#     cluster_rows = TRUE,
#     cluster_cols = TRUE,
#     show_rownames = TRUE,
#     show_colnames = TRUE,
#     main = "z-score for expression  GR interacting proteins",
#     color = my_colors_custom,
#     breaks = breaks,
#     labels_row = rownames(tpm_z_by_sample),   # symbols already set as rownames
#     fontsize_row = 5
#   )
# )
# 
# txHeatmap 

#Figure <- ggarrange(ncol=2,fullHM,pca123)


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




#scColors <- colorRampPalette(c("white", "thistle", "mediumpurple2", "purple4"))(length(breaks) - 1)
# 
# p <- grid::grid.grabExpr(
#   pheatmap(
#     cppt_z_by_sample,
#     cluster_rows = TRUE,
#     cluster_cols = TRUE,
#     show_rownames = TRUE,
#     show_colnames = TRUE,
#     main = "mean_logCPTT (all cells) by tissue",
#     color = scColors,
#     breaks = breaks,
#     fontsize_row = 5
#   )
# )
# p

#Combined both on same plot.

tissue_map_tx <- c(
  "blood" = "CD4+",
  "breast" = "Breast Epithelial",
  "ureter" = "Urothelial",
  "MCF7" = "MCF7",
  "M231" = "MDA-MB-231",
  "Jurkat" = "Jurkat",
  "KMBC2" = "KMBC2",
  "MCF10A" = "MCF10A",
  "BB7" = "PDX",
  "BB3RC31" = "PDX",
  "HBC34" = "PDX"
)

combined_plot<-cbind(cppt_z_by_sample[rownames(cppt_z_by_sample),],
      tpm_z_by_sample[rownames(cppt_z_by_sample),]
      )

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
txHeatmap 

#write.csv(combined_plot, file="Supplimentary Table S3/TableS3.csv")

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

# p<-grid::grid.grabExpr(pheatmap(
#   mat =DEResultsMat_z_by_sample,
#   cluster_rows = TRUE,
#   cluster_cols = TRUE,
#   show_rownames = TRUE,
#   show_colnames = TRUE,
#   main = "zScore RIME for correlation plot",
#   color = my_colors_custom,
#   breaks = breaks, # Uncomment to use custom breaks
#   labels_row=rownames(DEResultsMat),
#   #  scale = "row" # Uncomment to apply row scaling (Z-score)
#   fontsize_row = 5)
# )

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

  test <- cor.test(x, y, method = "pearson")
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

plot

#ggsave("FigureN.svg",plot,units="mm",dpi=300,width=180, height=215)



#HIPPO


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


Pa<-scatter_supp


#Linage for lympoid


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




ggarrange(Pa  + theme(legend.position="none") ,scatter_supp + theme(legend.position="none"),widths = c(0.66, 1), heights=c(1,1))
#ggsave("Merged.svg",units="mm",dpi=300,width=375, height=75)
