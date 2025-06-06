ofDDA<-cbind(read.csv("FragPipe Output/E033_DDA_Fusion.csv"),
             acquisition="DDA",
             MS="Orbitrap Fusion")
ttDDA<-cbind(read.csv("FragPipe Output/E033_DDA_TimsTOF.csv"),
             acquisition="DDA",
             MS="timsTOF")
ttDIA<-cbind(read.csv("FragPipe Output/E033_DIA_TimsTOF.csv"),
             acquisition="DIA",
             MS="timsTOF")
colnames(ttDIA)[1]<-colnames(ttDDA)[1]
#Set DIA NA to zero?

ProteinCounts <-rbind(
                     c(nrow(ofDDA),"Orbitrap Fusion DDA","Total Protein Count"),
                     c(nrow(ttDDA),"timsTOF DDA","Total Protein Count"),
                     c(nrow(ttDIA),"timsTOF DIA","Total Protein Count"),
                     c(nrow(ofDDA[ofDDA$significant==TRUE,]),"Orbitrap Fusion DDA","Significant over IgG"),
                     c(nrow(ttDDA[ofDDA$significant==TRUE,]),"timsTOF DDA","Significant over IgG"),
                     c(nrow(ttDIA[ttDIA$significant==TRUE,]),"timsTOF DIA","Significant over IgG")
                     )

ProteinCounts_df <- as.data.frame(ProteinCounts, stringsAsFactors = FALSE)
colnames(ProteinCounts_df) <- c("Count", "Acquisition", "Category")
ProteinCounts_df$Count <- as.numeric(ProteinCounts_df$Count)

ProteinCounts_df$Category <- factor(ProteinCounts_df$Category,
                                       levels = c("Total Protein Count", "Significant over IgG"))


library(ggplot2)
library(ggpubr)
p <- ggbarplot(
  ProteinCounts_df,
  x = "Acquisition",
  y = "Count",
  fill = "Category",
  palette = c("#1b9e77", "#d95f02"),
  position = position_dodge(0.8),
  ylab = "Protein Count",
  xlab = "Acquisition Method"
) 

p


#Test also individual reps
library(dplyr)
combined_df <- bind_rows(ofDDA, ttDDA, ttDIA)

library(tidyr)
long_df <- combined_df %>%
  pivot_longer(cols = c(GR_1.Intensity, GR_2.Intensity, GR_3.Intensity, 
                        IgG_1.Intensity, IgG_2.Intensity, IgG_3.Intensity),
               names_to = "Sample",
               values_to = "Intensity")

long_df <- long_df %>%
  mutate(Antibody = ifelse(grepl("^GR", Sample), "GR", "IgG"),
         Replicate = gsub(".*_(\\d).*", "\\1", Sample))

summary_df <- long_df %>%
  group_by(acquisition, MS, Antibody, Replicate) %>%
  summarise(Protein_Count = sum(Intensity > 0, na.rm=TRUE), .groups = "drop")
summary_df$acquisition<-paste(summary_df$MS, summary_df$acquisition)
summary_df$Antibody<-as.factor(summary_df$Antibody)


# p1 <- ggplot(summary_df,
#              aes(
#                x = acquisition,
#                y = Protein_Count,
#                fill = Antibody
#              )) +
#   geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
#   geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
#               size = 2) +
#   scale_fill_manual(values = c("#1b9e77", "#d95f02")) +
#   labs(x = "Acquisition Method", y = "Protein Count", title = "Protein Count per Condition") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1 <- ggbarplot(summary_df,
               x = "acquisition",
               y = "Protein_Count",
               fill = "Antibody",
               palette = "jco",
               add = "mean_sd", add.params = list(group = "Antibody"),
               position = position_dodge(0.8),
               ylab = "Protein Count",
               xlab = "Acquisition Method"
             )



p1

stat.test <- compare_means(
  Protein_Count ~ Antibody, data = summary_df,
  method = "t.test",
  group.by = "acquisition"
)

stat.test$y.position<-max(summary_df$Protein_Count*1.2)

#Uncomment for IgG vs GR p-values
p2 <- p1# + stat_pvalue_manual(stat.test, x="acquisition", label = "p = {p.adj}")

stat.test_GR <- compare_means(
  Protein_Count ~ acquisition,  summary_df[summary_df$Antibody == 'GR',],
  method = "t.test",
)
stat.test_GR$y.position=max(summary_df$Protein_Count)*c(1.0,1.11,1.05)

p3<-p2 + stat_pvalue_manual(stat.test_GR,  label = "p = {p.adj}")




combined_df
combined_df$MSacquisition<-paste(combined_df$MS,combined_df$acquisition)

groups <- unique(combined_df$MSacquisition)

venn_list <- lapply(groups, function(g) {
  unique(combined_df$Protein.ID[combined_df$MSacquisition == g])
})

names(venn_list) <- c(
  "Orbitrap Fusion DDA\n(60 min acquisition)",
  "timsTOF DDA\n(14.4 min acquisition)",
  "timsTOF DIA\n(14.4 min acquisition)"
)

library(ggvenn)
pv<-ggvenn(venn_list, set_name_size = 5, fill_color = c("#0073C2FF", "#EFC000FF","#D0D0D0"))
pv


##Volcano

filter_idx <- combined_df$MSacquisition == "timsTOF DIA"
subset_df <- combined_df[filter_idx, ]
subset_df$negLog10_padj <- -log10(subset_df$GR_vs_IgG_p.adj)

genes_to_label <- c("STAT5A;STAT5B", 
                    "MCM5", 
                    "STAT1",
                    "SMARCE1",
                    "DNAJA1",
                    "SMARCC1;SMARCC2", 
                    "HMGB1;HMGB1P1", 
                    "SMARCD2",
                    "TRIM28",
                    "NCOA3", 
                    "MED24",
                    "HDAC2",
                    "NR3C1")

# Create label column: gene name if in list, else blank
subset_df$label <- ifelse(subset_df$Genes %in% genes_to_label, subset_df$Genes, "")

subset_df$color <- ifelse(subset_df$Genes == "NR3C1", "#5f02d9",
                   ifelse(subset_df$Genes %in% genes_to_label, "#d95f02",
                   ifelse(subset_df$GR_vs_IgG_p.adj < 0.05, "#f9bfa2AA", "grey60")))

subset_df$label_short <- sub(";.*", "", subset_df$label)


library(ggrepel)
volp<- ggplot(subset_df, aes(x = GR_vs_IgG_log2.fold.change, y = negLog10_padj)) +
  geom_point(aes(color = color), size = 2, show.legend = FALSE) +  # use color column
  geom_text_repel(aes(label = label_short), size = 3, max.overlaps = 1000) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  scale_color_identity() +  
  labs(x = "Log2 Fold Change (GR vs IgG)", 
       y = "-Log10 Adjusted p-value") +
  theme_pubr()

volp
#Multiplot
mp<-ggarrange(p3,p,pv,volp)
mp

# Save as SVG for manuscript
ggsave("Figure2ABCD_ProteinCounts_plot.svg", mp, width = 7, height = 5, units = "in", dpi = 300)



