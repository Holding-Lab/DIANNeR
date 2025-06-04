library(FragPipeAnalystR)
library(ggplot2)
library(ggrepel)

load("FragPipe Output/Imputed_SE.RData")

p<-plot_pca(RData)

pca_data <- p$data
pca_data$Antibody <- ifelse(grepl("GR", pca_data$condition), "GR", 
                            ifelse(grepl("IgG", pca_data$condition), "IgG", 
                                   "Other"))

pca_data$Pool <- ifelse(grepl("Pool", pca_data$label), "Pool", "Indiviual")

pca_data$Prefix <- sub("_.*", "", pca_data$label)
pca_data$Prefix <- trimws(sub("_.*", "", pca_data$label)) #One sample has a tailing space

tissue_map <- c(
  "PrimaryCD4" = "CD4+",
  "PrimaryBreastEpis" = "Breast Epithelial",
  "PrimaryNHU" = "Urothelial",
  "MCF7" = "MCF7",
  "231" = "MDA-MB-231",
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
  "231" = "Breast",
  "MCF10A" = "Breast",
  "BB7" = "Breast",
  "BB3RC31" = "Breast",
  "HBC34" = "Breast",
  
  "KMBC2" = "Ureter",
  "PrimaryNHU" = "Ureter"
)

pca_data$Tissue <- tissue_map[pca_data$Prefix]
pca_data$Origin <- origin_map[pca_data$Prefix]




pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Origin,
                                 shape = Antibody,
                                 fill = Pool,
                                 label = Tissue)) +
  geom_point(size = 2, stroke = 1) +
  geom_text_repel(
    size = 2.8,           
    max.overlaps = 10
  ) +
  scale_shape_manual(values = c("GR" = 22, "IgG" = 21, "Other" = 23)) +
  scale_fill_manual(values = c("white", "grey50")) +
  theme_minimal(base_size = 10) +  # Set base font size to 12 pt
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10)
  ) +
  labs(
    title = "PCA of Protein Intensities by Antibody",
    x = p$labels$x,
    y = p$labels$y
  )

pca_plot

ggsave(
  filename = "Figure 3B - PCA plot.pdf",
  plot = pca_plot,
  width = 120,
  height = 100,
  units = "mm",
  dpi = 300,
  device = "pdf"
)



#Impute done in preprocessing? 
library(FragPipeAnalystR)
load("FragPipe Output/Imputed_SE.RData")
FPRime<-RData
#FPRime <- impute(FPRime, fun = "mixed")
FPRime <- test_limma(FPRime, type="all")
FPRime <- add_rejections(FPRime, alpha = 0.01, lfc = 1)

FPRime$Antibody <- ifelse(grepl("GR", colnames(FPRime)), "GR", 
                            ifelse(grepl("IgG", colnames(FPRime)), "IgG", 
                                   "Other"))

get_cluster_heatmap(FPRime, type="centered", 
                    indicate = c("Antibody"),
                    show_row_names = FALSE,
                    show_heatmap_legend = FALSE)

