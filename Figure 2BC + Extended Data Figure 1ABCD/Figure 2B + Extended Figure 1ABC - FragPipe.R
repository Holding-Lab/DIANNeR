#####
# Figure 2B, Extended Data Figure 1A-C
# 
# 23/7/2026
####


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




plotRimePCAs<-function(Rdata,pca_data,PCx=3,PCy=4) {

  
  pca_plot <- ggplot(pca_data, aes(x = .data[[paste0("PC", PCx)]], y = .data[[paste0("PC", PCy)]], color = Origin,
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
      x = p$labels$x,
      y = p$labels$y
    )
  
  pca_plot$data <- subset(
    pca_plot$data, 
    !grepl("^(BB7|BB3RC31|HBC34)", sample)
  )
  
  return(pca_plot)
}
  
#PCA 1,2 - Figure 2B
plotRimePCAs(Rdata, pca_data, 1, 2)


#PCA 3+ - Extended Data Figure 1 A-C
library(patchwork)

plots <- list()
pc_pairs <- list(c(3,4), c(5,6), c(7,8), c(9,10))

for (pair in pc_pairs) {
  plots[[length(plots) + 1]] <- plotRimePCAs(Rdata, pca_data, pair[1], pair[2])
}

grid <- wrap_plots(plots, ncol = 2) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
grid


#ggsave(
#  filename = "Figure 3 - Supplimentary - PCA3-10.svg",
#  plot = grid,
#  width = 190,
#  height = 260,
#  units = "mm",
#  dpi = 300,
#  device = "svg"
#)
