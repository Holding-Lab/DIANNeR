library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)

library(tidyverse)
library(dplyr)
library(ggpubr)
library(EnhancedVolcano)

####

DEResultsCombinedPDX <- read.csv('FragPipe PDX Combined/DE_results - PDX combined.csv')

DEResultsPDX <- read.csv('FragPipe PDX Combined/DE_results - PDX seperated.csv')

cols_to_add <- DEResultsCombinedPDX |> 
  dplyr::select(Gene.Name, Protein.ID,
         PrimaryCD4_plus_T_GR_vs_PDX_GR_IP_log2.fold.change, 
         PrimaryCD4_plus_T_GR_vs_PDX_GR_IP_p.adj,
         PDX_IgG_IP_vs_PDX_GR_IP_log2.fold.change,
         PDX_IgG_IP_vs_PDX_GR_IP_p.adj)

DEResultsPDX <- DEResultsPDX |> 
  left_join(cols_to_add, by = "Gene.Name")

CoreGRComplex<- c("NR3C1", "HOXA5" ,"ARNT","AHR","EP300","BCL11B","FOXP3")#c( "BCL11B", "SMARCD3", "ESR1", "PGR", "AR","FOXP3", "FOXA1")

filtered_genes <- DEResultsPDX |> 
  filter(Gene.Name %in% CoreGRComplex)

EpiComplex<-c("SMARCD3", "ESR1", "PGR", "AR",  "FOXA1","HOXA5")
LymphComplex<-c("FOXP3","BCL11B")

filtered_genes_colour <- filtered_genes |>
    mutate(
    "Associated Sample" = case_when(
      Gene.Name %in% EpiComplex  ~ "Breast Cancer PDX",
      Gene.Name %in% LymphComplex ~ "Normal CD4+ T Cell",
      TRUE                        ~ "GR Shared"
    )
  )



plotProteins<-function(filtered_genes_colour, Complex, contrast="PrimaryCD4_plus_T_GR_vs_PDX_GR_IP_log2.fold.change", rotated=TRUE, paletteColors= c("#decbe4bc", "#a0a0a0","#09546bff"), ylab_text=NULL) {
  filtered_genes_ordered <- filtered_genes_colour # |> arrange(desc(-.data[[contrast]]))
  filtered_genes_ordered <-   filtered_genes_ordered[filtered_genes_ordered$Gene.Name %in% Complex,]
  filtered_genes_ordered[contrast] <- -filtered_genes_ordered[contrast]
  
  p_adj_col <- paste0(gsub('.{16}$', '', contrast), "p.adj")
  
  filtered_genes_sig <- filtered_genes_ordered |> 
    arrange(desc(.data[[contrast]])) |>
    mutate(
      Significant = case_when(
        .data[[p_adj_col]] < 0.001 ~ "***",
        .data[[p_adj_col]] < 0.01  ~ "**",
        .data[[p_adj_col]] < 0.05  ~ "*",
        TRUE                       ~ ""
      )
    )
  
  filtered_genes_sig$Gene.Name <- factor(
    filtered_genes_sig$Gene.Name, 
    levels = rev(Complex)
  )
  if (is.null(ylab_text)) { ylab_text <- contrast}
  
  ggbarplot(filtered_genes_sig,
            "Gene.Name",
            contrast,
            orientation = "horizontal",
            fill = "Associated Sample",
            palette = paletteColors,
            xlab = "Protein Name",
            ylab = ylab_text,
            #ylab = "LFC Breast Cancer PDX (left) versus Primary CD4 + T Cell (right)",
  ) +coord_flip(ylim = c(-0.5, 8.5)) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "gray")+
    geom_hline(yintercept = c(0), linetype = "solid", color = "black")+
  geom_text(
    aes(
      label = Significant,
      y = .data[[contrast]] + (sign(.data[[contrast]]) * 1) # Shifts star slightly past the end of the bar
    ),
    vjust = 0.75,   # Centers the asterisks vertically with the bar
    hjust = 0.5     # Centers horizontally on its position
  )
}


p1<-plotProteins(filtered_genes_colour,contrast="PDX_IgG_IP_vs_PDX_GR_IP_log2.fold.change",  c("NR3C1", "EP300", "HOXA5","BCL11B","FOXP3"), ylab_text="LFC Enrichment  of GR interactors vs matched IgG control (PDX) ")
p2<-plotProteins(filtered_genes_colour,contrast="PrimaryCD4_plus_T_IgG_vs_PrimaryCD4_plus_T_GR_log2.fold.change",  c("NR3C1",  "EP300","HOXA5","BCL11B","FOXP3"), ylab_text="LFC Enrichment of GR interactors vs matched IgG contro ( CD4+ T Cell)")
p3<-plotProteins(filtered_genes_colour,contrast="PrimaryBreastEpis._IgG_vs_PrimaryBreastEpis._GR_log2.fold.change",  c("NR3C1",  "EP300","HOXA5","BCL11B","FOXP3"), ylab_text="LFC Enrichment of GR interactors vs matched IgG control (Breast Epithelium)")




#plotPDX<-function(filtered_genes_colour, contrast=c("HBC34_IgG_IP_vs_HBC34_GR_IP_log2.fold.change",
#                                                         "BB7_GR_IP_vs_BB3RC31_GR_IP_log2.fold.change",
#                                                         "BB3RC31_IgG_IP_vs_BB3RC31_GR_IP_log2.fold.change"),
#                                         genes) {
#   
#   #filtered_genes_ordered <- filtered_genes_colour |> arrange(desc(.data[[contrast[1]]]))
# 
#   
#   filtered_genes_colour_epi <- filtered_genes_colour |> filter(Gene.Name %in% genes)
#   #filtered_genes_colour_epi<-filtered_genes_colour[filtered_genes_colour$Gene.Name %in% EpiComplex,]
#   
#   
#   long_data <- pivot_longer(
#     filtered_genes_colour_epi[,c("Gene.Name",contrast)],
#     cols = all_of(contrast),       # The columns you want to pivot
#     names_to = "Contrast",         # The name for the new column containing contrast names
#     values_to = "Log2FC"           # The name for the new column containing the numbers
#   )
#   
# 
# 
#   
#   
#   
#   long_data$Log2FC    <- -long_data$Log2FC
#   long_data$Contrast <- gsub("_.*$","",long_data$Contrast)
#   
#   long_data <- long_data %>%
#     mutate(Gene.Name = fct_reorder(Gene.Name, Log2FC, .fun = mean, .desc = TRUE))
#   
#   ggbarplot(long_data ,
#             "Contrast",
#             "Log2FC",
#             position = position_dodge(0.8),
#             fill="Gene.Name",
#             xlab = "Contrast",
#             #ylab = "LFC Breast Cancer PDX (left) versus Primary CD4 + T Cell (right)",
#             palette = "jco",
#   ) 
# }
# 
# q4<-plotPDX(filtered_genes_colour, genes=c("AR","ESR1","PGR","SMARCD3","FOXA1","BCL11B","FOXP3"))
# q6<-ggarrange(q1, q4, ncol = 2, nrow = 1, labels = c("B", "C", "D"))
# ggarrange(p1, q6, ncol = 1, nrow = 2, labels = c("A"))

plotPDX <- function(DEResultsPDX, 
                    contrast = c("HBC34_IgG_IP_vs_HBC34_GR_IP_log2.fold.change",
                                 "BB3RC31_IgG_IP_vs_BB3RC31_GR_IP_log2.fold.change",
                                 "BB7_IgG_IP_vs_BB7_GR_IP_log2.fold.change"),
                    genes = c("NR3C1", "AR", "ESR1", "PGR", "SMARCD3")) { # Added SMARCD3 to check!
  
  p_adj_cols <- paste0(gsub('log2.fold.change$', '', contrast), "p.adj")
  
  filtered_genes_colour_epi <- DEResultsPDX[DEResultsPDX$Gene.Name %in% genes, ]
  
  long_data <- pivot_longer(
    filtered_genes_colour_epi[, c("Gene.Name", contrast, p_adj_cols)],
    cols = all_of(c(contrast, p_adj_cols)),
    names_to = c("Contrast", ".value"), 
    names_pattern = "(.*)_(log2.fold.change|p.adj)" 
  )
  
  # Clean up names to represent the specific PDX background models
  long_data$Contrast <- gsub("_GR.*$", "", long_data$Contrast)
  long_data$Contrast <- gsub("^.*_", "", long_data$Contrast)
  
  # Force dplyr::rename to avoid package conflicts
  long_data <- long_data %>% 
    dplyr::rename(
      Log2FC = `log2.fold.change`, 
      P.Adj = `p.adj`
    )
  long_data$Log2FC <- -long_data$Log2FC
  
  
  long_data <- long_data %>%
    mutate(Significance = case_when(
      P.Adj < 0.001 ~ "***",
      P.Adj < 0.01  ~ "**",
      P.Adj < 0.05  ~ "*",
      TRUE          ~ ""
    ))
  
  # Group by Gene so they cluster nicely together
  long_data <- long_data %>%
    mutate(Gene.Name = fct_reorder(Gene.Name, Log2FC, .fun = max, .desc = TRUE))
   

  # X-axis becomes the Genes, filled by the individual PDX Models (Contrast)
  ggplot(long_data, aes(x = Gene.Name, y = Log2FC, fill = Contrast)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, color = "black") +
    geom_text(aes(label = Significance, y = Log2FC + (sign(Log2FC) * 0.2)), 
              position = position_dodge(0.8), vjust = 0.5, fontface = "bold") +
    labs(x = "Protein / Gene Name", 
         y = "Log2FC (Reversed)", 
         fill = "PDX Model Group") +
    scale_fill_manual(values = c("HBC34" = "#8da0cb", "BB7" = "#fc8d62", "BB3RC31" = "#66c2a5")) + 
    ggpubr::theme_pubr() +
    theme(legend.position = "top")
} 

ggarrange(p1,p2,p3, ncol=1,nrow=3,common.legend = TRUE, 
          legend = "bottom")

#ggsave(file="PDXvsTCell.svg",width=140, height=160, units="mm")
