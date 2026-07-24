#####
# Figure 1H 
# 
# 23/7/2026
####

library(gprofiler2)
library(tools)
library(stringr)

# Load DIA data
ttDIA <- read.csv("../Figure 1C-K/FragPipe Output/E033_DIA_TimsTOF.csv")

filterDIA<- ttDIA[ ttDIA$GR_vs_IgG_log2.fold.change > 2 &
                    ttDIA$GR_vs_IgG_p.adj < 0.01,]



protein_ids_raw <- unique(filterDIA[[1]])
protein_ids <- unlist(strsplit(as.character(protein_ids_raw), split = "_"))


protein_ids

# Curated GO and pathway terms to highlight, only used on GOST plot 
highlight_ids_curated <- c(
  "GO:0003712", "GO:0003682", "GO:0140297", "GO:0016922", # MF
  "GO:0031981", "GO:0016514", "GO:0000785",               # CC
  "GO:0141193", "GO:0006357", "GO:0006325", "GO:0030518", # BP
  "KEGG:03082", "REAC:R-HSA-4839726", "REAC:R-HSA-9006931"
)

# Run gProfiler enrichment
gost_results <- gost(query = protein_ids,
                      organism = "hsapiens",    
                      correction_method = "fdr",
                      user_threshold = 0.05,
                      evcodes = TRUE,
                      sources = c("GO:MF","GO:CC","GO:BP","KEGG")
                    )

#sources = c("GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC"))


p<-gostplot(gost_results, capped = TRUE, interactive = FALSE)

#pdf("FigureS6_TCell_DIA_Gprofiler.pdf",width=10,height=12,pointsize=8)
publish_gostplot(p,
         highlight_terms = highlight_ids_curated)  # Set TRUE for interactive plot in RStudio viewer
#dev.off()


nuclear_related <- subset(gost_results$result,
                          grepl("nuclear|steroid|chromatin|transcription", term_name, ignore.case = TRUE))

#Bar plot
library(ggplot2)
library(tidyr)


#GG bar plot - Fig 1H
nuclear_related <- subset(gost_results$result,
                          grepl("nuclear receptor|chromatin|transcription", term_name, ignore.case = TRUE))
head(nuclear_related)

nuclear_filtered <- nuclear_related[nuclear_related$term_size >= 100, ]
nuclear_filtered <- nuclear_filtered[nuclear_filtered$term_size <= 1000, ]
nuclear_filtered$gene_ratio <- (nuclear_filtered$intersection_size / nuclear_filtered$query_size) 
nuclear_filtered$background_ratio <- (nuclear_filtered$term_size / 20400 ) 
nuclear_filtered$fold_enrichment <- nuclear_filtered$gene_ratio / nuclear_filtered$background_ratio
nuclear_filtered$log_fold_enrichment <- log2(nuclear_filtered$fold_enrichment)
nuclear_filtered$gene_ratio_text <- paste(round(nuclear_filtered$intersection_size/ nuclear_filtered$query_size *100),"%") 

plot_data<-head(nuclear_filtered [order(-nuclear_filtered$fold_enrichment),], 20)

plot_data$stars <- cut(
  plot_data$p_value,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "")
)

plot_data$term_name<-str_wrap(toTitleCase(as.character(plot_data$term_name)), width=40)


source_labels <- c(
  "GO:BP" = "GO Biological Process",
  "GO:CC" = "GO Cellular Component",
  "GO:MF" = "GO Molecular Function",
  "KEGG"  = "KEGG Pathway",
  "REAC"  = "Reactome Pathway"
)



termBarPlot<-ggplot(plot_data, aes(fill = source, x = reorder(term_name, gene_ratio), y = gene_ratio)) +
  geom_col() +
  geom_text(
    aes(label = stars),
    hjust = -0.2,       
    vjust = 0.75,       
    size = 4.5,
    fontface = "bold"
  ) +
  # geom_text(
  #   aes(label = gene_ratio_text,y=1),
  #   vjust = 0.75,      
  #   size = 4.5,
  #   fontface = "bold"
  # ) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(
    x = NULL, # Removes redundant Y-axis label since term names are self-explanatory
    y = expression("Gene Ratio") # Formats subscript properly
  )  +
  scale_fill_brewer(
    palette = "Set2", 
    labels = source_labels, 
    name = "Terms"
  ) + 
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal" 
  )
  
termBarPlot


#ggsave(termBarPlot,file="1h_barplot.svg", width=120, height=100, units="mm")

