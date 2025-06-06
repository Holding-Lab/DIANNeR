library(gprofiler2)

# Load DIA data
ttDIA <- read.csv("FragPipe Output/E033_DIA_TimsTOF.csv")

filterDIA<- ttDIA[ ttDIA$GR_vs_IgG_log2.fold.change > 2 &
                    ttDIA$GR_vs_IgG_p.adj < 0.01,]



protein_ids_raw <- unique(filterDIA[[1]])
protein_ids <- unlist(strsplit(as.character(protein_ids_raw), split = "_"))


protein_ids

# Curated GO and pathway terms to highlight
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
                      sources = c("GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC"))


p<-gostplot(gost_results, capped = TRUE, interactive = FALSE)


pdf("FigureS6_TCell_DIA_Gprofiler.pdf",width=10,height=12,pointsize=8)
publish_gostplot(p,
         highlight_terms = highlight_ids_curated)  # Set TRUE for interactive plot in RStudio viewer
dev.off()


nuclear_related <- subset(gost_results$result,
                          grepl("nuclear|steroid|chromatin|transcription", term_name, ignore.case = TRUE))

#Bar plot
library(ggplot2)
library(tidyr)


#GG bar plot - Fig 2E
nuclear_related <- subset(gost_results$result,
                          grepl("nuclear|steroid|chromatin|transcription", term_name, ignore.case = TRUE))


termBarPlot<-ggplot(head(nuclear_related[order(nuclear_related$p_value),], 20), aes(x = reorder(term_name, -p_value), y = -log10(p_value))) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top enrichend Nuclear, Chromatin and Transcription Terms",
       x = "GO/Reactome/KEGG Pathway Term",
       y = "-log10(p-value)") +
  theme_minimal()
ggsave("Figure2E_Terms.svg",device="svg", termBarPlot,units="mm",dpi=300,width=175, height=75)
