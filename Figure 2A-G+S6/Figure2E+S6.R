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
                      evcodes = TRUE,
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

p

#How many proteins are not nuclear, about 12 in total.
library(dplyr)
chromatin<-gost_results$result %>%
  filter(term_id == "GO:0006357") %>%
  pull(intersection) %>%
  strsplit(",") %>%
  unlist()

nucleoplasm <- gost_results$result %>%
  filter(term_id %in% c("GO:0005634", "GO:0005654","GO:0031981")) %>%
  pull(intersection) %>%
  strsplit(",") %>%
  unlist() %>%
  unique()

cytosol<-gost_results$result %>%
  filter(term_id == "GO:0005829") %>%
  pull(intersection) %>%
  strsplit(",") %>%
  unlist()

library(ggVennDiagram)

venn_list <- list(
  Chromatin = chromatin,
  Nucleoplasm = nucleoplasm,
  Cytosol = cytsol
)

ggVennDiagram(venn_list) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "GO Cellular Component Overlap")

cytosol_only <- setdiff(cytosol, c(chromatin, nucleoplasm))

library(org.Hs.eg.db)
library(AnnotationDbi)

AnnotationDbi::select(org.Hs.eg.db,
                      keys = c("Q13561","Q9NV70","O43303","P43403","Q14204",
                               "P46778","P51149","O75083","P24534","P53621",
                               "P04843","P78371"),
                      columns = "SYMBOL",
                      keytype = "UNIPROT")

###Are cytosol proteins depleted

all_cytosol <- AnnotationDbi::select(org.Hs.eg.db,
                                     keys = "GO:0005829",
                                     columns = "SYMBOL",
                                     keytype = "GO")

n_cytosol_universe <- length(unique(all_cytosol$SYMBOL))
n_universe <- 20400
n_detected <- length(unique(unlist(venn_list)))
n_cytosol_detected <- length(cytosol_only)

phyper(n_cytosol_detected,
       n_cytosol_universe,
       n_universe - n_cytosol_universe,
       n_detected,
       lower.tail = TRUE)  # lower.tail = TRUE tests for depletio
#[1] 1.903274e-05