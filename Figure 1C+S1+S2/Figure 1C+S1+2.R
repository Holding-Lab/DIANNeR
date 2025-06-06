#D676 Antibody validation
#Exported from D676.sf3 
#Protein Threshold 1% FDR,
#Min # of peptides 2,
#Peptide 1% FDR
#Exported Total Spectral Count 
df<-read.delim('Scafold/Samples Report With Clusters for D676.tsv')

dfNoClusters<-df[grep("\\.", df$X., invert=TRUE ),]


TFs<-dfNoClusters[dfNoClusters$D676_AT..AT. > 0,]$Alternate.ID

# Code to remove any protein from with IgG sample, however, this removes a lot
#  of proteins that are in BioGrid
#    [
#    !(dfNoClusters[dfNoClusters$D676_AT..AT. > 0,]$Alternate.ID %in% 
#      dfNoClusters[dfNoClusters$D676_IgG..IgG. > 0,]$Alternate.ID)
#        ]

length(TFs)
#991  TFs


#BiogGrid
#Downloaded Currated Data from BioGrid for NR3C1 30/5/2025
BioGrid<-read.delim("BIOGRID-GENE-109165-4.4.245.DOWNLOADS/BIOGRID-GENE-109165-4.4.245.tab3.txt")

BGa<-BioGrid[  (BioGrid$Experimental.System == "Affinity Capture-MS" | 
             BioGrid$Experimental.System ==  "Proximity Label-MS") & 
            BioGrid$Organism.Name.Interactor.A == "Homo sapiens" &
            BioGrid$Organism.Name.Interactor.B == "Homo sapiens",
            ]$Official.Symbol.Interactor.A

BGb<-BioGrid[  (BioGrid$Experimental.System == "Affinity Capture-MS" | 
            BioGrid$Experimental.System ==  "Proximity Label-MS") & 
            BioGrid$Organism.Name.Interactor.A == "Homo sapiens" &
            BioGrid$Organism.Name.Interactor.B == "Homo sapiens",
            ]$Official.Symbol.Interactor.B

BGTFs<-unique(c(BGa,BGb))
length(BGTFs)
#787 protein retrieved from Biogrid

length(TFs[TFs %in% BGTFs])
#106 TFs from RIME in Biogrid

#StringDB
#Low confidence < 0.15
#Max interactions 1st shell < 250
#Max interactions 2nd shell < 50
#active interaction sources - only experiments.
#Downloaded 30/5/2025
String<-read.delim("StringDB/string_interactions_short.tsv")

length(unique(String$X.node1,String$X.node2))
#499 proteins retrived from String database

StringTFs<-unique(String$X.node1)
length(TFs[TFs %in% StringTFs])
#99 TFs from RIME in String

##Hypergeometric tests 
#BioGrid
phyper(length(TFs[TFs %in% BGTFs]) - 1,
       length(unique(BGTFs)), 
       20400 - length(unique(BGTFs)) ,
       #Estimate of total number of proteins in Uniprot
       length(TFs),
       lower.tail = FALSE)
#[1] 6.070497e-22

#String
phyper(length(TFs[TFs %in% StringTFs]) - 1,
               length(unique(StringTFs)), 
               20400 - length(unique(StringTFs)) ,
                          #Estimate of total number of proteins in Uniprot
               length(TFs),
               lower.tail = FALSE)
#[1] 4.753852e-34

##StringDB images
library(STRINGdb)
string_db <- STRINGdb$new(version = "11.5", species = 9606,
                          score_threshold = 400, input_directory = "")

SpecificTFs<-dfNoClusters[dfNoClusters$D676_AT..AT. > 0,]$Alternate.ID[
    !(dfNoClusters[dfNoClusters$D676_AT..AT. > 0,]$Alternate.ID %in% 
      dfNoClusters[dfNoClusters$D676_IgG..IgG. > 0,]$Alternate.ID)
        ]
length(SpecificTFs)
#169

genes_df <- data.frame(gene = SpecificTFs)
mapped_genes <- string_db$map(genes_df, "gene", removeUnmappedRows = TRUE)
png("FigureS1_StringNetwork.png",width=180, height=180, units="mm",res = 300)
string_db$plot_network(mapped_genes$STRING_id)
dev.off()

#GProfiler on SpecificTFs
library(gprofiler2)

gost_results <- gost(
  query = SpecificTFs,
  organism = "hsapiens",  # Homo sapiens
  sources = c("GO:BP", "GO:CC", "GO:MF", "REAC", "KEGG"),
  user_threshold = 0.05,
  correction_method = "fdr"
)

#Extract nuclear related
nuclear_related <- subset(gost_results$result,
                          grepl("nuclear|steroid|chromatin|transcription", term_name, ignore.case = TRUE))

#Bar plot
library(ggplot2)
library(tidyr)

#GG bar plot - Fig 1C
termBarPlot<-ggplot(head(nuclear_related[order(nuclear_related$p_value),], 20), aes(x = reorder(term_name, -p_value), y = -log10(p_value))) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top enriched nuclear/chromatin/transcription terms",
       x = "GO/Reactome/KEGG Pathway Term",
       y = "-log10(p-value)") +
  theme_minimal()
ggsave("Figure1C_Terms.pdf",termBarPlot,units="mm",dpi=300,width=175, height=90)

#Gost plot 
p<-gostplot(gost_results, capped = TRUE, interactive = FALSE)
p
term_names<-grep("transcript|chromatin|SWI/SNF|nuclear|steroid", gost_results$result$term_name, ignore.case = TRUE, value = TRUE)
highlighted_terms <- gost_results$result[gost_results$result$term_name %in% term_names, ]
highlight_ids <- highlighted_terms$term_id
#Completed dominated by GO:BP and so big R hangs if we just take a top cut so have to be selective.
highlight_ids_curated<-c(
      "GO:0003712", "GO:0003682","GO:0140297","GO:0016922", #MF
      "GO:0031981", "GO:0016514","GO:0000785", #CC
      "GO:0141193","GO:0006357","GO:0006325","GO:0030518", #BP
      "KEGG:03082","REAC:R-HSA-4839726", "REAC:R-HSA-9006931"
      )      


pdf("FigureS2_Gprofiler.pdf",width=10,height=12,pointsize=8)
p2<-publish_gostplot(p, highlight_terms = highlight_ids_curated)
dev.off()

