#D714 Method Comparison / Low cell number data
#Sample 1 - RU486 - Thermo Kit - Unused for this paper
#Sample 2 - RU486 - Overnight  - Unused for this paper   
#Sample 3 - RU486 - 2 hours    - Unused for this paper   
#Sample 4 - Dex   - Thermo Kit - Unused for this paper   
#Sample 5 - Dex   - Overnight  
#Sample 6 - Dex   - 2 hours    - Unused for this paper   
#Sample 7 - DMSO  - Thermo Kit - Unused for this paper   
#Sample 8 - DMSO  - Overnight  
#Sample 9 - DMSO  - 2 hours    - Unused for this paper   

#Exported from D714.sf3 
#Protein Threshold 1% FDR,
#Min # of peptides 2,
#Peptide 1% FDR
#Selected Total Spectral Count 
#Exported to XLS using toolbar

df<-read.delim('Scafold/Samples Report With Clusters for D714.tsv')
dfSamples<-df[,c(1:9,15,18)]
#colnames(dfSamples)[c(10,11)]<-c("Dex","DMSO") # 2 hours
colnames(dfSamples)[c(9,10)]<-c("Dex","DMSO") #Overnight
dfNoClusters<-dfSamples[grep("\\.", dfSamples$X., invert=TRUE ),]


#Dex Results
TFs<-dfNoClusters[dfNoClusters$Dex > 0,]$Alternate.ID

length(TFs)
#overnight : 1635 
#2 hours: 623

BioGrid<-read.delim("../Figure S2+3/BIOGRID-GENE-109165-4.4.245.DOWNLOADS/BIOGRID-GENE-109165-4.4.245.tab3.txt")

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
#787 protein retrieved from Biogrid (Same as Figure S2+3)

length(TFs[TFs %in% BGTFs])
#2 hours: 60 TFs from RIME in Biogrid
#Overnight: 183 TFs from RIME in Biogrid


#StringDB - From Figure S2+3
#Low confidence < 0.15
#Max interactions 1st shell < 250
#Max interactions 2nd shell < 50
#active interaction sources - only experiments.
#Downloaded 30/5/2025
String<-read.delim("../Figure S2+3/StringDB/string_interactions_short.tsv")

length(unique(String$X.node1,String$X.node2))
#499 proteins retrived from String database

StringTFs<-unique(String$X.node1)
length(TFs[TFs %in% StringTFs])
#2 hours: 57 TFs from RIME in String
#Overnight: 144 TFs from RIME in String 

##Hypergeometric tests 
#BioGrid
phyper(length(TFs[TFs %in% BGTFs]) - 1,
       length(unique(BGTFs)), 
       20400 - length(unique(BGTFs)) ,
       #Estimate of total number of proteins in Uniprot
       length(TFs),
       lower.tail = FALSE)
#2 hours: 7.068341e-11
#Overnight: 2.116724e-41

#String
phyper(length(TFs[TFs %in% StringTFs]) - 1,
       length(unique(StringTFs)), 
       20400 - length(unique(StringTFs)) ,
       #Estimate of total number of proteins in Uniprot
       length(TFs),
       lower.tail = FALSE)
#2 hours 7.840464e-18
#OVernight: 3.739006e-44

#GProfiler for Dex
library(gprofiler2)

gost_results <- gost(
  query = TFs,
  organism = "hsapiens",  # Homo sapiens
  sources = c("GO:BP", "GO:CC", "GO:MF", "REAC", "KEGG"),
  user_threshold = 0.05,
  correction_method = "fdr"
)



p<-gostplot(gost_results, capped = TRUE, interactive = FALSE)
highlight_ids_curated<-c(
  "GO:0003712", "GO:0003682","GO:0140297","GO:0016922", #MF
  "GO:0031981", "GO:0016514","GO:0000785", #CC
  "GO:0141193","GO:0006357","GO:0006325","GO:0030518", #BP
  "KEGG:03082","REAC:R-HSA-4839726", "REAC:R-HSA-9006931"
)      


pdf("FigureS4_Gprofiler.pdf",width=10,height=12,pointsize=8)
publish_gostplot(p, highlight_terms = highlight_ids_curated)
dev.off()

#GProfiler for Dex
#Dex Results
DMSO_TFs<-dfNoClusters[dfNoClusters$DMSO > 0,]$Alternate.ID

gost_results_DMSO <- gost(
  query = DMSO_TFs,
  organism = "hsapiens",  # Homo sapiens
  sources = c("GO:BP", "GO:CC", "GO:MF", "REAC", "KEGG"),
  user_threshold = 0.05,
  correction_method = "fdr"
)

p_DMSO<-gostplot(gost_results_DMSO, capped = TRUE, interactive = FALSE)

pdf("FigureS5_Gprofiler_DMSO.pdf",width=10,height=12,pointsize=8)
publish_gostplot(p_DMSO, highlight_terms = highlight_ids_curated)
dev.off()

##Compare results

highlight_ids_curated_compared<-c(
  "GO:0003712", "GO:0003682","GO:0140297","GO:0016922", #MF
  "GO:0031981", "GO:0016514","GO:0000785", #CC
  "GO:0141193","GO:0006357","GO:0006325","GO:0030518", #BP
  "KEGG:03082","REAC:R-HSA-4839726", "REAC:R-HSA-9006931",
  "GO:0031072", "GO:0051082", "GO:0140662",
  "GO:0101031",
  "REAC:R-HSA-3371497"
)      


df_Dex<- cbind(gost_results$result[gost_results$result$term_id %in% 
                               highlight_ids_curated_compared,],
                condition="Dex")
df_DMSO<- cbind(gost_results_DMSO$result[gost_results_DMSO$result$term_id %in%
                                highlight_ids_curated_compared,],
                condition="DMSO")


combined_results <- rbind(df_Dex, df_DMSO)
combined_results$log10p <- -log10(combined_results$p_value)

combined_results$scaledPValue=1
for (termID in combined_results[combined_results$condition == 'DMSO',]$term_id ) {
  combined_results[combined_results$condition == 'DMSO' & 
                     combined_results$term_id == termID,]$scaledPValue <-
    combined_results[combined_results$condition == 'DMSO' & 
                       combined_results$term_id == termID,]$log10p /
    combined_results[combined_results$condition == 'Dex' & 
                       combined_results$term_id == termID,]$log10p
}


library(ggplot2)

pdf("FigureS6_term_comparison.pdf",width=10,height=12,pointsize=8)
ggplot(combined_results, aes(x = reorder(term_name, log10p), y = log10p, fill =condition)) +
  geom_bar(stat = "identity", position = position_dodge(preserve='single', width = -0.7), width = 0.6) +
  coord_flip() +
  scale_fill_manual(values = c("Dex" = "#E69F00", "DMSO" = "#56B4E9")) +
  labs(
    x = NULL, y = expression(-log[10]~"p-value"),
    title = "Comparison of Key Terms (Dex vs DMSO)"
  ) +
  theme_minimal(base_size = 12)



ggplot(combined_results, aes(x = reorder(term_name, log10p), y = scaledPValue, fill =condition)) +
  geom_bar(stat = "identity", position = position_dodge(preserve='single', width = -0.7), width = 0.6) +
  coord_flip() +
  scale_fill_manual(values = c("Dex" = "#E69F00", "DMSO" = "#56B4E9")) +
  labs(
    x = NULL, y = expression( Scaled -log[10]~"p-value"),
    title = "Comparison of Key Terms (Dex vs DMSO)"
  ) +
  theme_minimal(base_size = 12)

dev.off()
