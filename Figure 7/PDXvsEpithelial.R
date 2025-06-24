library(ggplot2)
library(ggpubr)


DEResults <- read.csv('../Figure 4+S7/FragPipe Output/DE_results.csv')
contrasts <- colnames(DEResults)


#Set up conditions and LFC contrasts
conditions <- unique(
  read.csv(
    '../Figure 4+S7/FragPipe Input/E733_All_FP-A_experimental_annotation_pool_removed.tsv',
    sep = "\t"
  )$condition
)
conditionsNoPool <- conditions[!grepl('Pool', conditions)]
conditionsGRonly <- conditionsNoPool[grepl('_GR$', conditionsNoPool)]
tissues <- gsub('.{3}$', '', conditionsGRonly)
tissues[tissues == '231'] <- "X231" #Safe naming
contrastsLFC <- paste0(tissues, "_IgG_vs_", tissues, "_GR_log2.fold.change")


"MCF7_GR_vs_X231_GR_log2.fold.change" %in%   contrasts

##NR3C1 #Downloaded Currated Data from BioGrid for NR3C1 30/5/2025
BioGrid <- read.delim(
  "../Figure 1C+S1+S2/BIOGRID-GENE-109165-4.4.245.DOWNLOADS/BIOGRID-GENE-109165-4.4.245.tab3.txt"
)

BGa <- BioGrid[(
  BioGrid$Experimental.System == "Affinity Capture-MS" |
    BioGrid$Experimental.System ==  "Proximity Label-MS"
) &
  BioGrid$Organism.Name.Interactor.A == "Homo sapiens" &
  BioGrid$Organism.Name.Interactor.B == "Homo sapiens", ]$Official.Symbol.Interactor.A

BGb <- BioGrid[(
  BioGrid$Experimental.System == "Affinity Capture-MS" |
    BioGrid$Experimental.System ==  "Proximity Label-MS"
) &
  BioGrid$Organism.Name.Interactor.A == "Homo sapiens" &
  BioGrid$Organism.Name.Interactor.B == "Homo sapiens", ]$Official.Symbol.Interactor.B

grComplex <- unique(c(BGa, BGb))

##

colVector <- rep("grey30", length(DEResults$'Gene.Name'))
colVector[na.omit(match(grComplex,DEResults$'Gene.Name'))] <- "steelblue"
names(colVector)[colVector == 'steelblue'] <- 'GR'
names(colVector)[colVector == 'grey30'] <- ''

#Plot grey points first
DEResultsReorder <- DEResults[order(!colVector=="grey30"),]
colVectorReorder <- colVector[order(!colVector=="grey30")]

library(EnhancedVolcano)

GRbiogrid<- DEResults$'Gene.Name'[colVector == 'steelblue']



CoreGRComplex<-c("NR3C1") #GRbiogrid

DEResultsReorder$Negative_BB3RC31_IgG_IP_vs_BB3RC31_GR_IP_log2.fold.change <- -DEResultsReorder$BB3RC31_IgG_IP_vs_BB3RC31_GR_IP_log2.fold.change

p1<-EnhancedVolcano(DEResultsReorder,
                lab = DEResultsReorder$'Gene.Name',
                x =  'Negative_BB3RC31_IgG_IP_vs_BB3RC31_GR_IP_log2.fold.change', #I prefer IgG on the left
                y = 'BB3RC31_IgG_IP_vs_BB3RC31_GR_IP_p.adj',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                title = NULL,
                subtitle = NULL,
                caption = "BB3RC31 IgG (left) vs BB3RC31 GR (right)",
                drawConnectors = TRUE,
                widthConnectors = 0.4,
                boxedLabels = TRUE,
                pointSize = 1.0,
                labSize = 2.5,
                max.overlaps = 50,
                colCustom  = colVectorReorder,
                selectLab = CoreGRComplex,
                colAlpha = 1
) + theme_classic()+ theme(legend.position = "none")

DEResultsReorder$Negative_BB7_IgG_IP_vs_BB7_GR_IP_log2.fold.change <- -DEResultsReorder$BB7_IgG_IP_vs_BB7_GR_IP_log2.fold.change

p2<-EnhancedVolcano(DEResultsReorder,
                    lab = DEResultsReorder$'Gene.Name',
                    x =  'Negative_BB7_IgG_IP_vs_BB7_GR_IP_log2.fold.change', #I prefer IgG on the left
                    y = 'BB7_IgG_IP_vs_BB7_GR_IP_p.adj',
                    pCutoff = 0.05,
                    FCcutoff = 1.0,
                    title = NULL,
                    subtitle = NULL,
                    caption = "BB7 IgG (left) vs BB7 GR (right)",
                    drawConnectors = TRUE,
                    widthConnectors = 0.4,
                    boxedLabels = TRUE,
                    pointSize = 1.0,
                    labSize = 2.5,
                    max.overlaps = 50,
                    colCustom  = colVectorReorder,
                    selectLab = CoreGRComplex,
                    colAlpha = 1
) + theme_classic()+ theme(legend.position = "none")

DEResultsReorder$Negative_HBC34_IgG_IP_vs_HBC34_GR_IP_log2.fold.change <- -DEResultsReorder$HBC34_IgG_IP_vs_HBC34_GR_IP_log2.fold.change

p3<-EnhancedVolcano(DEResultsReorder,
                    lab = DEResultsReorder$'Gene.Name',
                    x =  'Negative_HBC34_IgG_IP_vs_HBC34_GR_IP_log2.fold.change', #I prefer IgG on the left
                    y = 'HBC34_IgG_IP_vs_HBC34_GR_IP_p.adj',
                    pCutoff = 0.05,
                    FCcutoff = 1.0,
                    title = NULL,
                    subtitle = NULL,
                    caption = "HBC34 IgG (left) vs HBC34 GR (right)",
                    drawConnectors = TRUE,
                    widthConnectors = 0.4,
                    boxedLabels = TRUE,
                    pointSize = 1.0,
                    labSize = 2.5,
                    max.overlaps = 50,
                    colCustom  = colVectorReorder,
                    selectLab = CoreGRComplex,
                    colAlpha = 1
) + theme_classic()+ theme(legend.position = "none")


row1<-ggarrange(p1,p2,p3,ncol=3)
row1

#Violin plot
library(tidyverse)
library(ggpubr)

lfc_cols <- c("Negative_BB3RC31_IgG_IP_vs_BB3RC31_GR_IP_log2.fold.change",
              "Negative_BB7_IgG_IP_vs_BB7_GR_IP_log2.fold.change",
              "Negative_HBC34_IgG_IP_vs_HBC34_GR_IP_log2.fold.change")

sample_names <- c("BB3RC31", "BB7", "HBC34")


lfc_long <- DEResultsReorder %>%
  dplyr::select(Gene.Name, all_of(lfc_cols)) %>%
  pivot_longer(cols = -Gene.Name,
               names_to = "Comparison",
               values_to = "LFC") %>%
  mutate(Sample = case_when(
    str_detect(Comparison, "BB3RC31") ~ "BB3RC31",
    str_detect(Comparison, "BB7") ~ "BB7",
    str_detect(Comparison, "HBC34") ~ "HBC34"
  ),
  Group = ifelse(Gene.Name %in% grComplex, "GR Complex", "Other")
  )


violin_plot <- ggviolin(
  lfc_long[lfc_long$Group=="GR Complex",],
  x = "Sample",
  y = "LFC",
  fill = "steelblue",
  #palette = c("GR Complex" = "steelblue", "Other" = "grey80"),
  add.params = list(width = 0.1),
  trim = FALSE
)  +
  theme_classic() +
  labs(y = "log2 Fold Change (GR vs IgG)", x = NULL) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  labs(
    y = expression(Log[2]~Fold~Change),
    x = "PDX"
  ) +
  # Add dotted zero line
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,           # solid diamond shape
    size = 2.5,
    color = "black",
    position = position_dodge(width = 0.8)
  ) 

violin_plot

### ClusterProfiler
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(reactome.db)

enrichTerms<-function(pdx="BB7"){
    #Terms from Figure
    terms_of_interest <- c(
      "GO:0003712",     # transcription co-regulator activity
      "GO:0016922",     # nuclear receptor binding
    #  "R-HSA-9006931",  # signalling by nuclear receptors (react to see blow)
    #  "GO:0141193",     # nuclear receptor-mediated signaling pathway not in org.Hs.eg.db::org.Hs.egGO2ALLEGS
      "GO:0030518",      # steroid hormone signaling
      "GO:0042921" #nuclear receptor-mediated glucocorticoid signaling pathway
    )
    
    go2genes <- AnnotationDbi::mget(terms_of_interest, org.Hs.eg.db::org.Hs.egGO2ALLEGS, ifnotfound = NA)
    term2gene_go <- stack(go2genes)
    colnames(term2gene_go) <- c("gene", "term")
    term2gene_go <- term2gene_go[, c("term", "gene")]
    term2gene_go <- dplyr::filter(term2gene_go, !is.na(gene))
    
    
    
    xx <- as.list(reactomePATHID2EXTID)
    term2gene_React<-xx[["R-HSA-9006931"]]
    
    term2gene<-rbind(term2gene_go,data.frame(gene=term2gene_React,term="R-HSA-9006931"))
    
    lfc_col <- paste0("Negative_", pdx, "_IgG_IP_vs_", pdx, "_GR_IP_log2.fold.change")
    pval_col <- paste0(pdx, "_IgG_IP_vs_", pdx, "_GR_IP_p.val")
    
    
    filtered_proteins <- DEResultsReorder$Gene.Name[
                            DEResultsReorder[[lfc_col]] > 2 &
                            DEResultsReorder[[pval_col]] < 0.01
                            ]
          
    # Convert gene symbols to ENTREZ IDs
    gene_entrez <- bitr(filtered_proteins, fromType = "SYMBOL",
                        toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    
    # Convert gene symbols to ENTREZ IDs
    gene_entrez_bg <- bitr(DEResultsReorder$Gene.Name, fromType = "SYMBOL",
                        toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    
    
    
    head(gene_entrez)
    head(term2gene)
    gene_entrez$SYMBOL[gene_entrez$ENTREZID %in%  term2gene$gene]
    
    
    background_genes <- keys(org.Hs.eg.db, keytype = "ENTREZID")
    
    #https://github.com/YuLab-SMU/clusterProfiler/issues/283
    #Without this the background is the geneset, which means that there is 100% 
    #of finding thing within the geneset
    options(enrichment_force_universe=TRUE)
    
    es <- enricher(
      gene = gene_entrez$ENTREZID,
      TERM2GENE = term2gene,
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      universe = gene_entrez_bg$ENTREZID, #background_genes,
      minGSSize=0
    )
    
    return(as.data.frame(es))
   
}

df<-rbind(cbind(enrichTerms("BB3RC31")[c('ID','p.adjust')],pdx="BB3RC31"),
          cbind(enrichTerms("BB7")[c('ID','p.adjust')],pdx="BB7"),
          cbind(enrichTerms("HBC34")[c('ID','p.adjust')],pdx="HBC34")
          )

id_descriptions <- c(
  "GO:0003712" = "Transcription coregulator\nactivity (GO:0003712)",
  "GO:0030518" = "Nuclear receptor-mediated\nsteroid hormone signaling (GO:0030518)",
  "GO:0042921" = "Nuclear receptor-mediated\nglucocorticoid signaling (GO:0042921)",
  "GO:0016922" = "Nuclear receptor binding (GO:0016922)",
  "R-HSA-9006931" = "Signalling by Nuclear\nReceptors (R-HSA-9006931)"
)

df$negLogP <- -log10(df$p.adjust)
df$Description <- id_descriptions[df$ID]
rownames(df) <- NULL

barplot<-ggbarplot(df, x="pdx",y="negLogP",
          fill="Description",
          position= position_dodge(-0.75),
          palette = get_palette(c("#DD5577", "#00B8E7", "#FCAE07"), 5) )+
          rotate() +
          theme(legend.position = "right",
                legend.key.height = unit(1.0, "cm"),
                legend.spacing.y = unit(0.5, "cm"),
                legend.text = element_text(lineheight = 0.9))+ 
          geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
          labs(
            x = "PDX",
            y = expression(-log[10]~"(adjusted p-value)"),
            fill = "Term Description"
          )
          
row2<-ggarrange(violin_plot,barplot, widths=c(1,3),ncol=2)


ggarrange(row1,row2,nrow=2)


#Final blot is combined PDX vs IgG and PDX vs epithelial..