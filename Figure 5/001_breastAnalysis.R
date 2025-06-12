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

library(EnhancedVolcano)

#BioGrid ESR1 downloaded 1/05/2025
ESR1BioGrid <- read.csv("esr1-biogrid/BIOGRID-GENE-108403-4.4.245.tab3.txt", sep =
                          "\t")

ESR1BioGrid[(ESR1BioGrid$Experimental.System== "Affinity Capture-MS" |
              ESR1BioGrid$Experimental.System == "Proximity Label-MS") &
              ESR1BioGrid$Organism.ID.Interactor.A == "9606" &
              ESR1BioGrid$Organism.ID.Interactor.B == "9606", c(8, 9)]
ESR1interactions <- ESR1BioGrid[(ESR1BioGrid$Experimental.System == "Affinity Capture-MS" |
                                   ESR1BioGrid$Experimental.System. == "Proximity Label-MS") &
                                  ESR1BioGrid$Organism.ID.Interactor.A == "9606" &
                                  ESR1BioGrid$Organism.ID.Interactor.B == "9606", c(8, 9)]
ESR1interactor <- unique(
  c(
    ESR1interactions[ESR1interactions$Official.Symbol.Interactor.A == 'ESR1', ]$Official.Symbol.Interactor.B,
    ESR1interactions[ESR1interactions$Official.Symbol.Interactor.B ==
                       'ESR1', ]$Official.Symbol.Interactor.A
  )
)

esr1Complex <- ESR1interactor
CoreEsr1Complex<-c("ESR1", "FOXA1", "GATA3", "GREB1", "NCOA3", 'GRHL2','GRHL1',
                'PGR','FOS','JUN', 'NCOA1','MED1','NRIP1','SRC','ZEB1','ZEB2')  #includes ZEB1/2 though these not 'core'
#esr1Complex <- ESR1interactorsGenes$preferred_name

#NR3C1 #Downloaded Currated Data from BioGrid for NR3C1 30/5/2025
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

esr1ComplexUnique<-esr1Complex[!esr1Complex  %in% grComplex]
grComplexUnique<-grComplex[!grComplex %in% esr1Complex]
bothComplexes<-grComplex[grComplex %in% esr1Complex]

esr1Labs<-DEResults$'Gene.Name'
esr1Labs[!(DEResults$'Gene.Name' %in% CoreEsr1Complex)]<-""

###


colVector <- rep("grey30", length(DEResults$'Gene.Name'))
colVector[na.omit(match(esr1ComplexUnique,DEResults$'Gene.Name'))] <- "hotpink"
colVector[na.omit(match(grComplexUnique,DEResults$'Gene.Name'))] <- "steelblue"
colVector[na.omit(match(bothComplexes,DEResults$'Gene.Name'))] <- "seagreen"
names(colVector)[colVector == 'grey30'] <- 'Neither'
names(colVector)[colVector == 'steelblue'] <- 'GR'
names(colVector)[colVector == 'seagreen'] <- 'Both'
names(colVector)[colVector == 'hotpink'] <- 'ERa'

p1<-EnhancedVolcano(DEResults,
                lab = DEResults$'Gene.Name',
                x = 'MCF7_GR_vs_X231_GR_log2.fold.change',
                y = 'MCF7_GR_vs_X231_GR_p.val',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                title = NULL,
                subtitle = NULL,
                caption = "MDA-MB-231 (left) vs MCF7 (right)",
                drawConnectors = TRUE,
                widthConnectors = 0.4,
                boxedLabels = TRUE,
                pointSize = 1.0,
                labSize = 2.5,
                max.overlaps = 50,
                colCustom = colVector,
                selectLab = esr1Labs,
                #selectLab = c(""),
                colAlpha = 1
) + theme_classic()+ theme(legend.position = "none")

DEResults$MCF7_GR_vs_X231_GR_log2.fold.change
DEResults$'Gene.Name'

violinPlot<-data.frame(LFC=DEResults$MCF7_GR_vs_X231_GR_log2.fold.change,BioGrid=as.factor(names(colVector)))
violinPlot<-violinPlot[!violinPlot$BioGrid=='Both',]


library(ggpubr)
p2<-ggplot(violinPlot, aes(y=LFC, x=BioGrid,fill=BioGrid)) +
      geom_violin() + 
      theme_classic() + 
      scale_fill_manual(values=c( "hotpink","steelblue","grey30")) +
      labs(x="BioGrid Interactors", y= expression(Log[2]~fold~change) )

my_comparisons <- list( c("ERa", "Neither"), c("GR", "Neither") )
  
p2<-p2 +  stat_compare_means(comparisons = my_comparisons)

ptop<-ggarrange(p1,p2,ncol=2)

#####
##MCF7 vs EPIs

p3<-EnhancedVolcano(DEResults,
                lab = DEResults$'Gene.Name',
                x = 'MCF7_GR_vs_PrimaryBreastEpis_GR_log2.fold.change',
                y = 'MCF7_GR_vs_PrimaryBreastEpis_GR_p.val',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                title = NULL,
                subtitle = NULL,
                caption = "Primary Breast Epithelium (left) - MCF7 (right)",
                drawConnectors = TRUE,
                widthConnectors = 0.4,
                boxedLabels = TRUE,
                pointSize = 1.0,
                labSize = 2.5,
                max.overlaps = 50,
                colCustom = colVector,
                selectLab = esr1Labs,
                #selectLab = c(""),
                colAlpha = 1
) + theme_classic()+ theme(legend.position = "none")



violinPlot2<-data.frame(LFC=DEResults$MCF7_GR_vs_PrimaryBreastEpis_GR_log2.fold.change,BioGrid=as.factor(names(colVector)))
violinPlot2<-violinPlot2[!violinPlot2$BioGrid=='Both',]

library(ggpubr)
p4<-ggplot(violinPlot2, aes(y=LFC, x=BioGrid,fill=BioGrid)) +
  geom_violin() + 
  theme_classic() + 
  scale_fill_manual(values=c("hotpink","steelblue","grey30")) +
  labs(x="BioGrid Interactors", y= expression(Log[2]~fold~change) )

my_comparisons <- list( c("ERa", "Neither"), c("GR", "Neither") )

p4<-p4 +  stat_compare_means(comparisons = my_comparisons)


pbottom<-ggarrange(p3,p4,ncol=2)




##Epis vs 231
'%ni%' <- Negate('%in%')

HOX_genes <- grep("HOX", DEResults$Gene.Name) %ni% grep(";", DEResults$Gene.Name) 
DEResultsHOX <- DEResults[c(setdiff(1:nrow(DEResults), HOX_genes), HOX_genes), ]


Epis231Labels<-c("ESR1", "FOXA1", "GATA3", "GREB1",'GRHL2','GRHL1','ZEB1','ZEB2')  #includes ZEB1/2 though these not 'core'
HOX<-DEResults$Gene.Name[grep("HOX", DEResults$Gene.Name)]

colVector <- rep('#D0D0D010', length(DEResultsHOX$'Gene.Name'))
colVector[na.omit(match(DEResultsHOX$Gene.Name[grepl("HOX", DEResultsHOX$Gene.Name)], DEResultsHOX$Gene.Name))] <- '#FF0000FF'
names(colVector)[colVector == '#D0D0D05'] <- 'Neither'
names(colVector)[colVector == '#FF0000FF'] <- 'HOX'


p5<-EnhancedVolcano(DEResultsHOX,
                    lab = DEResultsHOX$'Gene.Name',
                    x = 'X231_GR_vs_PrimaryBreastEpis_GR_log2.fold.change',
                    y = 'X231_GR_vs_PrimaryBreastEpis_GR_p.val',
                    pCutoff = 0.05,
                    FCcutoff = 1.0,
                    title = NULL,
                    subtitle = NULL,
                    caption = "Primary Breast Epithelium (Left) vs MDA-MB-231 (right)",
                    drawConnectors = TRUE,
                    widthConnectors = 0.4,
                    boxedLabels = TRUE,
                    pointSize = 1.0,
                    labSize = 2.5,
                    max.overlaps = 50,
                    colCustom = colVector,
                    selectLab = c(Epis231Labels,HOX),
                    #selectLab = c("")
) + theme_classic()+ theme(legend.position = "none")


p6<-EnhancedVolcano(DEResultsHOX,
                    lab = DEResultsHOX$'Gene.Name',
                    x = 'MCF7_GR_vs_PrimaryBreastEpis_GR_log2.fold.change',
                    y = 'MCF7_GR_vs_PrimaryBreastEpis_GR_p.val',
                    pCutoff = 0.05,
                    FCcutoff = 1.0,
                    title = NULL,
                    subtitle = NULL,
                    caption = "Primary Breast Epithelium (Left) vs MCF7 (right)",
                    drawConnectors = TRUE,
                    widthConnectors = 0.4,
                    boxedLabels = TRUE,
                    pointSize = 1.0,
                    labSize = 2.5,
                    max.overlaps = 50,
                    colCustom = colVector,
                    selectLab = c(HOX),
                    #selectLab = c("")
) + theme_classic()+ theme(legend.position = "none")

#Flip to match previous 2 plots.
DEResultsHOX$PrimaryBreastEpis_GR_vs_MCF10A_GR_minuslog2.fold.change <- -DEResultsHOX$PrimaryBreastEpis_GR_vs_MCF10A_GR_log2.fold.change
p7<-EnhancedVolcano(DEResultsHOX,
                    lab = DEResultsHOX$'Gene.Name',
                    x = 'PrimaryBreastEpis_GR_vs_MCF10A_GR_minuslog2.fold.change',
                    y = 'PrimaryBreastEpis_GR_vs_MCF10A_GR_p.val',
                    pCutoff = 0.05,
                    FCcutoff = 1.0,
                    title = NULL,
                    subtitle = NULL,
                    caption = "Primary Breast Epithelium (Left) vs MCF10A (right)",
                    drawConnectors = TRUE,
                    widthConnectors = 0.4,
                    boxedLabels = TRUE,
                    pointSize = 1.0,
                    labSize = 2.5,
                    max.overlaps = 50,
                    colCustom = colVector,
                    selectLab = c(HOX),
                    #selectLab = c("")
) + theme_classic()+ theme(legend.position = "none")




#HOXA5 HOXA5 determines cell fate transition and impedes tumor initiation and progression in breast cancer through regulation of E-cadherin and CD24 
#https://pubmed.ncbi.nlm.nih.gov/27157614/



phox<-ggarrange(p5,p6,p7,ncol=3)

Figure5<-ggarrange(ptop,pbottom,phox,nrow=3)
ggsave("Figure5.svg",Figure5,units="mm",dpi=300,width=180, height=215)
