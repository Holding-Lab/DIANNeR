DEResults <- read.csv('../Figure 2DEFGH + Extended Data Figure 2ABC/FragPipe Output/DE_results.csv')
contrasts <- colnames(DEResults)


#Set up conditions and LFC contrasts
conditions <- unique(
  read.csv(
    '../Figure 2DEFGH + Extended Data Figure 2ABC/FragPipe Input/E733_All_FP-A_experimental_annotation_pool_removed.tsv',
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







##Epis vs 231
'%ni%' <- Negate('%in%')

HOX_genes <- grep("HOX", DEResults$Gene.Name) %ni% grep(";", DEResults$Gene.Name) 
DEResultsHOX <- DEResults[c(setdiff(1:nrow(DEResults), HOX_genes), HOX_genes), ]


Epis231Labels<-c("GATA3", 'GRHL2','GRHL1','ZEB1','ZEB2')  
HOX<-DEResults$Gene.Name[grep("HOX", DEResults$Gene.Name)]
HOX<-HOX[grep(";", HOX, invert=TRUE) ]

colVector <- rep('#D0D0D010', length(DEResultsHOX$'Gene.Name'))
colVector[na.omit(match(DEResultsHOX$Gene.Name[grepl("HOX", DEResultsHOX$Gene.Name)], DEResultsHOX$Gene.Name))] <- '#FF0000FF'
colVector[na.omit(match(DEResultsHOX$Gene.Name[grepl("HOX", DEResultsHOX$Gene.Name)], DEResultsHOX$Gene.Name))] <- '#FF0000FF'
names(colVector)[colVector == '#D0D0D05'] <- 'Neither'



p4<-EnhancedVolcano(DEResultsHOX,
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
                    selectLab = c(Epis231Labels),
                    #selectLab = c("")
) + theme_classic()+ theme(legend.position = "none")





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
                    selectLab = c(HOX),
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



Figure4ABCD<-ggarrange(p4,p5,p6,p7,ncol=4)
Figure4ABCD

