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

library(EnhancedVolcano)

#Curated list 
#https://www.genenames.org/data/genegroup/#!/group/1605

#Relevant references
#BCL11B and FOXP3 https://pmc.ncbi.nlm.nih.gov/articles/PMC6685721/



BAFComplex<-read.csv("HUGO/group-1604-baf.csv",skip=1)$Approved.symbol 
pBAFComplex<-read.csv("HUGO/group-1605-pbaf.csv",skip=1)$Approved.symbol 

SharedComplex<-BAFComplex[(BAFComplex %in% pBAFComplex)]
BAFComplexFiltered<-BAFComplex[!(BAFComplex %in% pBAFComplex)]
pBAFComplexFiltered<-pBAFComplex[!(pBAFComplex %in% BAFComplex)]


labels<-c(SharedComplex,BAFComplex,pBAFComplexFiltered)

BAFLabs<-DEResults$'Gene.Name'
BAFLabs[!(DEResults$'Gene.Name' %in% labels)]<-""

###

colVector <- rep("#86868699", length(DEResults$'Gene.Name'))
colVector[na.omit(match(SharedComplex,DEResults$'Gene.Name'))] <- "#0073C2FF" 
colVector[na.omit(match(BAFComplexFiltered,DEResults$'Gene.Name'))] <- "#EFC000FF"
colVector[na.omit(match(pBAFComplexFiltered,DEResults$'Gene.Name'))] <- "#CD534CFF"
names(colVector)[colVector == "#86868699"] <- ''
names(colVector)[colVector == "#0073C2FF" ] <- 'Shared'
names(colVector)[colVector == "#EFC000FF"] <- 'BAF'
names(colVector)[colVector == "#CD534CFF"] <- 'pBAF'


BAFLabsEpi<- DEResults$'Gene.Name'[DEResults$PrimaryCD4_plus_T_GR_vs_PrimaryBreastEpis_GR_p.adj < 0.05][
              DEResults$'Gene.Name'[DEResults$PrimaryCD4_plus_T_GR_vs_PrimaryBreastEpis_GR_p.adj < 0.05]  %in% BAFLabs]

#Plot grey points first
DEResultsReorder <- DEResults[order(!colVector=="#86868699"),]
colVectorReorder <- colVector[order(!colVector=="#86868699")]

p1<-EnhancedVolcano(DEResultsReorder,
                lab = DEResultsReorder$'Gene.Name',
                x = 'PrimaryCD4_plus_T_GR_vs_PrimaryBreastEpis_GR_log2.fold.change',
                y = 'PrimaryCD4_plus_T_GR_vs_PrimaryBreastEpis_GR_p.adj',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                title = NULL,
                subtitle = NULL,
                caption = "Primary Breast Epithelium (left) vs Primary CD4+T Cells (right)",
                drawConnectors = TRUE,
                widthConnectors = 0.4,
                boxedLabels = TRUE,
                pointSize = 1.0,
                labSize = 2.5,
                max.overlaps = 50,
                colCustom = colVectorReorder,
                selectLab = BAFLabsEpi,
                #selectLab = c(""),
                colAlpha = 1
) + theme_classic()+ theme(legend.position = "none")

p1



##Tcell vs NHU

BAFLabsNHU<- DEResults$'Gene.Name'[DEResults$PrimaryCD4_plus_T_GR_vs_PrimaryNHU_GR_p.adj < 0.05][
  DEResults$'Gene.Name'[DEResults$PrimaryCD4_plus_T_GR_vs_PrimaryNHU_GR_p.adj < 0.05]  %in% BAFLabs]


p2<-EnhancedVolcano(DEResultsReorder,
                lab = DEResultsReorder$'Gene.Name',
                x = 'PrimaryCD4_plus_T_GR_vs_PrimaryNHU_GR_log2.fold.change',
                y = 'PrimaryCD4_plus_T_GR_vs_PrimaryNHU_GR_p.adj',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                title = NULL,
                subtitle = NULL,
                caption = "Primary Urothelium (left) vs Primary CD4+T Cells (right)",
                drawConnectors = TRUE,
                widthConnectors = 0.4,
                boxedLabels = TRUE,
                pointSize = 1.0,
                labSize = 2.5,
                max.overlaps = 50,
                colCustom = colVectorReorder,
                selectLab = BAFLabsNHU,
                #selectLab = c(""),
                colAlpha = 1
) + theme_classic()+ theme(legend.position = "none")


p2

row1<-ggarrange(p1,p2,ncol=2)

###FOX proteins
labelsFOX <- DEResults$Gene.Name[grep("^FOX",DEResults$Gene.Name)]
labelsFOX 
colVectorFOX <- rep("#86868699", length(DEResults$'Gene.Name'))
colVectorFOX[na.omit(match(labelsFOX,DEResults$'Gene.Name'))] <- "#7300C2FF" 
names(colVectorFOX)[colVectorFOX == "#86868699"] <- ''
names(colVectorFOX)[colVectorFOX == "#7300C2FF" ] <- 'ForkHead proteins'

DEResultsFOXreorder <- DEResults[order(!colVectorFOX=="#86868699"),]
colVectorFOXreorder <- colVectorFOX[order(!colVectorFOX=="#86868699")]

labelsFOXEpi<- DEResults$'Gene.Name'[DEResults$PrimaryCD4_plus_T_GR_vs_PrimaryBreastEpis_GR_p.adj < 0.05][
  DEResults$'Gene.Name'[DEResults$PrimaryCD4_plus_T_GR_vs_PrimaryBreastEpis_GR_p.adj < 0.05]  %in% labelsFOX]


p3<-EnhancedVolcano(DEResultsFOXreorder,
                    lab = DEResultsFOXreorder$'Gene.Name',
                    x = 'PrimaryCD4_plus_T_GR_vs_PrimaryBreastEpis_GR_log2.fold.change',
                    y = 'PrimaryCD4_plus_T_GR_vs_PrimaryBreastEpis_GR_p.adj',
                    pCutoff = 0.05,
                    FCcutoff = 1.0,
                    title = NULL,
                    subtitle = NULL,
                    caption = "Primary Breast Epithelium (left) vs Primary CD4+T Cells (right)",
                    drawConnectors = TRUE,
                    widthConnectors = 0.4,
                    boxedLabels = TRUE,
                    pointSize = 1.0,
                    labSize = 2.5,
                    max.overlaps = 50,
                    colCustom = colVectorFOXreorder,
                    selectLab = labelsFOXEpi,
                    #selectLab = c(""),
                    colAlpha = 1
) + theme_classic()+ theme(legend.position = "none")

p3


labelsFOXNHU<- DEResults$'Gene.Name'[DEResults$PrimaryCD4_plus_T_GR_vs_PrimaryNHU_GR_p.adj < 0.05][
  DEResults$'Gene.Name'[DEResults$PrimaryCD4_plus_T_GR_vs_PrimaryNHU_GR_p.adj < 0.05]  %in% labelsFOX]



p4<-EnhancedVolcano(DEResultsFOXreorder,
                    lab = DEResultsFOXreorder$'Gene.Name',
                    x = 'PrimaryCD4_plus_T_GR_vs_PrimaryNHU_GR_log2.fold.change',
                    y = 'PrimaryCD4_plus_T_GR_vs_PrimaryNHU_GR_p.adj',
                    pCutoff = 0.05,
                    FCcutoff = 1.0,
                    title = NULL,
                    subtitle = NULL,
                    caption = "Primary Urothelium (left) vs Primary CD4+T Cells (right)",
                    drawConnectors = TRUE,
                    widthConnectors = 0.4,
                    boxedLabels = TRUE,
                    pointSize = 1.0,
                    labSize = 2.5,
                    max.overlaps = 50,
                    colCustom = colVectorFOXreorder,
                    selectLab = labelsFOXNHU,
                    #selectLab = c(""),
                    colAlpha = 1
) + theme_classic()+ theme(legend.position = "none")

p4


colVectorFOX_BRAF<-colVectorFOXreorder

colVectorFOX_BRAF[na.omit(match(SharedComplex,DEResultsFOXreorder$'Gene.Name'))] <- "#0073C2FF" 
colVectorFOX_BRAF[na.omit(match(BAFComplexFiltered,DEResultsFOXreorder$'Gene.Name'))] <- "#EFC000FF"
colVectorFOX_BRAF[na.omit(match(pBAFComplexFiltered,DEResultsFOXreorder$'Gene.Name'))] <- "#CD534CFF"
names(colVectorFOX_BRAF)[colVectorFOX_BRAF == "#0073C2FF" ] <- 'Shared'
names(colVectorFOX_BRAF)[colVectorFOX_BRAF == "#EFC000FF"] <- 'BAF'
names(colVectorFOX_BRAF)[colVectorFOX_BRAF == "#CD534CFF"] <- 'pBAF'

labelsFOXJurkat<-c(labelsFOX, BAFComplex, pBAFComplex)

labelsFOXJurkat<- DEResults$'Gene.Name'[DEResults$PrimaryCD4_plus_T_GR_vs_Jurkat_GR_p.adj < 0.05][
  DEResults$'Gene.Name'[DEResults$PrimaryCD4_plus_T_GR_vs_Jurkat_GR_p.adj < 0.05]  %in% labelsFOXJurkat]

DEResultsFOX_BRAFreorder <- DEResultsFOXreorder[order(!colVectorFOX_BRAF=="#86868699"),]
colVectorFOX_BRAFreorder <- colVectorFOX_BRAF[order(!colVectorFOX_BRAF=="#86868699")]

p5<-EnhancedVolcano(DEResultsFOX_BRAFreorder,
                    lab = DEResultsFOX_BRAFreorder$'Gene.Name',
                    x = 'PrimaryCD4_plus_T_GR_vs_Jurkat_GR_log2.fold.change',
                    y = 'PrimaryCD4_plus_T_GR_vs_Jurkat_GR_p.adj',
                    pCutoff = 0.05,
                    FCcutoff = 1.0,
                    title = NULL,
                    subtitle = NULL,
                    caption = "Jurkat Cells (left) vs Primary CD4+T Cells (right)",
                    drawConnectors = TRUE,
                    widthConnectors = 0.4,
                    boxedLabels = TRUE,
                    pointSize = 1.0,
                    labSize = 2.5,
                    max.overlaps = 50,
                    colCustom = colVectorFOX_BRAFreorder,
                    selectLab = labelsFOXJurkat,
                    #selectLab = c(""),
                    colAlpha = 1
) + theme_classic()+ theme(legend.position = "none")

p5

row2<-ggarrange(p3,p4,p5,ncol=3)


ggarrange(row1,row2,nrow=2)

#Figure6<-ggarrange(ptop,pbottom,phox,nrow=3)
#ggsave("Figure6.svg",Figure5,units="mm",dpi=300,width=180, height=215)
