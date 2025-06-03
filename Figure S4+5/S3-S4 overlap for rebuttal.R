df<-read.delim('../Figure S2+3/Scafold/Samples Report With Clusters for D676.tsv')
dfNoClusters<-df[grep("\\.", df$X., invert=TRUE ),]
TFs_676<-dfNoClusters[dfNoClusters$D676_AT..AT. > 0,]$Alternate.ID


df<-read.delim('Scafold/Samples Report With Clusters for D714.tsv')
dfSamples<-df[,c(1:9,15,18)]
#colnames(dfSamples)[c(10,11)]<-c("Dex","DMSO") # 2 hours
colnames(dfSamples)[c(9,10)]<-c("Dex","DMSO") #Overnight
dfNoClusters<-dfSamples[grep("\\.", dfSamples$X., invert=TRUE ),]


#Dex Results
TFs_714<-dfNoClusters[dfNoClusters$Dex > 0,]$Alternate.ID


length(TFs_714)
length(TFs_676)

TFs_714 <- TFs_714[!is.na(TFs_714)]
TFs_676 <- TFs_676[!is.na(TFs_676)]


length(TFs_714[TFs_714 %in% TFs_676])

library(ggvenn)
venn_list <- list(`D714` = TFs_714, `D676` = TFs_676)

# Plot
ggvenn(venn_list,
       fill_color = c("#E69F00", "#56B4E9"),
       stroke_size = 0.5,
       set_name_size = 4)


# D714 specific 837
# Intersection = 781
# D676 specitic = 202
781/(781+202)
#[1] 0.7945066
#Low cell method (D714) detects 80% of the proteins found in the original method (D676)