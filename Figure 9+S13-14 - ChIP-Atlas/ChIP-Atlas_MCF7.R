#SRX lookup http://chip-atlas.org/view?id=SRX8520805
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

df<-read.delim("csv/NR3C1.Breast.tsv")
rowMetaData<-read.delim("csv/ChIPAtlasRowMetaData.csv", sep=",")

tfs<-c("AR", #TFs which we find good chip data for
       "PGR",
       "ESR1")

#Manual annotation from ChIP-atlas website - samples were removed if unclear timing and/or treatment
sample_id<-c("SRX1161165", "SRX1161166"   , "SRX1161167"   , "SRX19024866"  , "SRX19024868"  , "SRX19024870"  , "SRX19024878"       , "SRX19024879"       , "SRX19024880"       )  
treatment <-c("Untreated", "Dexamethasone", "Dexamethasone", "Dexamethasone", "Dexamethasone", "Dexamethasone", "Dexamethasone", "Dexamethasone" , "Dexamethasone" )                          
timepoint<-c("30 mins"   , "30 mins"      , "30 mins"      , "2 hours"      , "2 hours"      , "2 hours"      , "2 hours", "2 hours", "2 hours")
treatment2<-c("None", "None"         , "None"         , "None"         , "None"         , "None"         , "Fulvestrant"       , "Fulvestrant" , "Fulvestrant" )                          
antibody <-c("sc-1003"   , "sc-1003"      , "sc-1003"      , "12041 CST"    , "12041 CST"    , "12041 CST"    , "12041 CST"         , "12041 CST"         , "12041 CST"         )

col_experimental_meta<-data.frame(sample_id=sample_id,
                                  treatment=treatment,
                                  treatment2=treatment2,
                                  timepoint=timepoint,
                                  antibody=antibody)

meta_cols <- c("Experiment","Cell_subclass","Protein")
sample_cols <- setdiff(names(df),c(meta_cols,"NR3C1.Average", "STRING"))

parse_col <- function(s) {
  # Split on the first dot only
  parts <- strsplit(s, "\\.", fixed = FALSE)[[1]]
  if (length(parts) >= 2) {
    list(id = parts[1], cell = paste(parts[-1], collapse = "."))  # keep any extra dots in the cell subclass
  } else {
    # No dot present: try to treat as ID + Unknown cell class
    list(id = s, cell = "Unknown")
  }
}

### Filter 
# 1) Filter to MCF7
# 2) Filter to just AR/PR/ER
# 3) Filter on Filter on ChIP-Atlas Ranks/NR3C1.Average as too many ER chip-seqs in MCF7 to show
# 4) Filter conditions other then steroid hormones (+ some samples with unclear metadata)

#Filter 1: MCF7 only
cellsOfInterest<-c("MCF.7")   # MCF7 matches our proteomics
subclassOfInterest<-c("MCF-7")

##Filter Rows
df<-df[df$Cell_subclass %in% subclassOfInterest,]

##Filter Columns
col_meta <- tibble(sample = sample_cols) |>
  rowwise() |>
  mutate(parsed = list(parse_col(sample)),
         SRX = parsed$id,
         CellClass = parsed$cell) |>
  ungroup() |>
  select(sample, SRX, CellClass)

col_meta_filtered<-col_meta[col_meta$CellClass %in% cellsOfInterest & col_meta$SRX %in% sample_id,] 
sample_cols_filtered <- col_meta_filtered$sample
df<-df[,c(colnames(df)[1:4],sample_cols_filtered) ]

#Filter 2: AR/PR/ER rows only
df<-df[df$Protein %in% tfs,]


#Filter 3: Cut off as so many ER ChIP-seq samples for MCF7 and the need auditing for reagents
df<-df[ df$NR3C1.Average > 0.425,] 

#Filter 4:
steroidConditions<-c(#"R5020" , removed to tidy figure
                     "Progesterone",
                     "DHT",
                     "Tamoxifen resitance",
                     "Estradiol",
                     "Full media",
                     "Dexamethasone",
                     "Fulvestrant",
                     "Tamoxifen",
                     "No estrogen",   
                     "Long-term estrogen deprived",
                     "No progesterone or estrogen"
                     )




rowMetaData$Treatment1[rowMetaData$Treatment1=="No estrogen"]<-"No progesterone or estrogen" #to simplify legend

df<-df[df$Experiment %in% rowMetaData$SRX[rowMetaData$Treatment1 %in% steroidConditions & rowMetaData$Heatmap==TRUE],]

rowMetaData<-rowMetaData[rowMetaData$Treatment1 %in% steroidConditions,]


df_num <- df

#Create matrix of numerical data
mat <- as.matrix(df_num[, sample_cols_filtered, drop = FALSE])

# Give rows readable names; originally these contained factor and cell line. For clarity we did
# 1 heatmap for MCF7 rather than combining cell lines as overlap between cell lines didn't answer
# the reviewers question. Factor/Protein/TF was moved to the metadata annotation. 
#rownames(mat) <- make.unique(paste(df$Experiment, df$Cell_subclass, df$Protein, sep = " | "))
rownames(mat) <- make.unique(paste(df$Experiment))

# Column annotation: Cell subclass
cell_classes <- unique(col_meta_filtered$CellClass)
# Generate cell line colours - not so useful now we focus on 1 cell line. 
if (length(cell_classes) <= 12) {
  col_cell <- setNames(brewer.pal(max(3, length(cell_classes)), "Set3")[seq_along(cell_classes)], cell_classes)
} else {
  col_cell <- setNames(hcl.colors(length(cell_classes), palette = "Dynamic"), cell_classes)
}

#Heatmap column annotation - will probably remove from heatmap, but useful if we do go back to 
#looking at multiple cell lines in this figure panel.
ha_cols <- HeatmapAnnotation(
  CellClass = col_meta_filtered$CellClass,
  col = list(CellClass = col_cell),
  annotation_legend_param = list(CellClass = list(title = "Cell subclass"))
)


# Heatmap row annotation - Protein/TF bait for ChIP
proteins <- unique(df$Protein)
if (length(proteins) <= 12) {
  col_prot <- setNames(brewer.pal(max(3, length(proteins)), "Paired")[seq_along(proteins)], proteins)
} else {
  col_prot <- setNames(hcl.colors(length(proteins), palette = "Zissou 1"), proteins)
}

ra_rows <- rowAnnotation(
  Protein = df$Protein,
  col = list(Protein = col_prot),
  annotation_legend_param = list(Protein = list(title = "Protein"))
)



#ChIP atlas uses a score of 10 for the sample being the 'same' as the one being compared to.
df_num[sample_cols_filtered] <- lapply(df_num[sample_cols_filtered], function(x) {
  x <- as.character(x)
#  x[x %in% c("N.D.","ND","n.d.","na","NA")] <- "0"   # N.D. -> 0
  x[x %in% c("Same","same")] <- "10"                 # Same -> 10
  out <- suppressWarnings(as.numeric(x))
  out[is.na(out)] <- 0
  out
})

#mat <- as.matrix(df_num[, sample_cols_filtered, drop = FALSE])

# Colour map as used by chip atlas
# Pairs:
# 9=H-H, 6=H-M/M-H, 4=M-M, 3=H-L/L-H, 2=M-L/L-M, 1=L-L
# Optional: 0=N.D. (grey), 10=Same (black)
score_cols <- c(
  `0`  = "#CeCeCe",  # N.D. (grey)
  `1`  = "#2a71f5",  # L-L (blue)
  `2`  = "#64defa",  # M-L / L-M (light blue)
  `3`  = "#73f9b0",  # H-L / L-H (mint)
  `4`  = "#74f75c",  # M-M (green)
  `6`  = "#bcfa4f",  # H-M / M-H (yellow-green)
  `9`  = "#e83222",  # H-H (red)
  `10` = "#000000"   # Same (black)
)

# Pretty legend labels 
present_vals <- sort(unique(as.numeric(mat)))
legend_at     <- intersect(c(9,6,4,3,2,1,0,10), present_vals)
legend_labels <- c(
  `9`="High–High",
  `6`="High–Medium / Medium–High",
  `4`="Medium–Medium",
  `3`="High–Low / Low–High",
  `2`="Medium–Low / Low–Medium",
  `1`="Low–Low",
  `0`="Not Detected",
  `10`="Same"
)
legend_labels <- legend_labels[as.character(legend_at)]

df$Protein <- factor(df$Protein, levels = unique(df$Protein))  # keep appearance order
row_order <- order(df$Protein)

cell_order <- order(factor(col_meta_filtered$CellClass))

mat_sorted <- mat[row_order, cell_order]
col_meta_sorted <- col_meta_filtered[cell_order, , drop = FALSE]


#Tidy col names
colnames(mat_sorted)<-col_meta$SRX[match(colnames(mat_sorted), col_meta$sample)]

# Simple reproducible color maps 
make_cols <- function(x,palette) {
  ux <- unique(as.character(sort(x)))
  stats::setNames(palette, ux)
}

paletteTreatmentH<-c(  "#FF8811", "#FFBB66")
palettePreTreatment<-c("#CFFF04","#FFFFDF")
paletteTimePoint<-c("#6060FF","#C0C0FF")


#Hand curated meta data from ChIP-atlas
ha_cols <- {
  # Align metadata rows to your heatmap matrix columns (replace `mat` if your matrix has a different name)
  stopifnot(all(colnames(mat_sorted) %in% col_experimental_meta$sample_id))
  meta <- col_experimental_meta
  rownames(meta) <- meta$sample_id
  meta <- meta[colnames(mat_sorted), c("treatment2", "treatment", "timepoint"), drop = FALSE]#, "antibody"
  colnames(meta) <- c("Pre-treatment", "Treatment", "Timepoint") # "Antibody"
  
  
  HeatmapAnnotation(
    df = meta,
    annotation_name_side = "left",
    col = list(
      "Pre-treatment" = make_cols(meta[["Pre-treatment"]],palettePreTreatment),
      "Treatment"     = make_cols(meta[["Treatment"]],paletteTreatmentH),
      "Timepoint"     = make_cols(meta[["Timepoint"]],paletteTimePoint)
     # "Antibody"      = make_cols(meta[["Antibody"]])
    ),
    # Show annotation names on the top annotation, matching legend titles
    show_annotation_name = TRUE,
    annotation_name_gp   = grid::gpar(fontsize = 11, fontface = "bold"),
    # Legend titles consistent with annotation labels
    annotation_legend_param = list(
      "Pre-treatment" = list(title = "Pre-treatment"),
      "Treatment"     = list(title = "Treatment"),
      "Timepoint"     = list(title = "Timepoint")
      #"Antibody"      = list(title = "Antibody")
    )
  )
}


#Meta data for rows from ChIP-atlas

stopifnot(all(rownames(mat_sorted) %in% rowMetaData$SRX))

# Align to heatmap rows (which are SRX)
rmeta <- rowMetaData[rowMetaData$SRX %in% rownames(mat_sorted),]
rownames(rmeta) <- rmeta$SRX
rmata <- rmeta[rownames(mat_sorted), c("Antigen", "Treatment1"), drop = FALSE] #"Time", "Antibody"
rmata[] <- lapply(rmata, \(x) as.character(x))  # ensure discrete annotations



paletteAntigen<-c("#2E008B", "#EC008C", "#00B6ED")
paletteTreatment<-c("#FF8811FF", 
                    "#5544DDFF",
                    "#FFAFAFFF",
                    #"#F561C0FF",
                    "#CFFF04FF", 
                    "#DFDFDFFF",
                    "#B6DEE2FF", 
                   # "#AFFFFFFF",
                    "#168E92FF"
#                    "#E30022FF",
                #    "#36CE92FF",
               #     "#0FFF0FFF",
              #      "#AFFFAFFF"
                    )




ra_rows <- rowAnnotation(
  df  = rmata,
  col = list(
    Antigen    = make_cols(c("ER","PGR","AR"),paletteAntigen),
    Treatment1 = make_cols(rmata$Treatment1,paletteTreatment)
    #Time       = make_cols(rmata$Time),
    #Antibody   = make_cols(rmata$Antibody)
  ),
  show_annotation_name = TRUE,
  annotation_name_gp   = grid::gpar(fontsize = 11,fontface = "bold"),
  annotation_legend_param = list(
    Antigen    = list(title = "Antigen"),
    Treatment1 = list(title = "Treatment")
  #  Time       = list(title = "Time"),
  #  Antibody   = list(title = "Antibody")
  )
)

# ----- Heatmap  -----
htMCF7 <- Heatmap(
  mat_sorted,
  name = "Colocalization score",
  col = score_cols,
  top_annotation  = ha_cols,
  left_annotation = ra_rows,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = FALSE,     
  cluster_columns = FALSE, 
  row_order = NULL,        
  column_order = NULL,     
  split = rmata$Antigen[match(rownames(mat_sorted), rownames(rmata))],
  #column_title = "MCF7 GR datasets",
  column_title_gp = grid::gpar(fontsize = 11),
  row_title = "Overlap with AR, ER and PGR binding sites",
  row_title_gp =  grid::gpar(fontsize = 11),
  heatmap_legend_param = list(at = legend_at, labels = legend_labels),
  row_names_gp = grid::gpar(fontsize = 8),
  rect_gp = gpar(col = "white", lwd = 2)
)

draw(htMCF7, 
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
    )

