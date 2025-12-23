#SRX lookup http://chip-atlas.org/view?id=SRX8520805
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

df<-read.delim("csv/NR3C1.Breast.tsv")
rowMetaData<-read.delim("csv/ChIPAtlasRowMetaMCF10A.csv", sep=",")

tfs<-c("AR", #TFs which we find good chip data for
       "PGR",
       "ESR1")

#Manual annotation from ChIP-atlas website - samples were removed if unclear timing and/or treatment
#Samples SRX3070473-76 are treated with EGF - removed
# Rep 2 is missing SRX3070468 and 70 - checking QC both have ~ 200 and 500 peaks 
# respectively, cf 60 mins (SRX..71) has 14158 peaks which is typical for SR chip
sample_id<-c("SRX3070467", "SRX3070469", "SRX3070471") #Rep 2 isn't complete, and as doesn't seem to pass QC later & samples are +EGF
# Rep 2 was therefore removed
treatment <-c("Untreated", "Dexamethasone", "Dexamethasone")
timepoint <-c("1 hour",    "20 min",        "1 hour")

col_experimental_meta<-data.frame(sample_id=sample_id,
                                  treatment=treatment,
                                  timepoint=timepoint)




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
# 1) Filter to Cols = MCF10A, Rows = MCF10A + MCF7 (lack of 10A ChIPs for AR/PR/ER)
# 2) Filter to just AR/PR/ER
# 3) Force inclusion of AR samples, MCF10A rows,
# 4) Filter on ChIP-Atlas Ranks/NR3C1.Average as too many ER chip-seqs in MCF7 to show

# Filter 1 
cellsOfInterest<-c("MCF_10A")   
subclassOfInterest<-c("MCF_10A","MCF-7")

##Filter rows
df<-df[df$Cell_subclass %in% subclassOfInterest ,] 
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

#Filter 2
df<-df[df$Protein %in% tfs,]

#Filter 3 & 4
df<-df[ df$Protein %in% c("AR") | df$Cell_subclass %in% c("MCF_10A") |
          (df$Protein=="ESR1" & df$NR3C1.Average > 0.8) |
          df$NR3C1.Average > 0.8,] 


#Create matrix of numerical data
df_num <- df
mat <- as.matrix(df_num[, sample_cols_filtered, drop = FALSE])
rownames(mat) <- make.unique(paste(df$Experiment))

# Column annotation: Cell subclass
cell_classes <- unique(col_meta_filtered$CellClass)
# Generate cell line colours - not so useful now we focus on 1 cell line. 
if (length(cell_classes) <= 12) {
  col_cell <- setNames(brewer.pal(max(3, length(cell_classes)), "Set3")[seq_along(cell_classes)], cell_classes)
} else {
  col_cell <- setNames(hcl.colors(length(cell_classes), palette = "Dynamic"), cell_classes)
}


# Heatmap row annotation - Protein/TF bait for ChIP
proteins <- unique(df$Protein)
if (length(proteins) <= 12) {
  col_prot <- setNames(brewer.pal(max(3, length(proteins)), "Paired")[seq_along(proteins)], proteins)
} else {
  col_prot <- setNames(hcl.colors(length(proteins), palette = "Zissou 1"), proteins)
}



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
  `0`  = "#7d817e",  # N.D. (grey)
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

mat_sorted <- mat[, cell_order] #row_order
col_meta_sorted <- col_meta_filtered[cell_order, , drop = FALSE]


#Tidy col names
colnames(mat_sorted)<-col_meta$SRX[match(colnames(mat_sorted), col_meta$sample)]

 rmeta <- rowMetaData
 rownames(rmeta) <- rmeta$SRX
 rmata <- rmeta[rownames(mat_sorted), c("Antigen", "CellLine"), drop = FALSE] 
 rmata[] <- lapply(rmata, \(x) as.character(x))  # ensure discrete annotations

 paletteAntigen<-c("#2E008B", "#EC008C", "#00B6ED")
 paletteCellLine<-c("#0000AA", "#EEEE00" )


# Simple reproducible color maps 
make_cols <- function(x,palette) {
  ux <- unique(as.character(sort(x)))
  stats::setNames(palette, ux)
}

#rmata<-data.frame(Antigen=df$Protein)

ra_rows <- rowAnnotation(
  df  = rmata,
  col = list(
    CellLine    = make_cols(rmata$CellLine,paletteCellLine),
    Antigen    = make_cols(c("ER","PGR","AR"),paletteAntigen)
  ),
  show_annotation_name = TRUE,
  annotation_name_gp   = grid::gpar(fontsize = 11, fontface = "bold"),
  annotation_legend_param = list(
    CellLine    = list(title = "Cell line"),
    Antigen    = list(title = "Antigen")
  )
)


paletteTreatmentH<-c( "#FF8811", "#FFBB66")
palettePreTreatment<-c("#FFFFDF","#CFFF04")
paletteTimePoint<-c("#C0C0FF","#6060FF")
 
 ha_cols <- {
   # Align metadata rows to your heatmap matrix columns (replace `mat` if your matrix has a different name)
   stopifnot(all(colnames(mat_sorted) %in% col_experimental_meta$sample_id))
   meta <- col_experimental_meta
   rownames(meta) <- meta$sample_id
   meta <- meta[colnames(mat_sorted), c( "treatment", "timepoint"), drop = FALSE]#, "antibody"
   colnames(meta) <- c( "Treatment", "Timepoint") # "Antibody"
   
   
   HeatmapAnnotation(
     df = meta,
     annotation_name_side = "left",
     col = list(
       "Treatment"     = make_cols(meta[["Treatment"]],paletteTreatmentH),
       "Timepoint"     = make_cols(meta[["Timepoint"]],paletteTimePoint)
       # "Antibody"      = make_cols(meta[["Antibody"]])
     ),
     # Show annotation names on the top annotation, matching legend titles
     show_annotation_name = TRUE,
     annotation_name_gp   = grid::gpar(fontsize = 11, fontface = "bold"),
     # Legend titles consistent with annotation labels
     annotation_legend_param = list(
       "Treatment"     = list(title = "Treatment"),
       "Timepoint"     = list(title = "Timepoint")
     )
   )
 }
 



 ha_blank <- columnAnnotation(spacer = anno_empty(height = unit(14, "pt"), border = FALSE),
                              show_legend = FALSE)
# ----- Heatmap  -----
htMCF10A <- Heatmap(
  mat_sorted,
  name = "ChIP score",
  col = score_cols,
  top_annotation  = c(ha_cols , ha_blank),
  left_annotation = ra_rows  ,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = FALSE,     
  cluster_columns = FALSE, 
  row_order = NULL,        
  column_order = NULL,     
  split = rmata$Antigen[match(rownames(mat_sorted), rownames(rmata))],
  #column_title = "MCF10A\nGR datasets",
  row_title = "Overlap with AR, ER and PGR binding sites",
  row_title_gp =  grid::gpar(fontsize = 11),
  column_title_gp = grid::gpar(fontsize = 11),
  heatmap_legend_param = list(at = legend_at, labels = legend_labels),
  row_names_gp = grid::gpar(fontsize = 8),
  rect_gp = gpar(col = "white", lwd = 2)
)

draw(htMCF10A, 
     heatmap_legend_side = "right",
     annotation_legend_side = "right",

)

#Note the MCF10A ChIP is a MYC fusion, NOT ER. And removed in the figure.
