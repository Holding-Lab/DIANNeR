#SRX lookup http://chip-atlas.org/view?id=SRX8520805
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

df<-read.delim("csv/NR3C1.Breast.tsv")
rowMetaData<-read.delim("csv/ChIPAtlasRowMetaDataPatient.csv", sep=",")

tfs<-c("AR", #TFs which we find good chip data for
       "PGR",
       "ESR1")

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
# 1) Filter to Tumor/PDX
# 2) Filter to just AR/PR/ER
# 3) Filter samples with no overlap
cellsOfInterest<-c("Breast_cancer")   
subclassOfInterest<-c("Breast_cancer","Breast_cancer_cells")


meta_cols <- c("Experiment","Cell_subclass","Protein")
sample_cols <- setdiff(names(df),c(meta_cols,"NR3C1.Average", "STRING"))

col_meta <- tibble(sample = sample_cols) |>
  rowwise() |>
  mutate(parsed = list(parse_col(sample)),
         SRX = parsed$id,
         CellClass = parsed$cell) |>
  ungroup() |>
  select(sample, SRX, CellClass)

col_meta_filtered<-col_meta[col_meta$CellClass %in% cellsOfInterest ,] 
sample_cols_filtered <- col_meta_filtered$sample

#Filter 1: Selects only Tumor/PDX  columns and Rows
df<-df[,c(colnames(df)[1:4],sample_cols_filtered) ]
df<-df[df$Cell_subclass %in% subclassOfInterest,]

#Filter 2: AR/PR/ER rows
df<-df[df$Protein %in% tfs,]

#Filter 3: Samples with no overlap
#df<-df[ df$Experiment %in% rowMetaData$SRX[rowMetaData$HeatmapTRUE],]
df<-df[df$SRX3229187.Breast_cancer > 0,]


#Create numeric DF
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

mat_sorted <- mat[row_order, cell_order]
col_meta_sorted <- col_meta_filtered[cell_order, , drop = FALSE]


#Tidy col names
#colnames(mat_sorted)<-col_meta$SRX[match(colnames(mat_sorted), col_meta$sample)]

mat_sorted<-t(t(mat_sorted))
colnames(mat_sorted)<-"SRX3229187"


# Simple reproducible color maps 
make_cols <- function(x,palette) {
  ux <- unique(as.character(sort(x)))
  stats::setNames(palette, ux)
}

paletteSex<-c("#F4C2C2", "#89CFF0")
paletteDesc<-c("#09546BFF", "#FC7B52FF")


#Hand curated meta data from ChIP-atlas
 ha_cols <- {
   # Align metadata rows to your heatmap matrix columns (replace `mat` if your matrix has a different name)
   
   
   meta <- data.frame(Sex="Male",Description="Tumour")
   rownames(meta) <-"SRX3229187"
   
   
      HeatmapAnnotation(
     df = meta,
     annotation_name_side = "left",
     col = list(
       "Description" = make_cols(meta[["Description"]],paletteDesc[2]), #had to use [2] because there is only one but the pallete is for rows too
       "Sex"     = make_cols(meta[["Sex"]],paletteSex[2])
     ),
     # Show annotation names on the top annotation, matching legend titles
     show_annotation_name = TRUE,
     annotation_name_gp   = grid::gpar(fontsize = 11,fontface = "bold"),
     # Legend titles consistent with annotation labels
     annotation_legend_param = list(
       "Sex" = list(title = "Sex"),
       "Description"     = list(title = "Descriptions")
     )
   )
 }

# 
# #Meta data for rows from ChIP-atlas
# 
 stopifnot(all(rownames(mat_sorted) %in% rowMetaData$SRX))
 
# Align to heatmap rows (which are SRX)
rmeta <- rowMetaData[rowMetaData$SRX %in% rownames(mat_sorted),]
rownames(rmeta) <- rmeta$SRX
rmata <- rmeta[rownames(mat_sorted), c("Antigen", "Sex","Description"), drop = FALSE] 
rmata[] <- lapply(rmata, \(x) as.character(x))  # ensure discrete annotations
 
paletteAntigen<-c("#2E008B", "#EC008C", "#00B6ED")



 ra_rows <- rowAnnotation(
   df  = rmata,
   col = list(
     Antigen    = make_cols(rmata$Antigen,paletteAntigen),
     Description =  make_cols(rmata$Description,paletteDesc),
     Sex =  make_cols(rmata$Sex,paletteSex)
   ),
   show_annotation_name = TRUE,
   annotation_name_gp   = grid::gpar(fontsize = 11, fontface = "bold"),
   annotation_legend_param = list(
     Antigen    = list(title = "Antigen"),
     Description = list(title = "Description"),
     Sex = list(title = "Sex")
   )
 )

 
#Sorting for heatmap 

OrderBySex<-match(
      rownames(rmata)[order(rmata$Sex,rmata$Antigen)],
      rownames(mat_sorted)
      )



ha_blank <- columnAnnotation(spacer = anno_empty(height = unit(14, "pt"), border = FALSE),
                             show_legend = FALSE)
 
# ----- Heatmap  -----
htBRCA <- Heatmap(
  mat_sorted,
  name = "Colocalization score",
  col = score_cols,
  top_annotation  = c(ha_cols , ha_blank),
  left_annotation = ra_rows,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = FALSE,     
  cluster_columns = FALSE, 
  row_order = OrderBySex,        
  column_order = NULL,   
  split = rmata$Sex[match(rownames(mat_sorted), rownames(rmata))],
  #column_title = "PDX and Tumour\nGR datasets",
  column_title_gp = grid::gpar(fontsize = 11),
  row_title = "Overlap with AR, ER and PGR binding sites",
  row_title_gp =  grid::gpar(fontsize = 11),
  heatmap_legend_param = list(at = legend_at, labels = legend_labels),
  row_names_gp = grid::gpar(fontsize = 8),
  rect_gp = gpar(col = "white", lwd = 2)
)

draw(htBRCA, 
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
    )

