library(FragPipeAnalystR)
load("FragPipe Output/Imputed_SE.RData")
FPRime<-RData

#Impute done in preprocessing
#FPRime <- impute(FPRime, fun = "mixed")

FPRime <- test_limma(FPRime, type="all")
FPRime <- add_rejections(FPRime, alpha = 0.01, lfc = 1)

FPRime$Antibody <- ifelse(grepl("GR", FPRime$label), "GR", 
                          ifelse(grepl("IgG", FPRime$label), "IgG", 
                                 "Other"))

FPRime$Prefix <- trimws(sub("_.*", "", FPRime$label))

tissue_map <- c(
  "PrimaryCD4" = "CD4+T Cell",
  "PrimaryBreastEpis" = "Breast Epithelial",
  "PrimaryNHU" = "Urothelial",
  "MCF7" = "MCF7",
  "231" = "MDA-MB-231",
  "Jurkat" = "Jurkat",
  "KMBC2" = "KMBC2",
  "MCF10A" = "MCF10A",
  "BB7" = "ER+ PDX",
  "BB3RC31" = "ER+ PDX",
  "HBC34" = "ER+ PDX"
)


oncology_status_map <- c(
  "PrimaryCD4" = "Normal Primary",
  "PrimaryBreastEpis" = "Normal Primary",
  "PrimaryNHU" = "Normal Primary",
  
  "MCF7" = "Cancer Cell Line",
  "231" = "Cancer Cell Line",
  "Jurkat" = "Cancer Cell Line",
  "KMBC2" = "Cancer Cell Line",
  
  "MCF10A" = "Normal Cell Line",
  
  "BB7" = "Cancer PDX",
  "BB3RC31" = "Cancer PDX",
  "HBC34" = "Cancer PDX"
)

origin_map <- c(
  "PrimaryCD4" = "Lymphoid",
  "Jurkat" = "Lymphoid",
  
  "PrimaryBreastEpis" = "Epithelium",
  "MCF7" = "Epithelium",
  "231" = "Epithelium",
  "MCF10A" = "Epithelium",
  "BB7" = "Epithelium",
  "BB3RC31" = "Epithelium",
  "HBC34" = "Epithelium",
  
  "KMBC2" = "Urothelium",
  "PrimaryNHU" = "Urothelium"
)

pretty_names <- c(
  "PrimaryCD4_plus_T" = "CD4+",
  "PrimaryBreastEpis" = "Primary Breast",
  "PrimaryNHU" = "Primary Urothelium",
  "MCF7" = "MCF7",
  "231" = "MDA-MB-231",
  "Jurkat" = "Jurkat",
  "KMBC2" = "KMBC2",
  "MCF10A" = "MCF10A",
  "BB7" = "ER+ PDX (BB7)",
  "BB3RC31" = "ER+ PDX (BB3RC31)",
  "HBC34" = "ER+ PDX (HBC34)"
)

# Function to parse and reformat sample names
format_sample_name <- function(name) {
  # Extract components
  is_pool <- grepl("Pool", name)
  ab <- if (grepl("GR", name)) "GR" else if (grepl("IgG", name)) "IgG" else "Other"
  rep_match <- regmatches(name, regexpr("_\\d+$", name))
  rep <- if (is_pool) "Pool" else gsub("_", " Rep ", rep_match)
  
  # Extract prefix (this will match keys in pretty_names)
  prefix <- gsub("_.*", "", name)
  
  # Special cases: PrimaryCD4_plus_T and PrimaryBreastEpis need full prefix
  if (grepl("^PrimaryCD4_plus_T", name)) prefix <- "PrimaryCD4_plus_T"
  if (grepl("^PrimaryBreastEpis", name)) prefix <- "PrimaryBreastEpis"
  if (grepl("^PrimaryNHU", name)) prefix <- "PrimaryNHU"
  
  # Get display name
  label <- pretty_names[prefix]
  if (is.na(label)) label <- prefix
  
  paste(label, rep, ab)
}

FPRime$pretty_label <- vapply(FPRime$sample_name, format_sample_name, character(1))


FPRime$CellularIdentity <- tissue_map[FPRime$Prefix]
FPRime$Histology  <- origin_map[FPRime$Prefix]
FPRime$Pool <- ifelse(grepl("Pool", FPRime$label), "Pool", "Indiviual")
FPRime$ModelType <- oncology_status_map[FPRime$Prefix]



ht<- get_cluster_heatmap(FPRime, type="centered", 
                    indicate = c("Antibody", "ModelType", "CellularIdentity","Histology", "Pool"),
                    show_row_names = FALSE,
                    show_heatmap_legend = FALSE,
                    column_labels=FPRime$pretty_label
                    #use_raster=FALSE Made a 120MB Svg that wouldn't load.
                          )

svg("Figure 3C - cluster_heatmap.svg", width = 14, height = 14)
    draw(ht[[1]])
dev.off()


png("Figure 3C - cluster_heatmap.png", width = 14, height = 14, units="in", res=300)
  draw(ht[[1]])
dev.off()

#Attempt to redraw with  column_names_gp = grid::gpar(fontsize = 8) - but parameter
#is not accepted by draw function?
#library(ComplexHeatmap)
#draw(ht[[1]])
