source('Figure 5F - ChIPAtlas.R')
setwd("../Figure 3HIJ + Extended Data Figure 3C")
source('Figure 3H - ChIP-Atlas_MCF10A.R')
source('Figure 3H - ChIP-Atlas_MCF7.R')

htBRCA
htMCF10A
htMCF7

htMCF10A_small<-
    draw(htMCF10A, newpage = FALSE, 
     show_heatmap_legend = FALSE, 
     show_annotation_legend = FALSE)

htBRCA_small<-
  draw(htBRCA, newpage = FALSE, 
       show_heatmap_legend = FALSE, 
       show_annotation_legend = FALSE)

htMCF7_small<-
  draw(htMCF7, newpage = FALSE, 
       show_heatmap_legend = FALSE, 
       show_annotation_legend = FALSE)

library(grid)

#svg(filename="chipAtlas.svg", width=7, height=8)

grid.newpage()
lyt <- grid.layout(
  nrow = 3, ncol = 4,
  widths  = unit.c(unit(7*0.34, "in"), unit(7*0.19, "in"), unit(7*0.165, "in"), unit(7*0.4 , "in")),
  heights = unit.c(unit(4.9, "in"),unit(2.1, "in"), unit(0.1, "in"))
)
pushViewport(viewport(layout = lyt))

# Left: full height (rows 1:2), col 1
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(htMCF7_small,
     newpage = FALSE)
popViewport()

# Top-right: row 1, col 2
pushViewport(viewport(layout.pos.row = 1:2, layout.pos.col = 2))
draw(htMCF10A_small,
     newpage = FALSE,
     annotation_legend_side = "left")
popViewport()

# Bottom-right: row 2, col 2
pushViewport(viewport(layout.pos.row = 1:3, layout.pos.col = 3))
draw(htBRCA_small,
     newpage = FALSE)
popViewport()


Lg1<-Legend(title="Colocalization Score", labels= htMCF7@matrix_legend_param$labels,legend_gp = gpar(fill = rev(htMCF7@matrix_color_mapping@colors)))


drawAnno<-function(anno, title="",categories="") {
  print(title)
  if (title=="") {title=anno@label}
  if (any(categories=="")) {categories=anno@color_mapping@levels}
  Legend(title=title,
       labels=categories,
       legend_gp = gpar(fill=anno@color_mapping@colors)
       )
}

Lg2<-drawAnno(htMCF7@left_annotation@anno_list$Antigen)
Lg3<-drawAnno(htMCF7@left_annotation@anno_list$Treatment1, title="Treatment (Row)")
Lg4<-drawAnno(htMCF7@top_annotation@anno_list$`Pre-treatment`)
Lg5<-drawAnno(htMCF7@top_annotation@anno_list$Treatment, title="Treatment (Column)")
Lg6<-drawAnno(htMCF7@top_annotation@anno_list$Timepoint, categories=c("1-2 hours","20-30 mins"))
Lg7<-drawAnno(htBRCA@left_annotation@anno_list$Sex)
Lg8<-drawAnno(htBRCA@left_annotation@anno_list$Description)
Lg9<-drawAnno(htMCF10A@left_annotation@anno_list$CellLine, title="Cell line", categories=c("MCF10A","MCF7"))


pd<-packLegend(Lg1, Lg2, Lg3, Lg4, Lg5, Lg6, Lg7, Lg8, Lg9)



pushViewport(viewport(layout.pos.row = 1:3, layout.pos.col = 4))
draw(pd)
popViewport()

popViewport()

#dev.off()
