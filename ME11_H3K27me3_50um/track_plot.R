library(ArchR)
library(Seurat)
library(grid)
library(ggplot2)
library(patchwork)
library(dplyr)

project <- loadArchRProject("~/spatial_cuttag_project/ME11_H3K27me3_50um/Proj_with_terr")

# genome track plot
groups <- c("isolated",as.character(1:16)) # all of the territory labels

feature <- c('Hand2',"Gata6") # marker genes for heart

p <- plotBrowserTrack(
  ArchRProj = project, 
  groupBy = "territory",
  minCells = 0,
  normMethod ="nFrags",
  geneSymbol = feature, 
  upstream = 20000,
  downstream = 50000
)

grid.newpage()
grid.draw(p$Hand2)

ArchRBrowser(project)
