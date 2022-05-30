library(Matrix)
library(parallel)
library(imager)
library(imagerExtra)
library(Seurat)
library(tidyverse)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(tvR)
library(sp)
library(grid)
library(gridExtra)
library(ggplot2)

library(SuperpixelImageSegmentation)

library(Signac)
library(tidyverse)

# Standard libraries - just in case you are a command line fiend
library(utils)
library(stats)
library(graphics)
library(grDevices)

# load new version vesalius
setwd("~/spatial_cuttag_project/Vesalius/R")
files.sources = list.files()
files.sources = files.sources[-1]
sapply(files.sources, source)

# firstly build vesalius object from seurat object
# read in seurat object
object <- readRDS("~/Brain_H3K4me3_20um/seurat_object")
beads <- GetTissueCoordinates(object@images$slice1)
beads$barcodes <- rownames(beads)
colnames(beads) <- c("xcoord","ycoord","barcodes")

counts <- GetAssayData(object, slot = "counts")

vesalius <- buildVesaliusObject(beads,counts)

# building vesalius embeddings
ves <- buildVesaliusEmbeddings(vesalius,
                               norm = "TFIDF",
                               method ="LSI",
                               pcs = 30,
                               tensorResolution = 1,
                               cores = 5)

saveRDS(ves,file="~/Brain_H3K4me3_20um/new_vesalius/new_vesalius_embedding.Rds")
saveRDS(ves,file="~/Brain_H3K4me3_20um/new_vesalius/new_vesalius_embedding_allspots.Rds")

# make the grayscale plots for principal components

p <- list()

for(i in 1:30){
  p[[i]] <- imagePlot(ves, dims = i, cex = 12)+theme_void()
}
do.call(grid.arrange,p)


# equalize histogram
ves_eh <- equalizeHistogram(
  ves,
  embedding = "LSI",
  type = "ECDF",
  dims = 1:30,
  sleft = 5,
  sright = 5
)

p3 <- list()

for(i in 1:30){
  p3[[i]] <- imagePlot(ves_eh, dims = i, cex = 12)+theme_void()
}
do.call(grid.arrange,p3)

# image segmentation

seg <-
  imageSegmentation(
    ves,
    dims = 1:30,
    colDepth = 7,
    smoothIter = 5,
    smoothType = c("iso","median"),
    sigma = 2,
    box = 10,
    embedding = "LSI"
  )

saveRDS(seg,file="~/Brain_H3K4me3_20um/new_vesalius/new_vesalius_seg.Rds")

territoryPlot(seg,cex.pt = 2)

# territory isolation
# when captureRadius=1 there will be no isolation within a cluster
ter <- isolateTerritories(seg,
                     trial = "last",
                     captureRadius = 1)

territoryPlot(ter,cex.pt = 2)+theme(
  plot.title = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank())

saveRDS(ter,file="~/Brain_H3K4me3_20um/new_vesalius/new_vesalius_ter.Rds")

# add territory to seurat object
territory <- ter@territories[c("barcodes","Territory_Trial_1")]
rownames(territory) <- territory$barcodes
territory <- territory["Territory_Trial_1"]
object <- readRDS("~/Brain_H3K4me3_20um/seurat_object")
object <- AddMetaData(object = object, metadata = territory)

saveRDS(object,file="~/Brain_H3K4me3_20um/new_vesalius/seurat_object_ter")











