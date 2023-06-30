### contruct seurat object from tile matrix

library(ArchR)
library(Seurat)
library(grid)
library(Signac)

# set the working environment
setwd("~/spatial_cuttag_project/Brain_H3K4me3_20um/")

# get tile matrix using ArchR
addArchRGenome("mm10")
addArchRThreads(threads = 8)
# addArchRThreads(threads = 32)

inputFiles <- './GSM5622964_Brain_H3K4me3_20um.fragments.tsv.gz'
sampleNames <- 'Brain_H3K4me3_20um'

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  filterTSS = 0,
  filterFrags = 0,
  minFrags = 0,
  maxFrags = 1e+07,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  offsetPlus = 0,
  offsetMinus = 0,
  TileMatParams = list(tileSize = 5000)
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = sampleNames,
  copyArrows = TRUE
)

# Select pixels in tissue
meta.data <- as.data.frame(getCellColData(ArchRProj = proj))
meta.data['cellID_archr'] <- row.names(meta.data)

data.dir <- getwd()
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"
image <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), 
                       filter.matrix = filter.matrix)

new_row_names <- row.names(meta.data)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data) <- new_row_names

meta.data.spatial <- meta.data[row.names(image@coordinates), ]
proj_in_tissue <- proj[meta.data.spatial$cellID_archr, ]

saveArchRProject(ArchRProj = proj_in_tissue, outputDirectory = "./Save-Proj_in_tissue1", load = FALSE)

# transfer the tile matrix to a Seurat object
meta.data <- as.data.frame(getCellColData(ArchRProj = proj_in_tissue))
meta.data['cellID_archr'] <- row.names(meta.data)
new_row_names <- row.names(meta.data)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data) <- new_row_names

tile_matrix <- getMatrixFromProject(proj_in_tissue,useMatrix = "TileMatrix", binarize = TRUE)
matrix <- Matrix(tile_matrix@assays@data@listData$TileMatrix) # sparse matrix

# make the colnames of sparse matrix the same as image rownames
new_col_names <- colnames(matrix)
new_col_names <- unlist(lapply(new_col_names, function(x) gsub(".*#","", x)))
new_col_names <- unlist(lapply(new_col_names, function(x) gsub("-.*","", x)))
colnames(matrix) <- new_col_names
rownames(matrix) <- paste0("Tile", 1:nrow(matrix))

object <- CreateSeuratObject(counts = matrix, assay = assay, meta.data = meta.data)

image <- image[Cells(x = object)]

DefaultAssay(object = image) <- assay
object[[slice]] <- image

# save seurat object
saveRDS(object, file = "./seurat_object")


### perform vesalius

# source the dependency functions from Vesalius and library packages
setwd("~/spatial_cuttag_project/VesaliusDev-main/VesaliusDev-main/R")
files.sources = list.files()
files.sources <- files.sources[!endsWith(files.sources, 'Rda')]
# files.sources <- files.sources[!"objects.R"]
sapply(files.sources, source)

# library(vesalius)
library(imagerExtra)
library(ggplot2)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(tvR)
library(sp)
# Standard libraries - just in case you are a command line fiend
library(utils)
library(stats)
library(graphics)
library(grDevices)

# source the rgb_lsi function
source("~/rgb_LSI.R")


# perform vesalius 
vesaliusLSI <- rgbLSI(object)

vesalius <- buildImageArray(vesaliusLSI, cores = 5)

# Vesalius_image <- imagePlot(vesalius,as.cimg =FALSE) + theme_void() + scale_y_reverse()

# print(Vesalius_image)

# source the old version vesalius

source("~/old.version.vesalius.R")

vesalius <- equalizeHistogram(vesalius,type="ECDF")

# vesalius <- regulariseImage(vesalius, lambda = 10,niter = 200, normalise=T,verbose =F)

# vesalius <- smoothArray(vesalius,method = "box", box = 1,verbose =F)

vesalius <- iterativeSegmentation.array(vesalius,
                                        colDepth = 6,
                                        smoothIter = 10,
                                        method = c("iso","median"),
                                        sigma=1.5,
                                        box = 10,
                                        useCenter = F,
                                        verbose =F)

# vesalius <- isolateTerritories.array(vesalius, captureRadius = 0.01,minBar = 10,verbose =F)

vesalius <- isolateTerritories.array(vesalius, captureRadius = 0.1, minBar = 10,verbose =F)                                     
                                     
# imgTerritory <- territoryPlot(vesalius, randomise = TRUE,cex =15 , cex.pt=1.5) + scale_y_reverse() + scale_x_reverse()+ coord_flip()+ggtitle("H3K4me3(20μm)")

imgTerritory <- territoryPlot(vesalius, randomise = TRUE,cex =15 , cex.pt=1.5) + coord_flip()+ggtitle("H3K4me3(20μm)")
print(imgTerritory)


