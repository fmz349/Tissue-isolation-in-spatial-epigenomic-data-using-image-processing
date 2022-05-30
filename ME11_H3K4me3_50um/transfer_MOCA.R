library(ArchR)
library(Seurat)
library(grid)
library(Signac)
library(patchwork)
library(tidyverse)

addArchRThreads(threads = 5)
addArchRGenome("mm10")

MOCA_E11 <- readRDS("~/spatial_cuttag_project/ref_data/MOCA_E11_seurat_object")

# read in the ArchRProj in tissue
proj_in_tissue <- loadArchRProject("~/ME11_H3K4me3_50um/Save-Proj_in_tissue1")

proj_in_tissue <- addIterativeLSI(
  ArchRProj = proj_in_tissue,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list(
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force = TRUE
)

print("LSI has finished")

proj_in_tissue <- addImputeWeights(proj_in_tissue)

# estimate cell identity
proj_in_tissue <- addGeneIntegrationMatrix(
  ArchRProj = proj_in_tissue, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = MOCA_E11,
  addToArrow = TRUE,
  groupRNA = "Main_cell_type",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  force = TRUE
)

saveArchRProject(ArchRProj = proj_in_tissue,
                 outputDirectory = "~/ME11_H3K4me3_50um/ArchR_Proj_in_tissue_predicted",
                 load = FALSE)


print("prediction finished")

meta.data.integration <- as.data.frame(getCellColData(ArchRProj = proj_in_tissue))[, c('predictedCell', 'predictedGroup', 'predictedScore')]
new_row_names <- row.names(meta.data.integration)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data.integration) <- new_row_names

# add prediction to seurat object as metadata
spatial.obj <- readRDS("~/ME11_H3K4me3_50um/seurat_object")
spatial.obj <- AddMetaData(object = spatial.obj, metadata = meta.data.integration)

saveRDS(spatial.obj, file = "~/ME11_H3K4me3_50um/seurat_object_predicted")

print("integrated seraut object has been saved")

# select predict group and territory
iso_vesalius <- ter@territories
terr.df <- iso_vesalius %>% select(barcodes, Territory_Trial_1)
terr.df2 <- terr.df[!duplicated(terr.df$barcodes),]

spatial.obj@meta.data$barcodes <- rownames(spatial.obj@meta.data)

combine <- merge(terr.df2,spatial.obj@meta.data,by="barcodes")
rownames(combine) <- combine$barcodes

object2 <- subset(spatial.obj, cells = combine$barcodes)
object2 <- AddMetaData(object = object2, metadata = combine)

saveRDS(object2,"~/ME11_H3K4me3_50um/new_vesalius/seurat.predictgroup.ter")

# plot predicted group
source('~/spatial_cuttag_project/downloaded_Data_visualization_code/Data_visualization/SpatialDimPlot_new.R')

Idents(object2) <- 'predictedGroup'

# the cell group we want to show
ids.highlight <- "Limb mesenchyme"

p <- SpatialDimPlot_new(
  object2,
  cells.highlight = CellsByIdentities(object = object2, idents = ids.highlight),
  facet.highlight = TRUE,
  pt.size.factor = 3,
  alpha = c(1, 0),
  stroke = 0
)
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22)
p

# plot different ter
Idents(object2) <- 'Territory_Trial_1'
ids.highlight <- "20"

p <- SpatialDimPlot_new(
  object2,
  cells.highlight = CellsByIdentities(object = object2, idents = ids.highlight),
  facet.highlight = TRUE,
  pt.size.factor = 3,
  alpha = c(1, 0),
  stroke = 0
)
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22)
p

















