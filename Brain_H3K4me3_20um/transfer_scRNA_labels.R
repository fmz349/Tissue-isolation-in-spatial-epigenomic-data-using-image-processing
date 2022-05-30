library(ArchR)
library(Seurat)
library(grid)
library(Signac)
library(patchwork)
library(tidyverse)
library(abind)
library(scater)
# library(SeuratDisk)
library(loomR)
library(mclust)

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


addArchRThreads(threads = 5)
addArchRGenome("mm10")

######################################## cell type identification#######################################

## integrate scRNA mouse brain

# read in loom file 
# loom <- Connect(filename = "~/dev_all.loom", mode = "r")
lfile <- connect(filename = "~/scRNA_mousebrain/l5_all.loom", 
                 mode = "r+", 
                 skip.validate = TRUE)
# count matrix
matrix <- lfile[["matrix"]][,]
matrix <- t(matrix)         # gene*cell


# cell annotation and filtering
# filter repeat cells
CellID <- lfile[["col_attrs/CellID"]][]
colnames(matrix) <- CellID

tt <- table(CellID)
cell_rep <- names(which(tt > 1))

col_del_fun <- function(x){
  cols <- which(colnames(matrix) == x)
  return(cols[2:length(cols)])
}
col_del <- unlist(lapply(cell_rep, col_del_fun))
matrix2 <- matrix[,-col_del]

Age <- lfile[["col_attrs/Age"]][] # "p21"
Class <- lfile[["col_attrs/Class"]][]
ClusterName <- lfile[["col_attrs/ClusterName"]][]
Probable_location <- lfile[["col_attrs/Probable_location"]][]
Region <- lfile[["col_attrs/Region"]][]
AnalysisProject <- lfile[["col_attrs/AnalysisProject"]][] # "Adolescent" == "OK"
PassedQC <- lfile[["col_attrs/PassedQC"]][]  # "OK"
Species <- lfile[["col_attrs/Species"]][]  # "Mm"
Taxonomy_group <- lfile[["col_attrs/Taxonomy_group"]][]  
Developmental_compartment <- lfile[["col_attrs/Developmental_compartment"]][]  

cell.anno <-
  data.frame(
    CellID,
    Age,
    Class,
    ClusterName,
    Probable_location,
    Region,
    AnalysisProject,
    PassedQC,
    Species,
    Taxonomy_group,
    Developmental_compartment
  )

cell.anno <- cell.anno[-col_del,]
rownames(cell.anno) <- cell.anno$CellID

# gene annotations
Accession <- lfile[["row_attrs/Accession"]][]
Gene <- lfile[["row_attrs/Gene"]][]

gene.anno <- data.frame(Gene)
rownames(gene.anno) <- Accession

rownames(matrix2) <- Accession

scRNA_MB <- CreateSeuratObject(counts = matrix2, project = 'scRNA_MB')
scRNA_MB <- AddMetaData(object = scRNA_MB, metadata = cell.anno)

# filter
# filter age == "p21"
scRNA_p21 <- subset(scRNA_MB, Age == "p21")

delete_region <- c(
  "Enteric nervous system",
  "Pons",
  "Medulla",
  "Hypothalamus",
  "Midbrain ventral",
  "Midbrain dorsal",
  "Thalamus",
  "Spinal cord",
  "Subcommissural organ",
  "Pons,Medullae,Cerebellum",
  "Pallidum",
  "Dorsal root ganglion",
  "Dorsal root ganglion,Sympathetic ganglion",
  "Medulla,Thalamus",
  "Sympathetic ganglion"
)
region <- unique(Region)[!(unique(Region) %in% delete_region)]

scRNA_filtered <- subset(scRNA_p21, Region %in% region)

# delete repeat genes
scRNA_filtered.raw.data <- as.matrix(GetAssayData(scRNA_filtered, slot = 'counts'))
scRNA_filtered.raw.data <- as.data.frame(scRNA_filtered.raw.data)
scRNA_filtered.raw.data <- cbind(gene.anno, scRNA_filtered.raw.data)
which(is.na(scRNA_filtered.raw.data$Gene))

tt <- table(scRNA_filtered.raw.data$Gene)
name_rep <- names(which(tt > 1))
row_del_fun <- function(x){
  rows <- which(scRNA_filtered.raw.data$Gene == x)
  return(rows[2:length(rows)] )
}
row_del <- unlist(lapply(name_rep, row_del_fun))
scRNA_filtered.raw.data <- scRNA_filtered.raw.data[-row_del, ]

row.names(scRNA_filtered.raw.data) <- scRNA_filtered.raw.data$Gene
scRNA_filtered.raw.data <- scRNA_filtered.raw.data[, -1, drop=FALSE]

scRNA_filtered <-
  CreateSeuratObject(counts = scRNA_filtered.raw.data,
                     project = 'scRNA_MB_filtered',
                     meta.data = scRNA_filtered@meta.data)

saveRDS(scRNA_filtered, file = "~/scRNA_mousebrain/scRNA_filtered_seurat_object")
saveRDS(scRNA_p21, file = "~/scRNA_mousebrain/scRNA_p21_seurat_object")

# read in the ArchRProj in tissue
proj_in_tissue <- loadArchRProject("~/Brain_H3K4me3_20um/Save-Proj_in_tissue1")

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

saveArchRProject(ArchRProj = proj_in_tissue, outputDirectory = "~/Brain_H3K4me3_20um/ArchR_Proj_in_tissue_LSI", load = FALSE)

# estimate cell identity
proj_in_tissue <- addGeneIntegrationMatrix(
  ArchRProj = proj_in_tissue, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  dimsToUse = 2:30,
  seRNA = scRNA_filtered,
  addToArrow = TRUE,
  groupRNA = "ClusterName",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  force = TRUE
)

saveArchRProject(ArchRProj = proj_in_tissue, outputDirectory = "~/Brain_H3K4me3_20um/ArchR_Proj_in_tissue_predicted", load = FALSE)

meta.data.integration <- as.data.frame(getCellColData(ArchRProj = proj_in_tissue))[, c('predictedCell', 'predictedGroup', 'predictedScore')]
new_row_names <- row.names(meta.data.integration)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data.integration) <- new_row_names

# add prediction to seurat object as metadata
spatial.obj <- readRDS("~/Brain_H3K4me3_20um/seurat_object")
spatial.obj <- AddMetaData(object = spatial.obj, metadata = meta.data.integration)

saveRDS(spatial.obj, file = "~/Brain_H3K4me3_20um/seurat_object_predicted")  # MOCA
saveRDS(spatial.obj, file = "~/Brain_H3K4me3_20um/seurat_object_predicted2") # with scrna from brain

# plot the location of a predicted cell type
# filter out low quality cells
filteredGroup <- spatial.obj$predictedGroup
transfer_score <- spatial.obj$predictedScore
low_quality <- transfer_score <= 0.25 # the first quantile
filteredGroup[low_quality] <- "lowquality"

Idents(spatial.obj) <- 'filteredGroup'

# the cell group we want to show
ids.highlight <- names(table(spatial.obj$predictedGroup))[21] # Oligodendrocyte Progenitors; 21: Premature oligodendrocyte
ids.highlight

source('~/spatial_cuttag_project/downloaded_Data_visualization_code/Data_visualization/SpatialDimPlot_new.R')

# correct the direction of image
image <- spatial.obj@images$slice1@image
r <- image[,,1]
g <- image[,,2]
b <- image[,,3]
alpha <- image[,,4]

new_r <- t(r)
new_g <- t(g)
new_b <- t(b)
new_alpha <- t(alpha)

new_r <- apply(new_r, 2, rev)
new_g <- apply(new_g, 2, rev)
new_b <- apply(new_b, 2, rev)
new_alpha <- apply(new_alpha, 2, rev)
col_list <- list(new_r,new_g,new_b,new_alpha)
new_image <- abind(col_list, along = 3)

obj <- spatial.obj

obj@images$slice1@image <- new_image

a <- obj@images$slice1@coordinates$imagerow
obj@images$slice1@coordinates$imagerow <- obj@images$slice1@coordinates$imagecol
obj@images$slice1@coordinates$imagecol <- a

obj@images$slice1@coordinates$imagerow <- 1067-obj@images$slice1@coordinates$imagerow

ids.highlight <- "TEGLU3"
ids.highlight <- "TEGLU8"
ids.highlight <- "MFOL1"
ids.highlight <-"MOL2"
ids.highlight <-"ACTE2"
ids.highlight <-"PER1"
ids.highlight <-"PVM1"
ids.highlight <-"HBGLU7"

Idents(obj) <- 'predictedGroup'
ids.highlight <- "ACNT1"
ids.highlight <- "ACTE2"
ids.highlight <- "TEGLU8"
ids.highlight <- "TEGLU3"
ids.highlight <- "MGL1"
ids.highlight <- "MFOL1"
ids.highlight <- "VLMC2"
ids.highlight <- "MEINH2"
ids.highlight <- "EPEN"
ids.highlight <- "VLMC1"
ids.highlight <- "PVM1"
ids.highlight <- "PER1"

p <- SpatialDimPlot_new(
  obj,
  cells.highlight = CellsByIdentities(object = obj, idents = ids.highlight),
  facet.highlight = TRUE,
  pt.size.factor = 6,
  alpha = c(1, 0),
  stroke = 0
)
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22)
p

saveRDS(obj,"~/Brain_H3K4me3_20um/new_vesalius/seruatobject.integrate.scRNA")

# add region and probable location
repeat_cells <- names(which(table(obj@meta.data$predictedCell)>1))
cell.anno2 <- cell.anno[cell.anno$CellID %in% obj@meta.data$predictedCell,]

meta.data <- obj@meta.data

meta.data$Region <-
  unlist(lapply(meta.data$predictedCell, function(x)
    return(cell.anno2$Region[cell.anno2$CellID == x])))

meta.data$Probable_location <-
  unlist(lapply(meta.data$predictedCell, function(x)
    return(cell.anno2$Probable_location[cell.anno2$CellID == x])))

meta.data$Taxonomy_group <-
  unlist(lapply(meta.data$predictedCell, function(x)
    return(cell.anno2$Taxonomy_group[cell.anno2$CellID == x])))

obj <- AddMetaData(object = obj, metadata = meta.data)

saveRDS(obj, file = "~/Brain_H3K4me3_20um/seurat_object_predicted_final")

ter <- readRDS("~/Brain_H3K4me3_20um/new_vesalius/new_vesalius_ter.Rds")
territory <- ter@territories["Territory_Trial_1"]
rownames(territory) <- ter@territories$barcodes

obj <- obj[,rownames(territory)]

obj <- AddMetaData(obj,metadata = territory)

Idents(obj) <- 'Territory_Trial_1'
ids.highlight <- 6

p <- SpatialDimPlot_new(
  obj,
  cells.highlight = CellsByIdentities(object = obj, idents = ids.highlight),
  facet.highlight = TRUE,
  pt.size.factor = 6,
  alpha = c(1, 0),
  stroke = 0
)
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22)
p

Idents(obj2) <- 'Region'
ids.highlight <- "Striatum dorsal"
ids.highlight <- "Striatum ventral"
ids.highlight <- "Telencephalon"   #unclear
ids.highlight <- "Midbrain dorsal,Midbrain ventral" # Striatum
ids.highlight <- "Amygdala"   #unclear
ids.highlight <- "Hypothalamus,Thalamus,Midbrain dorsal,Midbrain ventral,Pons,Medulla,Spinal cord"
ids.highlight <- "Hippocampus,Cortex"
ids.highlight <- "Cortex"
ids.highlight <- "CNS"
ids.highlight <- "Cerebellum" #unclear
ids.highlight <- "Brain"

p <- SpatialDimPlot_new(
  obj2,
  cells.highlight = CellsByIdentities(object = obj2, idents = ids.highlight),
  facet.highlight = TRUE,
  pt.size.factor = 6,
  alpha = c(1, 0),
  stroke = 0
)
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22)
p

# adjusted rand index between region annotation and vesalius territory

region.anno <- meta.data[c("Region","predictedGroup","predictedScore")]

filter.region.anno <-
  region.anno[region.anno$predictedScore > summary(region.anno$predictedScore)[[2]], ]

ter <- readRDS("~/Brain_H3K4me3_20um/new_vesalius/new_vesalius_ter.Rds")
territory <- ter@territories[c("Territory_Trial_1","barcodes")]
rownames(territory) <- territory$barcodes
territory <- territory["Territory_Trial_1"]
# merge territory with raw region.anno
com <- merge(region.anno, territory, by=0, all=TRUE)

com2 <- com[!is.na(com$Territory_Trial_1),]

simi_score2 <- adjustedRandIndex(com2$Region,com2$Territory_Trial_1)
simi_score <- adjustedRandIndex(com2$predictedGroup,com2$Territory_Trial_1)

# merge territory with filtered region.anno 
filter.com <- merge(filter.region.anno, territory, by=0, all=TRUE)
filter.com <- filter.com[!is.na(filter.com$Territory_Trial_1),]
filter.com <- filter.com[!is.na(filter.com$Region),]

adjustedRandIndex(filter.com$Region,filter.com$Territory_Trial_1)                                
adjustedRandIndex(filter.com$predictedGroup,filter.com$Territory_Trial_1)

# add new annotations based on cell type and original region annotations
Cluster <- c("TEGLU8","ACTE2","PER1","TEGLU3","VLMC1","OPC","EPEN","MGL1","MFOL1","VLMC2","MOL2","VECC","VSMCA","TEGLU11","TEGLU10","MEINH2","PVM1","TEINH5","PER2")
Location <-c("Cortex","Cortex","CNS","Cortex","CNS","CNS","CNS","CNS","CNS","CNS","CNS","CNS","CNS","Cortex","Cortex","Striatum","CNS","Cortex","CNS")

rough.anno <- data.frame(Cluster,Location)

new.anno <-
  unlist(lapply(filter.com$predictedGroup, function(x)
    return(rough.anno$Location[rough.anno$Cluster==x])))

filter.com$new.anno <- new.anno

adjustedRandIndex(filter.com$new.anno,filter.com$Territory_Trial_1)                                





