# spatial_data_visualization
library(ArchR)
library(Seurat)
library(grid)
library(ggplot2)
library(patchwork)
library(dplyr)

source('~/spatial_cuttag_project/downloaded_Data_visualization_code/Data_visualization/getGeneScore_ArchR.R')
source('~/spatial_cuttag_project/downloaded_Data_visualization_code/Data_visualization/SpatialPlot_new.R')

project <- loadArchRProject("~/ME11_H3K27me3_50um/new_vesalius/Proj_with_terr")

## Prepare meta data
meta.data <- as.data.frame(getCellColData(ArchRProj = project))
meta.data['cellID_archr'] <- row.names(meta.data)
new_row_names <- row.names(meta.data)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data) <- new_row_names

## Create Seurat object from imputed gene score matrix
## generate imputed gene score matrix
project <- addIterativeLSI(
  ArchRProj = project,
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

project <- addImputeWeights(
  ArchRProj = project,
  reducedDims = "IterativeLSI",
  dimsToUse = NULL,
  scaleDims = TRUE,
  corCutOff = 0.75,
  td = 1,
  ka = 2,
  sampleCells = 500,
  nRep = 1,
  k = 5,
  epsilon = 1,
  useHdf5 = TRUE,
  randomSuffix = FALSE,
  threads = getArchRThreads(),
  seed = 1,
  verbose = TRUE,
  logFile = createLogFile("addImputeWeights")
)
saveArchRProject(ArchRProj = project,
                 outputDirectory = "~/ME11_H3K27me3_50um/new_vesalius/project_lsi_terr",
                 load = FALSE)


# load in marker genes
data_dir <- '~/ME11_H3K27me3_50um/new_vesalius/highexpressed.greymatter.csv'
markerList <- read.csv(data_dir, header = TRUE, stringsAsFactors = FALSE)
markers <- markerList$name

# get imputed gene score matrix
gene_score <-
  getGeneScore_ArchR(
    ArchRProj = project,
    name = markers,
    imputeWeights = getImputeWeights(project)
  )

# creat seurat object
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"

object <- CreateSeuratObject(counts = gene_score, assay = assay, meta.data = meta.data)

image <-
  Read10X_Image(image.dir = "~/ME11_H3K27me3_50um/spatial", 
                filter.matrix = filter.matrix)
image <- image[Cells(x = object)]
DefaultAssay(object = image) <- assay
object[[slice]] <- image

spatial.obj <- object

feature <- 'Zic4'
feature <- "Gata6"

p <-
  SpatialPlot_new(
    spatial.obj,
    features = feature,
    pt.size.factor = 3,
    image.alpha = 0,
    stroke = 0,
    min.cutoff ="q3",
    max.cutoff = "q90"
  ) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15)
  )

p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22)
p

min.cutoff <- mapply(
  FUN = function(cutoff, feature) {
    return(ifelse(
      test = is.na(x = cutoff),
      yes = min(data[,sfeature]),
      no = cutoff
    ))
  },
  cutoff = min.cutoff,
  feature = features
)



################################### log(raw_gene_score +1)########################################

genescore_matrix <- getMatrixFromProject(project,useMatrix = "GeneScoreMatrix", binarize = TRUE)
matrix <- Matrix(genescore_matrix@assays@data@listData$GeneScoreMatrix)
matrix <- log2(matrix+1)
matrix@Dimnames[[1]] <- genescore_matrix@elementMetadata@listData$name

# creat seurat object
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"

object <- CreateSeuratObject(counts = matrix, assay = assay, meta.data = meta.data)

image <-
  Read10X_Image(image.dir = "~/ME11_H3K27me3_50um/spatial", 
                filter.matrix = filter.matrix)
image <- image[Cells(x = object)]
DefaultAssay(object = image) <- assay
object[[slice]] <- image

spatial.obj <- object

# make gene score plot
feature <- 'Hand2'
p <-
  SpatialPlot_new(
    spatial.obj,
    features = feature,
    pt.size.factor = 3,
    image.alpha = 0,
    stroke = 0,
    min.cutoff = "q10",
    max.cutoff = "q90"
  ) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15)
  )

p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape = 22)
p

##################################################################################################



