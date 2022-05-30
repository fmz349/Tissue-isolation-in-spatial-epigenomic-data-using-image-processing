library(ArchR)
library(Seurat)
library(grid)
library(Signac)
library(patchwork)
library(tidyverse)

K4.spatial <- readRDS("~/Brain_H3K4me3_20um/K4.spatial.seurat.object")
K4.single_cell <- readRDS("~/Brain_H3K4me3_20um/K4.single_cell.seurat.object")

# select peaks in common
peaks.common <- table(c(rownames(K4.spatial),rownames(K4.single_cell))) == 2
peaks.common <- peaks.common[peaks.common]
# select common peaks that have more than 5 reads in both datasets
assay = 'bins_5000'
min_reads = 5
features.common.table <-
  table(c(rownames(K4.spatial)[rowSums(K4.spatial[[assay]]@data) > min_reads],
          rownames(K4.single_cell)[rowSums(K4.single_cell[[assay]]@counts) > min_reads]))
peaks.use <-
  names(features.common.table[features.common.table == 2])

anchors <-
  FindTransferAnchors(reference = K4.single_cell,
                      query = K4.spatial,
                      reference.assay = 'bins_5000',
                      query.assay = 'bins_5000',
                      reduction = "lsiproject",
                      reference.reduction = "lsi",
                      features = peaks.use,
                      dims = 2:30,
                      k.filter = NA
  )

print("the anchors have been found")

predictions <- TransferData(
  anchorset = anchors,
  refdata = K4.single_cell$cell_type,
  weight.reduction = "lsiproject"
)

print("predict finished")

write.csv(predictions,"~/Brain_H3K4me3_20um/predictions")

K4.spatial <- AddMetaData(object = K4.spatial, metadata = predictions)
print("the metadata has been added")

saveRDS(K4.spatial,"~/Brain_H3K4me3_20um/K4.spatial.withpreds.seurat.object")

obj <- readRDS(file = "~/Brain_H3K4me3_20um/seurat_object_predicted_final")
preds <- predictions["predicted.id"]
rownames(preds) <- predictions$X
obj <- AddMetaData(object = obj, metadata = preds) 

Idents(obj) <- 'predicted.id'
ids.highlight <- "OEC"
"Neurons_2"  "Astrocytes" "VLMC" "OPC" "Neurons_3"  "Microglia"  "mOL""Neurons_1" 
"OEC"
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




