library(ArchR)
library(Seurat)
library(grid)
library(Signac)
library(patchwork)
library(tidyverse)

addArchRThreads(threads = 5)
addArchRGenome("mm10")

K4.spatial <- readRDS("~/Brain_H3K4me3_20um/K4.spatial.seurat.object")

setwd("~/spatial_cuttag_project/scCUTandTAG_data/mouse-brain-H3K4me3")
K4.single_cell <- readRDS("./seurat.rds")

print("both datasets have been loaded")

assay = 'bins_5000'

DefaultAssay(K4.spatial)      <- assay
DefaultAssay(K4.single_cell)  <- assay

## select peaks
min_reads = 5
# select peaks in common
peaks.common <- table(c(rownames(K4.spatial),rownames(K4.single_cell))) == 2
peaks.common <- peaks.common[peaks.common]
# select common peaks that have more than 5 reads in both datasets
features.common.table <- table(c(rownames(K4.spatial)[rowSums(K4.spatial[[assay]]@data) > min_reads],
                                 rownames(K4.single_cell)[rowSums(K4.single_cell[[assay]]@counts) > min_reads]))
peaks.use <- names(features.common.table[features.common.table == 2])

anchors <- FindIntegrationAnchors(
  object.list = list(K4.single_cell,K4.spatial),
  anchor.features = peaks.use,
  assay = rep(assay,2),
  k.filter = NA,reference = 1
)
print("integration anchors have been found")


integrated <- IntegrateData(
  anchorset = anchors,
  preserve.order = TRUE
)
print("integrated data has been generated")

integrated <- RunSVD(
  object = integrated,
  n = 50,
  reduction.name = 'integratedLSI'
)
print("svd has finished")

integrated <- RunUMAP(
  object = integrated,
  dims = 2:40,
  reduction = 'integratedLSI')
print("umap has finished")

saveRDS(integrated,"~/Brain_H3K4me3_20um/new_integrated_scCUTandTAG.seurat.object")

ter <- readRDS("~/Brain_H3K27me3_20um/new_vesalius/new_vesalius_ter.Rds")
territory <- ter@territories["Territory_Trial_1"]
rownames(territory) <- ter@territories$barcodes

cell.type <- integrated@meta.data["cell_type"]

for (i in 1:nrow(territory)) {
  cellid <- rownames(territory)[i]
  cell.type[cellid,"cell_type"] <- territory[cellid,]
}

integrated <- AddMetaData(integrated,cell.type)

integrated2 <- integrated[,!is.na(integrated@meta.data$cell_type)]

Idents(integrated2) <- 'cell_type'

p1 <- DimPlot(integrated2[,integrated2$orig.ident != 'SeuratProject'],pt.size=0.1,label=TRUE)
p2 <- DimPlot(integrated2[,integrated2$orig.ident == 'SeuratProject'],pt.size=0.1,label=TRUE)

p1+p2









