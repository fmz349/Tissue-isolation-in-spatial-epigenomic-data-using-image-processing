library(Seurat)
library(Signac)
library(ArchR)
library(dplyr)

project <- loadArchRProject(path = "~/ME11_H3K27me3_50um/Save-Proj_in_tissue1", 
                            force = FALSE, 
                            showLogo = TRUE)

ter <- readRDS("~/ME11_H3K27me3_50um/new_vesalius/ter")

iso_vesalius <- ter@territories

## add territory label to archR project
# add barcodes to cellcoldata in archR project
cellcoldata <- getCellColData(project)
new_row_names <- row.names(cellcoldata)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
cellcoldata$barcodes <- new_row_names

# get territory labels
terr.df <- iso_vesalius %>% select(barcodes, Territory_Trial_1)
terr.df2 <- terr.df[!duplicated(terr.df$barcodes),]

# combine cellcoldata and territory labels
combine <- merge(terr.df2,cellcoldata,by="barcodes")

new_row_names <- combine$barcodes
new_row_names <- sprintf("%s-1", new_row_names)
new_row_names <- sprintf("ME11_H3K27me3_50um#%s", new_row_names)

rownames(combine) <- new_row_names

# add cellcoldata to archR project
project_terr <- addCellColData(
  ArchRProj = project,
  data = combine$Territory_Trial_1,
  name = "territory",
  cells = rownames(combine),
  force = FALSE)

# remove cell with territory NA
new_cellcoldata <- getCellColData(project_terr)
na <- which(is.na(new_cellcoldata$territory))
cells <- rownames(new_cellcoldata)[-na]
project_terr <- subsetCells(ArchRProj = project_terr, cellNames = cells)

# test (true)

a <- rownames(new_cellcoldata)
b <- new_cellcoldata$territory
c <- combine[a,]$territory
b<-b[!is.na(b)]
c<-c[!is.na(c)]
all(b==c)


saveArchRProject(ArchRProj = project_terr,
                 outputDirectory = "~/ME11_H3K27me3_50um/new_vesalius/Proj_with_terr",
                 load = FALSE)

# Identify the marker genes for each territory
markersGS <- getMarkerFeatures(
  ArchRProj = project_terr, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "territory",
  testMethod = "wilcoxon"
)

markerList_pos <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1")
markerList_neg <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC <= -1")

saveRDS(markerList_pos,"~/ME11_H3K27me3_50um/new_vesalius/markerGeneList_pos")
saveRDS(markerList_neg,"~/ME11_H3K27me3_50um/new_vesalius/markerGeneList_neg")


write.csv(
  markerList_neg[[11]],
  "~/ME11_H3K27me3_50um/new_vesalius/CSS_whitematter.csv",
  row.names = FALSE
)

write.csv(
  markerList_neg[[6]],
  "~/ME11_H3K27me3_50um/new_vesalius/CSS_blackmatter.csv",
  row.names = FALSE
)

write.csv(
  markerList_neg[[3]],
  "~/ME11_H3K27me3_50um/new_vesalius/CSS_innerliver.csv",
  row.names = FALSE
)

write.csv(
  markerList_neg[[2]],
  "~/ME11_H3K27me3_50um/new_vesalius/CSS_outerliver.csv",
  row.names = FALSE
)

write.csv(
  markerList_neg[[10]],
  "~/ME11_H3K27me3_50um/new_vesalius/CSS_unknown.csv",
  row.names = FALSE
)

write.csv(
  markerList_neg[[8]],
  "~/ME11_H3K27me3_50um/new_vesalius/CSS_heart.csv",
  row.names = FALSE
)
write.csv(
  markerList_neg[[14]],
  "~/ME11_H3K27me3_50um/new_vesalius/CSS_14.csv",
  row.names = FALSE
)
write.csv(
  markerList_neg[[12]],
  "~/ME11_H3K27me3_50um/new_vesalius/CSS_12.csv",
  row.names = FALSE
)
write.csv(
  markerList_neg[[15]],
  "~/ME11_H3K27me3_50um/new_vesalius/CSS_15.csv",
  row.names = FALSE
)
# compare white and grey matter
markersGS2 <- getMarkerFeatures(
  ArchRProj = project_terr, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "territory",
  useGroups = "11",
  bgdGroups = "6",
  testMethod = "wilcoxon"
)
markerList_pos2 <- getMarkers(markersGS2, cutOff = "FDR <= 0.05 & Log2FC >= 1")
markerList_neg2 <- getMarkers(markersGS2, cutOff = "FDR <= 0.05 & Log2FC <= -1")

saveRDS(markerList_pos2,"~/ME11_H3K27me3_50um/new_vesalius/markerGeneList_pos_forebrain")
saveRDS(markerList_neg2,"~/ME11_H3K27me3_50um/new_vesalius/markerGeneList_neg_forebrain")

write.csv(
  markerList_neg2[[1]],
  "~/ME11_H3K27me3_50um/new_vesalius/highexpressed.whitematter.csv",
  row.names = FALSE
)
write.csv(
  markerList_pos2[[1]],
  "~/ME11_H3K27me3_50um/new_vesalius/highexpressed.greymatter.csv",
  row.names = FALSE
)

# compare inner and outer liver
markersGS3 <- getMarkerFeatures(
  ArchRProj = project_terr, 
  useMatrix = "TileMatrix", 
  groupBy = "territory",
  useGroups = "3",
  bgdGroups = "2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList_pos3 <- getMarkers(markersGS3, cutOff = "FDR <= 0.2 & Log2FC >= 1")
markerList_neg3 <- getMarkers(markersGS3, cutOff = "FDR <= 0.2 & Log2FC <= -1")

saveRDS(markerList_pos3,"~/ME11_H3K27me3_50um/new_vesalius/markerGeneList_pos_liver")
saveRDS(markerList_neg3,"~/ME11_H3K27me3_50um/new_vesalius/markerGeneList_neg_liver")

write.csv(
  markerList_neg3[[1]],
  "~/ME11_H3K27me3_50um/new_vesalius/highexpressed.innerliver.csv",
  row.names = FALSE
)
write.csv(
  markerList_pos3[[1]],
  "~/ME11_H3K27me3_50um/new_vesalius/highexpressed.outerliver.csv",
  row.names = FALSE
)

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersGS3, 
  cutOff = "FDR <= 0.2 & Log2FC >= 0.5",
  transpose = TRUE,
  plotLog2FC = TRUE
)

draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

ComplexHeatmap::draw(heatmapPeaks,
                     heatmap_legend_side = "bot",
                     annotation_legend_side = "bot")





