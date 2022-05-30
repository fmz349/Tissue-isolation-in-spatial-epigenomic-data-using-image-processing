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
setwd("~/spatial_cuttag_project/new_Vesalius/R")
files.sources = list.files()
# files.sources = files.sources[-1]
sapply(files.sources, source)

# firstly build vesalius object from seurat object
# read in seurat object
object <- readRDS("~/ME11_H3K27me3_50um/seurat_object")
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
                               tensorResolution = 0.1,
                               cores = 5,
                               filterGrid = 0)

saveRDS(ves,file="~/ME11_H3K27me3_50um/new_vesalius/new_vesalius_embedding")

# make the grayscale plots

p <- list()

for(i in 1:30){
  p[[i]] <- imagePlot(ves, dims = i, cex = 12)+ theme_void() + scale_y_reverse()
}
do.call(grid.arrange,p)

# equalize histogram
ves_eh <- equalizeHistogram(ves, embedding = "LSI", dims = 1:30,sleft = 5,sright = 5)

# segmentation
ves <-
  imageSegmentation(
    ves_eh,
    dims = 1:30,
    colDepth = 13,
    smoothIter = 5,
    smoothType = "iso",
    sigma = 2,
    embedding = "LSI"
  )

territoryPlot(ves,cex.pt = 2)+coord_flip()+scale_x_reverse()

# territory isolation
ter <- isolateTerritories(ves,trial = "last")

territoryPlot(ter,cex.pt = 2)+coord_flip()+scale_x_reverse()

saveRDS(ter,"~/ME11_H3K27me3_50um/new_vesalius/ter")

## make plot
territoryPlot <- function(vesalius,
                          trial = "last",
                          split = FALSE,
                          randomise = TRUE,
                          cex=10,
                          cex.pt=0.25){
  #--------------------------------------------------------------------------#
  # Dirty ggplot - this is just a quick and dirty plot to show what it look
  # like
  # At the moment - I thinking the user should make their own
  # Not a prority for custom plotting functions
  #--------------------------------------------------------------------------#
  
  if(!is.null(vesalius@territories) & trial == "last"){
    trial <- colnames(vesalius@territories)[ncol(vesalius@territories)]
    ter <- vesalius@territories[,c("x","y",trial)]
    colnames(ter) <- c("x","y","territory")
    ter$territory <- as.factor(ter$territory)
    legend <- sapply(strsplit(trial,"_"),"[[",1)
  } else if(!is.null(vesalius@territories) & trial != "last") {
    if(length(grep(x = colnames(vesalius@territories),pattern = trial))==0){
      stop(paste(deparse(substitute(trial)),"is not in territory data frame"))
    }
    ter <- vesalius@territories[,c("x","y",trial)]
    colnames(ter) <- c("x","y","territory")
    ter$territory <- as.factor(ter$territory)
    legend <- sapply(strsplit(trial,"_"),"[[",1)
    
  } else {
    stop("No territories have been computed!")
  }
  
  #--------------------------------------------------------------------------#
  # Changing label order because factor can suck ass sometimes
  #--------------------------------------------------------------------------#
  
  #sorted_labels <- order(levels(ter$territory))
  #ter$territory <- factor(ter$territory) %>% fct_reorder(sorted_labels)
  
  
  #--------------------------------------------------------------------------#
  # My pure hatred for the standard ggplot rainbow colours has forced me
  # to use this palette instead - Sorry Hadely
  #--------------------------------------------------------------------------#
  ter_col <- length(levels(ter$territory))
  ter_pal <- colorRampPalette(brewer.pal(8, "Accent"))
  
  if(randomise){
    set.seed(7)
    ter_col <- sample(ter_pal(ter_col),ter_col)
  } else {
    ter_col <- ter_pal(ter_col)
  }
  if(split){
    terPlot <- ggplot(ter, aes(x,y,col = territory)) +
      geom_point(size= cex.pt, alpha = 0.65)+
      facet_wrap(~territory)+
      theme_classic() +
      scale_color_manual(values = ter_col)+
      theme(legend.text = element_text(size = cex *1.2),
            axis.text = element_text(size = cex *1.2),
            axis.title = element_text(size = cex *1.2),
            plot.tag = element_text(size=cex * 2),
            plot.title = element_text(size=cex * 1.5),
            legend.title = element_text(size=cex * 1.2)) +
      guides(colour = guide_legend(override.aes = list(size=cex * 0.3)))+
      labs(colour = legend, title = paste("Vesalius",trial),
           x = "X coordinates", y = "Y coordinates")
  } else {
    terPlot <- ggplot(ter, aes(x,y,col = territory)) +
      geom_point(size= cex.pt, alpha = 0.65)+
      theme_classic() +
      scale_color_manual(values = ter_col)+
      theme(legend.text = element_text(size = cex *1.2),
            axis.text = element_text(size = cex *1.2),
            axis.title = element_text(size = cex *1.2),
            plot.tag = element_text(size=cex * 2),
            plot.title = element_text(size=cex * 1.5),
            legend.title = element_text(size=cex * 1.2)) +
      guides(colour = guide_legend(override.aes = list(size=cex * 0.3)))+
      labs(colour = legend, title = paste("Vesalius",trial),
           x = "X coordinates", y = "Y coordinates")
  }
  
  return(terPlot)
}

# segmentation plot
territoryPlot(ves,cex.pt = 2)+coord_flip()+scale_x_reverse()+theme(
  plot.title = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank())

# territories plot
territoryPlot(ter,cex.pt = 2)+coord_flip()+scale_x_reverse()+theme(
  plot.title = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank())

territoryPlot(ter,cex.pt = 0.1,split = TRUE)+coord_flip()+scale_x_reverse()+theme(
  plot.title = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank())
