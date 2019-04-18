library(Seurat)
setwd("~/data/Gupta/") # set working directory

rep1.data <- Read10X(data.dir = "P1Cereb_Rep1_Data/")

rep1 <- CreateSeuratObject(raw.data = rep1.data, min.cells = 3, min.genes = 200,
                           project = "scISOr_Seq_P1_CB")
mito.genes <- grep(pattern = "^mt-", x = rownames(x = rep1@data), value = TRUE)
percent.mito <- Matrix::colSums(rep1@raw.data[mito.genes, ])/Matrix::colSums(rep1@raw.data)
rep1 <- AddMetaData(object = rep1, metadata = percent.mito, col.name = "percent.mito")

hist(percent.mito,breaks=100,main="Mitochondrial gene expression % per cell(Replicate 1)");

no.of.genes = colSums(rep1@raw.data>0)
hist( no.of.genes ,breaks=100,main="No of genes per cell (Replicate 1)");
abline(v=2500,col="red",lwd=1,lty=2)



rep1= FilterCells(object = rep1, subset.names = c("nGene"),
                  low.thresholds = c(200), high.thresholds = c(2500))

rep1 <- NormalizeData(object = rep1, normalization.method = "LogNormalize",
                      scale.factor = 10000)
rep1 <- FindVariableGenes(object = rep1, mean.function = ExpMean, dispersion.function = LogVMR,
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

rep1 <- ScaleData(object = rep1, vars.to.regress = c("nUMI", "percent.mito"))
rep1 <- RunPCA(object = rep1, pc.genes = rep1@var.genes, do.print = F, pcs.print = 1:5,
               genes.print = 2)

VizPCA(object = rep1, pcs.use = 1:2)
PCAPlot(object = rep1, dim.1 = 1, dim.2 = 2)
PCHeatmap(object = rep1, pc.use = 1:12, cells.use = 500, do.balanced = TRUE,  label.columns = FALSE, use.full = FALSE)
rep1 <- JackStraw(object =rep1, num.replicate = 20)
JackStrawPlot(object = rep1, PCs = 1:20)
PCElbowPlot(object = rep1)


rep1 <- FindClusters(object = rep1, reduction.type = "pca", dims.use = 1:20,
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)

rep1 <- RunTSNE(object = rep1, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(object = rep1)

save()

#Optional Cell annotation from the paper:
lapply( c("Neurod1","Zic1","Ptf1a","Atoh1","Hepacam","Apoe","Gdf10","Tfap2b","Igfbp7","Egfl7","Hbb-bs","C1qa","Pdgfra","Gad1","Gad2","Pnoc","Lhx9","Tcf7l2","Pcp2","Pcp4","Necab2"), function(x)FeaturePlot(object = rep1, features.plot =x, do.return = F,cols.use = c("grey", "indianred1","indianred3","indianred4"), reduction.use = "tsne",pt.size = 1) )
rep1.current.cluster.ids <- as.character(unique(rep1@ident))
rep1.new.cluster.ids = c( "IGL-DCN","BG",
                          "Unclassified","EGL",
                          "EGL","NPCs",
                          "Astrocytes","Endothelial",
                          "OPCs","IGL-DCN",
                          "Astrocytes","Precursors_Atoh1+",
                          "EGL","Precursors_Ptf1a+",
                          "Blood_cells","PCL",
                          "PCL","Microglia"
)

library(plyr)
rep1@ident <- plyr::mapvalues(x = rep1@ident,
                              from = rep1.current.cluster.ids,
                              to = rep1.new.cluster.ids)

TSNEPlot(object = rep1)

