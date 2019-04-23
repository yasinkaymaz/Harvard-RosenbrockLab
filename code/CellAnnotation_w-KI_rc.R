.libPaths("/n/home13/yasinkaymaz/biotools/Rlibs/")
source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/functions.R")

Hook2018gene <- get(load("~/data/Hook2018gene.seurat.Robj"))
Hook2018 <- get(load("~/data/Hook2018.seurat.Robj"))
dim(Hook2018gene@meta.data)


#Annotate cell types with KI-superset
zeisel.rank3.rfcv <- get(load("~/data/zeisel.rank3.rfcv.RF_model_notImproved.Robj"))
zeisel.rank4.rfcv <- get(load("~/data/zeisel.rank4.rfcv.RF_model_notImproved.Robj"))

Hook2018gene <- CellTyper2(SeuratObject = Hook2018gene, model = zeisel.rank4.rfcv, priorLabels = Hook2018gene@meta.data$region, outputFilename = paste("Hook2018gene",".rfcv.predictions.w-singleZeisel-rank4",sep = ""))
head(Hook2018gene@meta.data)

CrossCheck(PriorPostTable = Hook2018gene@meta.data[,c("subset","Intermediate")],outputprefix = "Hook.predictions.")


Hook2018gene <- SetAllIdent(Hook2018gene, id="Intermediate")
VlnPlot(object = Hook2018gene,x.lab.rot = 45, features = c("Th", "Lhx9","Slc6a3","Ldb2"),nCol = 2)
VlnPlot(object = Hook2018gene,x.lab.rot = 45, features = c("Crhr1", "Nsf", "Mapt"),nCol = 1)
VlnPlot(object = Hook2018gene,x.lab.rot = 45, features = c("Snap25"),nCol = 1)
VlnPlot(object = Hook2018,x.lab.rot = 45, features = c("Snap25"),nCol = 1)


#Gupta Annotate cells with KI superset
Gupta2017 <- CellTyper2(SeuratObject = Gupta2017, model = zeisel.rank4.rfcv, priorLabels = Gupta2017@meta.data$res.0.6, outputFilename = paste("Gupta2017",".rfcv.predictions.w-singleZeisel-rank4",sep = ""))
CrossCheck(PriorPostTable = Gupta2017@meta.data[,c("res.0.6","Intermediate")],outputprefix = "Gupta2017.predictions.")


#Double-check cell predictions with Tasic data - adult mice
load("~/data/Tasic2018.seurat.Robj")
Tasic2018 <- CellTyper2(SeuratObject = Tasic2018, model = zeisel.rank4.rfcv, priorLabels = Tasic2018@meta.data$cluster, outputFilename = paste("Tasic2018",".rfcv.predictions.w-singleZeisel-rank4",sep = ""))
CrossCheck(PriorPostTable = Tasic2018@meta.data[,c("cluster","Intermediate")],outputprefix = "Tasic2018.predictions.")


#Explore zeisel Rank4 annotations:
head(Hook2018gene@meta.data[, 10:48])
pheatmap::pheatmap(Hook2018gene@meta.data[, 10:48],show_rownames = F)

#Rank3 annotation
Hook2018gene.r3 <- Hook2018gene
Hook2018gene.r3@meta.data <- Hook2018gene@meta.data[, 1:9]

Hook2018gene.r3 <- CellTyper2(SeuratObject = Hook2018gene.r3,
                           model = zeisel.rank3.rfcv,
                           priorLabels = Hook2018gene.r3@meta.data$region,
                           outputFilename = paste("Hook2018gene.r3",".rfcv.predictions.w-singleZeisel-rank3",sep = ""))
head(Hook2018gene.r3@meta.data)
CrossCheck(PriorPostTable = Hook2018gene.r3@meta.data[,c("subset","Intermediate")], outputprefix = "Hook2018gene.r3.predictions.")

#Check Hook annotation with FPKM values from GEO
Hook2018fpkm <- get(load("~/data/Hook2018fpkm.seurat.Robj"))
dim(Hook2018fpkm@meta.data)
head(Hook2018fpkm@meta.data)

Hook2018fpkm <- CellTyper2(SeuratObject = Hook2018fpkm, model = zeisel.rank4.rfcv, priorLabels = Hook2018fpkm@meta.data$subset, outputFilename = paste("Hook2018fpkm",".rfcv.predictions.w-singleZeisel-rank4",sep = ""))
head(Hook2018fpkm@meta.data)

pheatmap::pheatmap(Hook2018fpkm@meta.data[, 12:50],show_rownames = F)

CrossCheck(PriorPostTable = Hook2018fpkm@meta.data[,c("subset","Intermediate")],outputprefix = "Hook2018fpkm.predictions.")



VlnPlot(object = Hook2018fpkm,x.lab.rot = 45, features = c("Hoxc4"),nCol = 1)
VlnPlot(object = Hook2018fpkm,x.lab.rot = 45, features = c("Lhx9"),nCol = 1)
VlnPlot(object = Hook2018fpkm,x.lab.rot = 45, features = c("En1"),nCol = 1)



#Annotate cells with scmap:
library(SingleCellExperiment)
library(scmap)

#Import training data:
zeisel.rank4.sub47 <- get(load("~/data/zeisel.rank4.sub47.seurat.Robj"))
ann <- as.data.frame(zeisel.rank4.sub47@meta.data[,c("TaxonomyRank4")])
rownames(ann) <- rownames(zeisel.rank4.sub47@meta.data)
colnames(ann) <- "cell_type1"
head(ann)
yan <- as.matrix(zeisel.rank4.sub47@data)
yan[1:3, 1:3]

sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(yan)), colData = ann)
logcounts(sce) <- normcounts(sce) ### Skipping this singe my data is in log2 already.
rowData(sce)$feature_symbol <- rownames(sce)
# remove features with duplicated names
sce <- sce[!duplicated(rownames(sce)), ]
sce

#Select features:
sce <- selectFeatures(sce, suppress_plot = FALSE)
table(rowData(sce)$scmap_features)
sce <- indexCluster(sce)
head(metadata(sce)$scmap_cluster_index)
heatmap(as.matrix(metadata(sce)$scmap_cluster_index))

Hook2018fpkm.sce <- as.SingleCellExperiment(from = Hook2018fpkm)
rowData(Hook2018fpkm.sce)$feature_symbol <- rownames(Hook2018fpkm.sce)
assays(Hook2018fpkm.sce)$logcounts <- as.matrix(assays(Hook2018fpkm.sce)$logcounts)
#Projection
scmapCluster_results <- scmapCluster(threshold = 0.5,
  projection = Hook2018fpkm.sce,
  index_list = list(
    yan = metadata(sce)$scmap_cluster_index
  )
)

plot(
  getSankey(
    colData(Hook2018fpkm.sce)$Prior,
    scmapCluster_results$scmap_cluster_labs[,'yan'],
    plot_height = 1400
  )
)

PriorPosttable <- data.frame(Prior=colData(Hook2018fpkm.sce)$Prior, Prediction=scmapCluster_results$scmap_cluster_labs[,'yan'])
head(PriorPosttable)
rownames(PriorPosttable) <- rownames(colData(Hook2018fpkm.sce))
CrossCheck(PriorPostTable = PriorPosttable, outputprefix = "output/Hook2sce")
save(PriorPosttable, file="~/data/PriorPosttable.Hook2scemap.Rdata")

Hook2018fpkm@meta.data$scemapPred <- scmapCluster_results$scmap_cluster_labs[,'yan']
head(Hook2018fpkm@meta.data)
save(Hook2018fpkm, file="~/data/Hook2018fpkm.seurat.Robj")

