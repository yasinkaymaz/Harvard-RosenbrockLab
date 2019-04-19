.libPaths("/n/home13/yasinkaymaz/biotools/Rlibs/")
source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/functions.R")

Hook2018gene <- get(load("~/data/Hook2018gene.seurat.Robj"))
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



