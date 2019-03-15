.libPaths("/n/home13/yasinkaymaz/biotools/Rlibs/")
source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/functions.R")

geneinfo <- read.csv("~/data/mouse_VISp_2018-06-14_genes-rows.csv", row.names = 4, header = T)
mt.genes <- droplevels(geneinfo[which(geneinfo$chromosome == "MT"),]$gene_symbol)
r.genes <- droplevels(geneinfo[grep("ribosomal",x=geneinfo$gene_name,value = F),]$gene_symbol)

#exp <- read.delim("~/data/Hook_rsem_genes_TPM_matrix.txt",row.names = 1,header = T)
exp <- read.delim("~/data/Hook_rsem_genes_Count_matrix.txt",row.names = 1,header = T)

colnames(exp) <- stringr::str_split_fixed(string = colnames(exp),pattern = "\\.",n = 2 )[,1]
rownames(exp) <- make.names(stringr::str_split_fixed(string = rownames(exp),pattern = "_",n = 2 )[,2],unique=TRUE)
#Filter genes
gm.genes <- grep(pattern = "^Gm", x = rownames(exp), value = TRUE)
genesTofilter <- c(gm.genes, mt.genes, r.genes)
exp <- exp[which(!rownames(exp) %in% genesTofilter),]
dim(exp)

meta.data <- read.delim("~/data/SRP126648_SraRunTable.txt", header=TRUE)
meta <- meta.data %>% select(BioSample, age, region, passed_qc) %>% filter(passed_qc == TRUE) %>% unique() %>% droplevels()
dim(meta)
head(meta)
rownames(meta) <- meta$BioSample

#Filter qc failed cells
exp2 <- droplevels(exp[,which(colnames(exp) %in% meta$BioSample)])
dim(exp2)
dim(meta)

Hook2018 <- SeuratWrapper(ExpData = exp2, perp = 5, ProjectLabel = "Hook2018", NewMeta = meta, Normalize = T, scale.only.var = F, PCs = 14, dump.files = F,min.cells = 20)
PCElbowPlot(object = Hook2018)


dim(Hook2018@data)
dim(Hook2018@meta.data)
zeisel.rank3.rfcv <- get(load("~/data/zeisel.rank3.rfcv.RF_model_notImproved.Robj"))
Hook2018 <- CellTyper2(SeuratObject = Hook2018, model = zeisel.rank3.rfcv, priorLabels = Hook2018@meta.data$region, outputFilename = paste("Hook2018",".rfcv.predictions.w-singleZeisel-rank3",sep = ""))
head(Hook2018@meta.data)



#ISOFORM
#Prepare isoform expression data. from Rsem-Star
iso.exp <- read.delim("~/data/Hook_rsem_isoforms_TPM_matrix.txt",row.names = 1,header = T)
iso.exp[1:10,1:3]
dim(iso.exp)
colnames(iso.exp) <- stringr::str_split_fixed(string = colnames(iso.exp),pattern = "\\.",n = 2 )[,1]

grep(pattern = "Grm4", x = rownames(iso.exp), value = TRUE)

rownames(iso.exp) <- make.names(stringr::str_split_fixed(string = rownames(iso.exp),pattern = "_",n = 2 )[,2],unique=TRUE)

iso.exp2 <- droplevels(iso.exp[,which(colnames(iso.exp) %in% meta$BioSample)])

gm.genes <- grep(pattern = "^Gm", x = rownames(iso.exp), value = TRUE)
genesTofilter <- c(gm.genes)
iso.exp2 <- iso.exp2[which(!rownames(iso.exp2) %in% genesTofilter),]

meta.for.iso <- cbind(meta,Hook2018@meta.data[,c("res.1","BestVotesPercent","Intermediate")])

Hook2018iso <- SeuratWrapper(ExpData = iso.exp2, perp = 5, ProjectLabel = "Hook2018", NewMeta = meta.for.iso, Normalize = F, scale.only.var = F, PCs = 10, dump.files = F,min.cells = 20)
head(Hook2018iso@meta.data)

Hook2018iso <- SetAllIdent(Hook2018iso, id = "Intermediate")
transcripts <- grep(pattern = "^Gria1", x = rownames(Hook2018iso@data), value = TRUE)
pdf("Isoplots.pdf",width = 20,height = 30)
VlnPlot(object = Hook2018iso, features.plot = transcripts, use.raw = TRUE, y.log = TRUE,nCol = 1,x.lab.rot = T)
dev.off()



PCAPlot(object = Hook2018,group.by="age",pt.size=4)
PCAPlot(object = Hook2018,group.by="region",pt.size=2)
PCAPlot(object = Hook2018,group.by="res.1")
VizPCA(object = Hook2018, pcs.use = 1:2)

TSNEPlot(object = Hook2018,group.by="res.1",pt.size=4)
TSNEPlot(object = Hook2018,group.by="age")
TSNEPlot(object = Hook2018,group.by="region")

TSNEPlot(object = Hook2018,group.by="Intermediate",pt.size=4)
VlnPlot(object = Hook2018, features.plot = c("Mapre3","Mpl"), use.raw = TRUE, y.log = TRUE,nCol = 1,x.lab.rot = T)


