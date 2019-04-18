
source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/functions.R")


source("code/functions.R")

geneinfo <- read.csv("~/data/mouse_VISp_2018-06-14_genes-rows.csv", row.names = 4, header = T)
mt.genes <- droplevels(geneinfo[which(geneinfo$chromosome == "MT"),]$gene_symbol)
r.genes <- droplevels(geneinfo[grep("ribosomal",x=geneinfo$gene_name,value = F),]$gene_symbol)

#Hook 2018
#FPKM data from GEO
exp <- read.delim("~/data/GSE108020_fpkm_table.txt",row.names = 1,header = T)
dim(exp)#473 cells
gds <- GEOquery::getGEO(filename = '~/data/GSE108020_series_matrix.txt.gz')
meta <- data.frame(row.names = colnames(exp),
                  Sample_Name=gds@phenoData@data[["geo_accession"]],
                  qc=gds@phenoData@data[["passed_qc:ch1"]],
                  age=gds@phenoData@data[["age:ch1"]],
                  region=gds@phenoData@data[["region:ch1"]],
                  subset=gds@phenoData@data[["subset.cluster:ch1"]])
head(meta)
setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("name", "age", "gender"))#usefull oneliner

#convert Ensembl ids to gene symbol ids:
genetable <- data.frame(row.names = rownames(exp),
                        ensembl_gene_id=stringr::str_split_fixed(string = rownames(exp),pattern = "\\.",n = 2 )[,1],
                        ensembl.id = rownames(exp))

head(genetable)

ensembl <- biomaRt::useMart("ensembl", dataset="mmusculus_gene_ensembl")

exp.genedatabase <- biomaRt::getBM(filters= "ensembl_gene_id",uniqueRows = TRUE, attributes= c("ensembl_gene_id","mgi_symbol", "description"),values=genetable$ensembl_gene_id,mart= ensembl)
missing.ids <- setdiff( rownames(genetable),exp.genedatabase$ensembl_gene_id)
duplicates <- names(table(exp.genedatabase$mgi_symbol)[table(exp.genedatabase$mgi_symbol)>1])
#exp.genedatabase[which(exp.genedatabase$mgi_symbol %in% duplicates),]
exp.genedatabase <- exp.genedatabase %>% mutate(mgi_symbol = if_else((mgi_symbol %in% duplicates),ensembl_gene_id, as.character(mgi_symbol) ))
rownames(exp.genedatabase) <- make.names(unique = T, names = exp.genedatabase$ensembl_gene_id)

test <- genetable %>% left_join(., exp.genedatabase, by= "ensembl_gene_id") %>% mutate(mgi_symbol = if_else(is.na(mgi_symbol),ensembl_gene_id, as.character(mgi_symbol) )) %>% filter(mgi_symbol != "4930553J12Rik")
rownames(test) <- test$ensembl.id
identical(rownames(test),rownames(exp))

rownames(exp) <- test$mgi_symbol
exp[1:10,1:3]

#Filter genes
gm.genes <- grep(pattern = "^Gm", x = rownames(exp), value = TRUE)
genesTofilter <- c(gm.genes, mt.genes, r.genes)
length(genesTofilter)

exp <- exp[which(!rownames(exp) %in% genesTofilter),]
dim(exp)

#Filter cells:
head(meta)
cells.pass <- rownames(meta[which(meta$qc == TRUE),])
exp <- exp[,which(colnames(exp) %in% cells.pass)]
meta <- meta[meta$qc == TRUE,]

Hook2018 <- SeuratWrapper(ExpData = exp, NewMeta = meta, ProjectLabel = "Hook2018", Normalize = F, scale.only.var = T, PCs = 20, dump.files = F,min.cells=20)
head(Hook2018@meta.data)

PCAPlot(object = Hook2018,group.by="age",pt.size=2)
PCAPlot(object = Hook2018,group.by="region",pt.size=2)
PCAPlot(object = Hook2018,group.by="res.1")
VizPCA(object = Hook2018, pcs.use = 1:2)

TSNEPlot(object = Hook2018,group.by="res.1",pt.size=2)
TSNEPlot(object = Hook2018,group.by="age")
TSNEPlot(object = Hook2018,group.by="region")
TSNEPlot(object = Hook2018,group.by="subset")


#Annotate cell types with KI-superset
zeisel.rank3.rfcv <- get(load("~/data/zeisel.rank3.rfcv.RF_model_notImproved.Robj"))
Hook2018 <- CellTyper2(SeuratObject = Hook2018, model = zeisel.rank3.rfcv, priorLabels = Hook2018@meta.data$region, outputFilename = paste("Hook2018",".rfcv.predictions.w-singleZeisel-rank3",sep = ""))
head(Hook2018@meta.data)

TSNEPlot(object = Hook2018,group.by="Intermediate")

pdf("Genes-Violins.pdf",width = 20,height = 30)
Hook2018 <- SetAllIdent(Hook2018, id = "subset")
VlnPlot(object = Hook2018, features.plot = c("Gria1","Gria4","Grm4", "Gpr83"), use.raw = TRUE, y.log = TRUE,nCol = 1,x.lab.rot = T)
dev.off()


#Prepare isoform expression data. from Rsem-Star
iso.exp <- read.delim("~/data/Hook_rsem_isoforms_TPM_matrix.txt",row.names = 1,header = T)
iso.exp[1:10,1:3]
dim(iso.exp)
colnames(iso.exp) <- stringr::str_split_fixed(string = colnames(iso.exp),pattern = "\\.",n = 2 )[,1]



gds <- GEOquery::getGEO(filename = '~/data/GSE108020_series_matrix.txt.gz')
meta <- data.frame(row.names = gds@phenoData@data[["geo_accession"]],
                   Sample_Name=gds@phenoData@data[["geo_accession"]],
                   qc=gds@phenoData@data[["passed_qc:ch1"]],
                   age=gds@phenoData@data[["age:ch1"]],
                   region=gds@phenoData@data[["region:ch1"]],
                   subset=gds@phenoData@data[["subset.cluster:ch1"]])
meta <- meta[meta$qc == TRUE,]
dim(meta)
head(meta)

meta.sra <- read.delim("~/data/SRP126648_SraRunTable.txt", header=TRUE)
meta.merged <- meta.sra %>% select(Sample_Name, BioSample, passed_qc) %>% filter(passed_qc == TRUE) %>% select(-passed_qc) %>% unique() %>% droplevels() %>% left_join(., meta, by="Sample_Name")
rownames(meta.merged) <- meta.merged$BioSample
head(meta.merged)

iso.exp <- droplevels(iso.exp[,which(colnames(iso.exp) %in% meta.merged$BioSample)])
dim(iso.exp)



Hook2018iso <- SeuratWrapper(ExpData = iso.exp, perp = 5, ProjectLabel = "Hook2018iso", NewMeta = meta.merged, Normalize = F, scale.only.var = F, PCs = 10, dump.files = F,min.cells = 20)
head(Hook2018iso@meta.data)



Hook2018iso <- SetAllIdent(Hook2018iso, id = "subset")
transcripts <- grep(pattern = "Gria1", x = rownames(Hook2018iso@data), value = TRUE)
pdf("output/Gria1.Isoplots.pdf",width = 20,height = 15)
VlnPlot(object = Hook2018iso, features.plot = transcripts, use.raw = TRUE, y.log = TRUE,nCol = 1,x.lab.rot = T)
dev.off()


transcripts <- grep(pattern = "Gria4", x = rownames(Hook2018iso@data), value = TRUE)
pdf("output/Gria4.Isoplots.pdf",width = 20,height = 20)
VlnPlot(object = Hook2018iso, features.plot = transcripts, use.raw = TRUE, y.log = TRUE,nCol = 1,x.lab.rot = T)
dev.off()


transcripts <- grep(pattern = "Grm4", x = rownames(Hook2018iso@data), value = TRUE)
pdf("output/Grm4.Isoplots.pdf",width = 20,height = 30)
VlnPlot(object = Hook2018iso, features.plot = transcripts, use.raw = TRUE, y.log = TRUE,nCol = 1,x.lab.rot = T)
dev.off()

transcripts <- grep(pattern = "Gpr83", x = rownames(Hook2018iso@data), value = TRUE)
pdf("output/Gpr83.Isoplots.pdf",width = 20,height = 10)
VlnPlot(object = Hook2018iso, features.plot = transcripts, use.raw = TRUE, y.log = TRUE,nCol = 1,x.lab.rot = T)
dev.off()








meta.data <- read.delim("~/data/SRP126648_SraRunTable.txt", header=TRUE)
meta <- meta.data %>% select(BioSample, age, region, passed_qc) %>% filter(passed_qc == TRUE) %>% unique() %>% droplevels()
dim(meta)
head(meta)
rownames(meta) <- meta$BioSample
exp2 <- droplevels(exp[,which(colnames(exp) %in% meta$BioSample)])
dim(exp2)
dim(meta)



Hook2018 <- SeuratWrapper(ExpData = exp, ProjectLabel = "Hook2018", Normalize = T, scale.only.var = T, PCs = 20, dump.files = F)

head(Hook2018@meta.data)
TSNEPlot(Hook2018, group.by="res.1",do.label = T)


Gria1 <- "ENSMUSG00000020524.16"
Gria4 <- "ENSMUSG00000025892.15"
VlnPlot(object = Hook2018, features.plot = c("Gria1","Gria4"), use.raw = TRUE, y.log = TRUE,nCol = 1,x.lab.rot = T)





sampnames <- paste(gds@phenoData@data[["project:ch1"]],"_0",gds@phenoData@data[["correct_source_plate:ch1"]],"_",gds@phenoData@data[["well:ch1"]], sep = "")

exp.fpkm <- read.delim("~/data/GSE108020_fpkm_table.txt",row.names = 1,header = T)

identical(sampnames,colnames(exp.fpkm))



