
source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/functions.R")
source("code/functions.R")

geneinfo <- read.csv("~/data/mouse_VISp_2018-06-14_genes-rows.csv", row.names = 4, header = T)
mt.genes <- droplevels(geneinfo[which(geneinfo$chromosome == "MT"),]$gene_symbol)
r.genes <- droplevels(geneinfo[grep("ribosomal",x=geneinfo$gene_name,value = F),]$gene_symbol)

#Hook 2018
#FPKM data from GEO
exp.fpkm <- read.delim("~/data/GSE108020_fpkm_table.txt",row.names = 1,header = T)
dim(exp.fpkm)#473 cells
gds <- GEOquery::getGEO(filename = '~/data/GSE108020_series_matrix.txt.gz')

meta <- data.frame(row.names = str_replace(gds@phenoData@data[["relation"]],"BioSample: https://www.ncbi.nlm.nih.gov/biosample/",""),
                   SAMid = str_replace(gds@phenoData@data[["relation"]],"BioSample: https://www.ncbi.nlm.nih.gov/biosample/",""),
                   sampnames = paste(gds@phenoData@data[["project:ch1"]],"_0",gds@phenoData@data[["correct_source_plate:ch1"]],"_",gds@phenoData@data[["well:ch1"]], sep = ""),
                   Sample_Name=gds@phenoData@data[["geo_accession"]],
                   qc=gds@phenoData@data[["passed_qc:ch1"]],
                   age=gds@phenoData@data[["age:ch1"]],
                   region=gds@phenoData@data[["region:ch1"]],
                   subset=gds@phenoData@data[["subset.cluster:ch1"]])

head(meta)
identical(as.character(meta$sampnames), colnames(exp.fpkm))

colnames(exp.fpkm) <- meta$SAMid

#Filter cells:
head(meta)
#Filter out low quality cell data as explained in the Hook et al.
meta <- meta[meta$qc == TRUE,]
#Exclude E15.5 cells as requested by CNSDR
meta <- meta[meta$age == "P7",]

cells.pass <- rownames(meta)
exp.fpkm <- exp.fpkm[,which(colnames(exp.fpkm) %in% cells.pass)]
dim(meta)
dim(exp.fpkm)
exp.fpkm[1:3,1:3]
#convert Ensembl ids to gene symbol ids:
genetable <- data.frame(row.names = rownames(exp.fpkm),
                        ensembl_gene_id=stringr::str_split_fixed(string = rownames(exp.fpkm),pattern = "\\.",n = 2 )[,1],
                        ensembl.id = rownames(exp.fpkm))

head(genetable)

ensembl <- biomaRt::useMart("ensembl", dataset="mmusculus_gene_ensembl")
exp.genedatabase <- biomaRt::getBM(filters= "ensembl_gene_id",uniqueRows = TRUE, attributes= c("ensembl_gene_id","mgi_symbol", "description"),values=genetable$ensembl_gene_id,mart= ensembl)
head(exp.genedatabase)
genetable2 <- merge(x = genetable, y = exp.genedatabase, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")

head(genetable2)
dim(exp.fpkm)
exp.fpkm <- exp.fpkm[which(rownames(exp.fpkm) %in% genetable2$ensembl.id),]
identical(as.character(rownames(exp.fpkm)), as.character(genetable2$ensembl.id))
rownames(exp.fpkm) <- make.names(genetable2$mgi_symbol,unique = T)



Hook2018fpkm <- SeuratWrapper(ExpData = exp.fpkm, perp = 10, NewMeta = meta, ProjectLabel = "Hook2018fpkm", Normalize = F, scale.only.var = T, PCs = 5, dump.files = F, min.cells=0)
head(Hook2018fpkm@meta.data)
dim(Hook2018fpkm@data)
Hook2018fpkm <- RunUMAP(object = Hook2018fpkm)
save(Hook2018fpkm, file="~/data/Hook2018fpkm.seurat.Robj")

Hook2018fpkm <- SetAllIdent(Hook2018fpkm, id="subset")
VlnPlot(object = Hook2018fpkm,x.lab.rot = 45, features = c("Th", "Lhx9","Slc6a3","Ldb2"),nCol = 2)
VlnPlot(object = Hook2018fpkm,x.lab.rot = 45, features = c("Crhr1", "Nsf", "Mapt"),nCol = 1)
VlnPlot(object = Hook2018fpkm,x.lab.rot = 45, features = c("Snap25"),nCol = 1)


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
CrossCheck(PriorPostTable = PriorPosttable, outputprefix = "output/Hook2sce")
save(PriorPosttable, file="~/data/PriorPosttable.Hook2scemap.Rdata")

Hook2018fpkm@meta.data$scemapPred <- scmapCluster_results$scmap_cluster_labs[,'yan']
head(Hook2018fpkm@meta.data)
save(Hook2018fpkm, file="~/data/Hook2018fpkm.seurat.Robj")
