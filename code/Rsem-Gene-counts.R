#RSEM Gene count processing
geneinfo <- read.csv("~/data/mouse_VISp_2018-06-14_genes-rows.csv", row.names = 4, header = T)
mt.genes <- droplevels(geneinfo[which(geneinfo$chromosome == "MT"),]$gene_symbol)
r.genes <- droplevels(geneinfo[grep("ribosomal",x=geneinfo$gene_name,value = F),]$gene_symbol)

#count data from Rsem runs:
exp <- read.delim("~/data/Hook_rsem_genes_Count_matrix.txt",row.names = 1,header = T)
dim(exp)
colnames(exp) <- str_replace(colnames(exp), ".genes.results", "")
#get gene symbol from id
rownames(exp) <- make.names(stringr::str_split_fixed(string = rownames(exp),pattern = "_", n = 2 )[,2],unique = T)
exp[1:10,1:10]

#Filter genes
gm.genes <- grep(pattern = "^Gm", x = rownames(exp), value = TRUE)
genesTofilter <- c(gm.genes, mt.genes, r.genes)
exp <- exp[which(!rownames(exp) %in% genesTofilter),]
dim(exp)



gds <- GEOquery::getGEO(filename = '~/data/GSE108020_series_matrix.txt.gz')
meta <- data.frame(row.names = str_replace(gds@phenoData@data[["relation"]],"BioSample: https://www.ncbi.nlm.nih.gov/biosample/",""),
                   Sample_Name=gds@phenoData@data[["geo_accession"]],
                   qc=gds@phenoData@data[["passed_qc:ch1"]],
                   age=gds@phenoData@data[["age:ch1"]],
                   region=gds@phenoData@data[["region:ch1"]],
                   subset=gds@phenoData@data[["subset.cluster:ch1"]])

#Filter cells:
head(meta)
#Filter out low quality cell data as explained in the Hook et al.
meta <- meta[meta$qc == TRUE,]
#Exclude E15.5 cells as requested by CNSDR
meta <- meta[meta$age == "P7",]

cells.pass <- rownames(meta)
exp <- exp[,which(colnames(exp) %in% cells.pass)]

dim(meta)
dim(exp)

Hook2018gene <- SeuratWrapper(ExpData = exp, perp = 10, NewMeta = meta, ProjectLabel = "Hook2018gene", Normalize = T, scale.only.var = T, PCs = 5, dump.files = F, min.cells=5)
head(Hook2018gene@meta.data)
dim(Hook2018gene@data)
Hook2018gene <- RunUMAP(object = Hook2018gene)
save(Hook2018gene, file="~/data/Hook2018gene.seurat.Robj")

#Diagnostic Plots
DimPlot(Hook2018gene, reduction.use = "umap",group.by = "region")
DimPlot(Hook2018gene, reduction.use = "umap",group.by = "subset")

PCAPlot(object = Hook2018gene, group.by="region", pt.size=2)

DimHeatmap(object = Hook2018gene, dims = 1:15, balanced = TRUE)

Hook2018gene <- SetAllIdent(Hook2018gene, id="subset")
VlnPlot(object = Hook2018gene, x.lab.rot = 45, features = c("Th", "Lhx9","Slc6a3","Ldb2"),nCol = 2)
VlnPlot(object = Hook2018gene, x.lab.rot = 45, features = c("Crhr1", "Nsf", "Mapt"),nCol = 1)

