require(reshape2)

#Prepare isoform expression data. from Rsem-Star
iso.exp <- read.delim("~/data/Hook_rsem_isoforms_Count_matrix.txt",row.names = 1,header = T)
colnames(iso.exp) <- str_replace(colnames(iso.exp), ".isoforms.results", "")
iso.exp[1:10,1:3]
dim(iso.exp)



gds <- GEOquery::getGEO(filename = '~/data/GSE108020_series_matrix.txt.gz')
meta <- data.frame(row.names = gds@phenoData@data[["geo_accession"]],
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
iso.exp <- iso.exp[,which(colnames(iso.exp) %in% cells.pass)]

dim(meta)
dim(iso.exp)


Hook2018iso <- SeuratWrapper(ExpData = iso.exp, perp = 10, ProjectLabel = "Hook2018iso", NewMeta = meta, Normalize = T, scale.only.var = F, PCs = 5, dump.files = F, min.cells = 0)
head(Hook2018iso@meta.data)
save(Hook2018iso, file="~/data/Hook2018iso.seurat.Robj")

#Percent isoform usage table for all genes and cells
pct.iso <- read.delim("~/data/Hook_rsem_isoforms_Percent_matrix.txt",row.names = 1,header = T)
colnames(pct.iso) <- str_replace(colnames(pct.iso), ".isoforms.results", "")
pct.iso[1:10,1:3]

pct.iso <- pct.iso[,which(colnames(pct.iso) %in% cells.pass)]
dim(pct.iso)
save(pct.iso, file="~/data/Hook_rsem_isoforms_Percent_matrix-Sub.Rdata")


#Flip - Flop isoforms
ensids <- c("ENSMUST00000094179.10",
            "ENSMUST00000036315.15",
            "ENSMUST00000075316.9",
            "ENSMUST00000107745.7",
            "ENSMUST00000165288.1",
            "ENSMUST00000076349.11",
            "ENSMUST00000027020.12",
            "ENSMUST00000063508.14")

transcripts <- NULL
for (id in ensids){
  print(id)
  transcripts <- c(transcripts,grep(pattern = id, x = rownames(Hook2018iso@data), value = TRUE))
  print(transcripts)
}

plotdata <- cbind(Hook2018iso@meta.data,
                  t(as.matrix(Hook2018iso@data)[transcripts,]),
                  t(as.matrix(Hook2018gene@data)[c("Gria1","Gria2","Gria3","Gria4"),]))
head(plotdata)

pdf("output/FlipFlop-gene-isoforms.pdf",width = 10,height = 5)
#Plot Gene expressions
plotdata[,c("subset","Gria1","Gria2","Gria3","Gria4")] %>%
  melt() %>%
  ggplot(aes(x=subset, y=value, fill=subset ))+
  geom_boxplot(aes(fill=variable),notch=FALSE,outlier.colour="red")+
  labs(y="Normalized Expression")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_fill_discrete(name = "Genes")+
  theme(legend.position="top")

plotdata[,c("subset","ENSMUST00000094179.10_Gria1-202","ENSMUST00000036315.15_Gria1-201")] %>%
  melt() %>%
  ggplot(aes(x=subset, y=value, fill=subset ))+
  geom_boxplot(aes(fill=variable),notch=FALSE,outlier.colour="red")+
  labs(y="Normalized Expression")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_fill_discrete(name = "Isoforms")+
  theme(legend.position="top")

plotdata[,c("subset","ENSMUST00000075316.9_Gria2-201","ENSMUST00000107745.7_Gria2-202")] %>%
  melt() %>%
  ggplot(aes(x=subset, y=value, fill=subset ))+
  geom_boxplot(aes(fill=variable),notch=FALSE,outlier.colour="red")+
  labs(y="Normalized Expression")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_fill_discrete(name = "Isoforms")+
  theme(legend.position="top")

plotdata[,c("subset","ENSMUST00000165288.1_Gria3-209","ENSMUST00000076349.11_Gria3-201")] %>%
  melt() %>%
  ggplot(aes(x=subset, y=value, fill=subset ))+
  geom_boxplot(aes(fill=variable),notch=FALSE,outlier.colour="red")+
  labs(y="Normalized Expression")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_fill_discrete(name = "Isoforms")+
  theme(legend.position="top")

plotdata[,c("subset","ENSMUST00000027020.12_Gria4-201","ENSMUST00000063508.14_Gria4-202")] %>%
  melt() %>%
  ggplot(aes(x=subset, y=value, fill=subset ))+
  geom_boxplot(aes(fill=variable),notch=FALSE,outlier.colour="red")+
  labs(y="Normalized Expression")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_fill_discrete(name = "Isoforms")+
  theme(legend.position="top")
dev.off()

pheatmap::pheatmap(pct.iso[c(grep(pattern = "Gria1", x = rownames(Hook2018iso@data), value = TRUE)),],
                   cluster_rows = F,height = 15,
                   show_colnames = F,cellheight = 10,
                   annotation_col = plotdata[,c("age", "region", "Gria1")]
)

pdf("output/Heatmap-Percent-isoforms.pdf",width = 12,height = 4)
pheatmap::pheatmap(pct.iso[c(grep(pattern = "Gria1", x = rownames(Hook2018iso@data), value = TRUE)),],
                   cluster_rows = F,show_colnames = F,cellheight = 10,
                   annotation_col = plotdata[,c("age", "region", "Gria1")]
                   )
pheatmap::pheatmap(pct.iso[c(grep(pattern = "Gria2", x = rownames(Hook2018iso@data), value = TRUE)),],
                   cluster_rows = F,show_colnames = F,cellheight = 10,
                   annotation_col = plotdata[,c("age", "region", "Gria2")]
)
pheatmap::pheatmap(pct.iso[c(grep(pattern = "Gria3", x = rownames(Hook2018iso@data), value = TRUE)),],
                   cluster_rows = F,show_colnames = F,cellheight = 10,
                   annotation_col = plotdata[,c("age", "region", "Gria3")]
)
pheatmap::pheatmap(pct.iso[c(grep(pattern = "Gria4", x = rownames(Hook2018iso@data), value = TRUE)),],
                   cluster_rows = F,show_colnames = F,cellheight = 10,
                   annotation_col = plotdata[,c("age", "region", "Gria4")]
)
dev.off()






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








