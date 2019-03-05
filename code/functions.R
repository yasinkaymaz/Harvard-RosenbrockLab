

DElimma <- function(SeuratObj ){
  library(limma)
  expr <- SeuratObj@data
  # Filter out genes that are 0 for every cell in this cluster
  bad <- which(rowSums(expr) == 0)
  expr <- expr[-bad,]
  meta <- SeuratObj@meta.data[,c("orig.ident", "subclass")]
  meta$subclass <- make.names(meta$subclass)
  head(meta)
  mm <- model.matrix(~0 + subclass, data = meta)

  fit <- lmFit(expr, mm)
  head(coef(fit)) # means in each sample for each gene
  conditions <- colnames(coef(fit))
  contr <- makeContrasts(paste(conditions[1]," - ",conditions[2],sep=""), levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contrasts = contr)
  tmp <- eBayes(tmp)
  topTable(tmp, sort.by = "P", n = 200)

  return(topTable(tmp, p.value = 0.01,lfc = 1,n=20000))
}

annotateGenes <- function(tableDE, geneinfoFile="~/data/mouse_ALM_2018-06-14_genes-rows.csv"){
  tableDE <- markers.L5IT.vs.Gaba_Glu
  geneinfo <- read.csv(geneinfoFile, row.names = 4, header = T)
  rownames(geneinfo) <- geneinfo$gene_symbol
  tableDE$Description <- as.character(geneinfo[rownames(tableDE),"gene_name"])

  return(tableDE)
}

FilterTable <- function(DEtable, allmarkers=F){
  if(allmarkers == T){
    FilteredDETable <- DEtable %>%
      select(-p_val)%>%
      filter(p_val_adj < 0.01) #%>%
      #filter(avg_logFC > 0.2) %>%
      #filter(pct.1 > 0.5) %>% filter( ((pct.1 - pct.2) > 0.25))
  }else{
  
  FilteredDETable <- DEtable %>%
    add_column(gene=rownames(DEtable)) %>%
    select(-p_val)%>%
    filter(p_val_adj < 0.01) #%>%
    #filter(avg_logFC > 0.2) %>%
    #filter(pct.1 > 0.5) %>% filter( ((pct.1 - pct.2) > 0.25))
  }
  return(FilteredDETable)
}

InteractiveTable <- function(DEtable){
  library(DT)
  DEtable %>%
    datatable( extensions = 'Scroller', options = list(
      deferRender = TRUE,
      scrollY = 200,
      scroller = TRUE
    ))
}


volcanoPlotly <- function(DEtable){
  library(plotly)
  p <- ggplot(DEtable, aes(avg_logFC, y=-1*log10(p_val_adj+10^-200), color=cluster ))+
    geom_point(aes(text=sprintf("Gene: %s", gene)))+
    ylab("-log10(padj)")+
    xlab("log Fold Change")
  ggplotly(p)
}
