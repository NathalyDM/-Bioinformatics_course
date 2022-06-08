library(BiocManager)
library(GEOquery)
library(limma)
library(affy)
library(DOSE)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(GOSemSim)
library(dplyr)
library(ggplot2)

input_ADvsCTL <- read.csv("alzheimer.csv", header=TRUE)
truefalse<-input_ADvsCTL$sig=='Up regulated'
load('alzheimer.rdata')

enrichment <- function(input,truefalse){
  input$total_contrast<-input$sig
  df_d7vs0<-input
  up<-na.omit(df_d7vs0$genes[truefalse])
  universe_genes<-na.omit(df_d7vs0$genes)

  up <- mapIds(org.Hs.eg.db, keys = as.character(up) , keytype = "SYMBOL", column="ENSEMBL")
  universe_genes2 <- mapIds(org.Hs.eg.db, keys = as.character(universe_genes) , keytype = "SYMBOL", column="ENSEMBL")

  up<-as.vector(na.omit(up))
  universe_genes<-as.character(as.vector(na.omit(universe_genes2)))

  ids<-bitr(up, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
  ids_un<-bitr(universe_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
  dedup_ids_un = ids[!duplicated(ids[c("ENSEMBL")]),]

  geneList2<-df_d7vs0$log2FoldChange[truefalse]
  names(geneList2)<-df_d7vs0$genes[truefalse]
  names(geneList2)<-mapIds(org.Hs.eg.db, keys = names(geneList2) , keytype = "SYMBOL", column="ENSEMBL")
  nombres<-bitr(names(geneList2), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
  names(geneList2)<-nombres$ENTREZID
  enriquecimiento<-c()

  ego <- enrichGO(gene          = up,
                  universe      = universe_genes,
                  ont = 'BP',
                  keyType       = "ENSEMBL",
                  OrgDb         = org.Hs.eg.db,
                  minGSSize     = 5,
                  maxGSSize     = 500,
                  pvalueCutoff  = 0.05,
                  readable      = TRUE)


  enriquecimiento$ego<-ego

  kk <- enrichKEGG(gene         = ids$ENTREZID,
                   organism     = 'hsa',
                   pvalueCutoff = 0.05)

  enriquecimiento$kk<-kk

  x <- enrichDO(gene          = ids$ENTREZID,
                universe      = ids_un$ENTREZID,
                minGSSize     = 5,
                maxGSSize     = 500,
                qvalueCutoff  = 0.05,
                readable      = FALSE)

  enriquecimiento$x<-x
  enriquecimiento$geneList2<-geneList2

  return(enriquecimiento)
}

significance <- function(input_ADvsCTL,p,logtrh){
  input_ADvsCTL$sig<-"No significant"
  for (i in 1:nrow(input_ADvsCTL)){
    if (input_ADvsCTL$log2FoldChange[i]>=logtrh & input_ADvsCTL$pvalue[i]>=log10(p)){
      input_ADvsCTL$sig[i]="Up regulated"
    }
    if (input_ADvsCTL$log2FoldChange[i]<=(-logtrh) & input_ADvsCTL$pvalue[i]>=-log10(p)){
      input_ADvsCTL$sig[i]="Down regulated"
    }
  }
  return(input_ADvsCTL)
}


#=================== Analyisis ========================
enriquecimiento<-enrichment(input_ADvsCTL,truefalse)
edox <- setReadable(enriquecimiento$x, 'org.Hs.eg.db', keyType='ENTREZID')

emapplot(edox,showCategory = 20, categorySize="pvalue", foldChange=enriquecimiento$geneList2)


#barplot(enriquecimiento$kk, showCategory=20)
venn<-upsetplot(enriquecimiento$x,n=20)
View(venn[["New_data"]])
df<-venn[["New_data"]]
vector<-df$`Alzheimer's disease`+df$dementia
gene_name<-mapIds(org.Hs.eg.db, keys = as.character(rownames(df[vector>0,])) , keytype = "ENTREZID", column="SYMBOL")


truefalse<-input_ADvsCTL$genes %in% as.vector(gene_name)
enriquecimiento2<-enrichment(input_ADvsCTL,truefalse)
edox2 <- setReadable(enriquecimiento2$x, 'org.Hs.eg.db', keyType='ENTREZID')
p1 <- cnetplot(edox2,showCategory = 10, categorySize="pvalue", foldChange=enriquecimiento$geneList2)
p1

#En caso no se coloree azul rojo
min.value = 0
max.value = 0.4
p1 + scale_colour_gradient2(name = "fold change",
                            low = "blue", mid = "blue", high = "red",
                            limits= c(min.value, max.value),
                            breaks=c(min.value , 0, max.value) )


df_names<-data.frame(nombre=as.vector(gene_name))
write.csv(df_names,'nombre_genes_dementia_alzheimer.csv')
