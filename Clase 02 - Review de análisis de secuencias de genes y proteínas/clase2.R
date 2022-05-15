colores1 <- c(rgb( 121/255, 159/255,167/255,alpha = 1),
              rgb( 206/255, 232/255,221/255,alpha = 1),
              rgb( 250/255, 240/255,230/255,alpha = 1),
              rgb( 245/255, 213/255,202/255,alpha = 1),
              rgb( 219/255, 154/255,150/255,alpha = 1)
)

colores2 <- c(rgb( 121/255, 159/255,167/255,alpha = .5),
              rgb( 206/255, 232/255,221/255,alpha = .5),
              rgb( 250/255, 240/255,230/255,alpha = .5),
              rgb( 245/255, 213/255,202/255,alpha = .5),
              rgb( 219/255, 154/255,150/255,alpha = .5)
)

## Instalando paquetes con dependencias

install.packages("seqinr", dependencies = T)
install.packages("ape", dependencies = T)
library(seqinr)
library(ape)

## 1. Cargando una secuencia manualmente ====

path <- "C:/Users/pitca/Documents/2021-1/Bioinformatica_Avanzada/Marcelo/bioinf_avanza/Clases/Clase2y3/data/sequence.fasta"
actin.necator <- seqinr::read.fasta(path, seqtype = "AA")


attr(actin.necator$ETN76030.1,"Annot")
tabla <- table(actin.necator$ETN76030.1)
sort(table(actin.necator$ETN76030.1))
sort(table(actin.necator$ETN76030.1),decreasing = T)
plot(sort(table(actin.necator$ETN76030.1),decreasing = T), col=colores1)
barplot(sort(table(actin.necator$ETN76030.1),decreasing = T), col=colores1)
longitud <- length(actin.necator$ETN76030.1)
freqre <- tabla/longitud
barplot(freqre)

#2. Cargando multiples secuencias ====

path <- "C:/Users/pitca/Documents/2021-1/Bioinformatica_Avanzada/Marcelo/bioinf_avanza/Clases/Clase2y3/data/app_Hs_multiple_seq.txt"
app.Hs <- seqinr::read.fasta(path, as.string = T, seqtype = "AA")

#3. Importando secuencias desde R ====

seqinr::choosebank(infobank = T)
seqinr::choosebank("swissprot")
app2 <- query("app2", "K=Amyloid AND sp=Homo sapiens")
summary(app2)
app2$call
app2$nelem
app2$req
app2$req[[2]]
app2$req[2]
seqinr::getSequence(app2$req[[2]])

#4. Usando ape ====

acs.num <- c("MT419843.1","MT450923.1","MT419850.1",
             "MT419838.1","MT419849.1","MT419840.1")
library(ape)
cov2 <- ape::read.GenBank(acs.num, as.character = T)
cov2

colores1 <- c(rgb( 121/255, 159/255,167/255,alpha = 1),
              rgb( 206/255, 232/255,221/255,alpha = 1),
              rgb( 250/255, 240/255,230/255,alpha = 1),
              rgb( 245/255, 213/255,202/255,alpha = 1),
              rgb( 219/255, 154/255,150/255,alpha = 1)
)

colores2 <- c(rgb( 121/255, 159/255,167/255,alpha = .5),
              rgb( 206/255, 232/255,221/255,alpha = .5),
              rgb( 250/255, 240/255,230/255,alpha = .5),
              rgb( 245/255, 213/255,202/255,alpha = .5),
              rgb( 219/255, 154/255,150/255,alpha = .5)
)
acs.num <- c("MT419843.1","MT450923.1","MT419850.1",
             "MT419838.1","MT419849.1","MT419840.1")
library(ape)
library(seqinr)
cov2 <- ape::read.GenBank(acs.num, as.character = T)
names(cov2)
summary(cov2)
cov2[[1]]
round(table(cov2$MT419843.1)*100/length(cov2$MT419843.1),1)

#Crear pdf con los graficos
pdf(file = paste0(getwd(),"/comparacion.pdf"), paper = "US")
par(mfrow = c(3,2))

for (gen in 1:length(cov2)) {
  print(table(cov2[gen]))
  barplot(round(table(cov2[[gen]])*100/length(cov2[[gen]]), 1), main = names(cov2)[gen], las = 1, col=colores1)
}
dev.off()

#Indice

#Analisis estadistico ====
#analisis entre probabilidad estocastica que suceda dos nucleotiodos/entre la verdadera presencia del gen
seqinr::rho(cov2$MT419843.1) 
plot(seqinr::rho(cov2$MT419843.1)-1,ylim = c(-0.5,0.5))
#Otra forma -> escalar cualquier rango de datos a una escala conocida (Z)
zbase <- seqinr::zscore(sequence = cov2$MT419843.1, modele = "base", exact = T) 
zbase
plot(zbase)
#rho pero con trinucleotidos
seqinr::rho(cov2$MT419843.1,wordsize = 3)
plot(seqinr::rho(cov2$MT419843.1,wordsize = 3)-1,ylim = c(-0.5,0.5))

seqinr::GC(cov2[[1]]) #Porcentaje GC de la secuencia
dev.new(name = "GC%") #Generar nueva ventana
par(mfrow=c(2,2)) #Dividir pantalla de grafico
#Porcentaje de GC en grupos (win.size)
for (win.size in c(20,50,100,500)){
  gc.perc<-vector()
  k<-1
  for (i in seq(from=1,to=length(cov2[[1]])-win.size,by=win.size/10)){
    j<-i+win.size-1
    gc.perc[k]<- seqinr::GC(cov2[[1]][i:j])
    k<-k+1
  }
  gc.perc <- gc.perc[!is.na(gc.perc)]
  plot(gc.perc*100,type="l", ylim = c(0,100),
       main = win.size, ylab="GC%", xlab="", col = colores1[5])
}


#Instalar GenomeGraph
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
install.packages("GenomeGraphs", dependencies = T)

library(biomaRt)
library(GenomeGraphs)

# http://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?g=ENSG00000142192;r=21:25880550-26170770;t=ENST00000346798
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
listDatasets(mart)
#Descargar el gen desde la base de datos
gene <- makeGene(id = "ENSG00000136244",type = "ensembl_gene_id",
                 biomart = mart)
#Muestra intrones(amarillo) y exones (negro)
gdPlot(gene)
#Descargar el transcrito
transcript <- makeTranscript(id = "ENSG00000136244",type =
                               "ensembl_gene_id", biomart = mart)
#Grafica gen con su transcrito
gdPlot(list(gene, transcript)) #numero isoformas = numero de lineas

gdPlot(list(makeTitle("IL6"),
            makeIdeogram(chromosome = 7), gene, transcript,
            makeGenomeAxis()))

#Agregas anotaciones para se?alar cierta regiones de un gen
customann<-makeAnnotationTrack(start=c(25980550,26020650,26080850),
                               end=c(26000550,26025650,26085850),
                               feature=c('bind','bind','del'), dp=DisplayPars(bind = 'gray',del='black'))
gdPlot(list(makeTitle("Human amyloid beta precursor gene"),
            makeIdeogram(chromosome = 21), gene, transcript,
            customann, makeGenomeAxis()))
