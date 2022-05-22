
#Instalacion Biomart====
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")
BiocManager::install("biomaRt")
library(biomaRt)

#Grupo de datos en Biomart====

# Ver que base de datos hay disponibles
biomaRt::listMarts()

# Conexión con una base de datos y ver su grupo de datos del mismo (subset)
ensembl <- biomaRt::useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL")
datasets <- biomaRt::listDatasets(ensembl)
head(datasets)
dim(datasets)
str(datasets)

#Para encontrar donde se encuentra el registro correspondinete a " Homo sapiens " (hsapiens_gene_ensembl), podemos usar la función str_detect() del paquete "stringr".
ensembl <- useDataset(mart=ensembl, dataset = "hsapiens_gene_ensembl") #Otra manera de hacerlo
install.packages("stringr", dependencies = T)
library(stringr)
datasets[stringr::str_detect(string = "hsapiens_gene_ensembl",pattern = as.character(datasets$dataset) ),]

#Para seleccinar la base de datos de seres humanos, podemos resescribir nuestro objeto "ensembl" usando la función useDataset():
ensembl <- biomaRt::useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
str(ensembl) #clase Mart = clase 4 creada

#Ahora el objeto "ensembl" generado es tipo lista, creado por el paquete (Objeto S4), que tiene 7 slots:

#@ biomart: vector tipo caracter que tiene el nombre del datasets ("ENSEMBL_MART_ENSEMBL").
#@ host: vector tipo caracter que tiene la dirección de donde se descargó el dataset "https://www.ensembl.org:443/biomart/martservice")
#@ vschema: vector tipo caracter que tiene el detalle de como fue importado el dataset ("default")
#@ version: vector tipo caracter que debe especificar la versión del dataset.
#@ dataset: vector tipo caracter que especifica el nombre del dataset ("hsapiens_gene_ensembl").
#@ filters: objeto tipo data.frame que tiene las varibles por las que podemos filtra nuestra información (432)
#@ attributes: objeto tipo data.frame que especifica los atributos de las variables (3061)

#Para darnos una idea del contenido de los filtros que podemos usar veamos una parte del data.frame "@filters"
ensembl@filters
head(ensembl@filters)[,1:2]

#El paquete Biomart nos ofrece una función que nos permite ver que filtros podríamos utilizar en nuestro dataset con la función "listFilters()".
filtros <- biomaRt::listFilters(ensembl)
head(filtros)
#Como verán esta función es muy parecida a la que usamos anteriormente. De la misma manera podemos explorar el data.frame "Attrributes"
head(ensembl@attributes)
#Y el paquete BiomaRt también nos ofrece una función para listar los atributos que podemos utilizar con nuestro dataset, con la función listAtrributes()
atributos <- biomaRt::listAttributes(ensembl)
head(atributos)

#La lista de atributos es muy extensa. Biomart agrupo estos atributos en grupos, usando la columna "page" como factor. Para saber cuales son las categorias disponibles para este factor, Biomart nos ofrece la función attributePages
biomaRt::attributePages(mart = ensembl)

#Vemos que se consideran 6 categorías. Nosotros podemos llamar a una de las categorías en específico y mostrar en pantalla la lista de atributos que están considerado dentro "snp" por ejemplo.
head(biomaRt::listAttributes(mart = ensembl, page = "snp"))
head(biomaRt::listAttributes(mart = ensembl, page = "sequences"))

#Hacer un query en Biomart====

#Para generar el "query" usamos la función getBM(), que tiene 3 parámetros:
  
  #attributes: Es un objeto vector con la lista de atributos que quiero importar
  #filters: Es un objeto vector con la lista de filtros para seleccionar mi "query".
  #values: Es un objeto vector con valores para los filtros.
  #mart: es el objeto Mart que vamos a analizar, en nuestro caso es ensembl.

affyids=c("202763_at","209310_s_at","207500_at")
getBM(attributes=c('affy_hg_u133_plus_2', 'entrezgene_id', 'external_gene_name'), 
      filters = 'affy_hg_u133_plus_2', 
      values = affyids, 
      mart = ensembl)

geneID = c("3064","9370","6855", "351", "7124")
getBM(attributes=c('affy_hg_u133_plus_2', 'entrezgene_id', 'external_gene_name'), 
      filters = 'entrezgene_id', 
      values = geneID, 
      mart = ensembl)

getBM(attributes=c('affy_hg_u133_plus_2', 'entrezgene_id', 'external_gene_name'), 
      filters = 'entrezgene_id', 
      values = geneID, 
      mart = ensembl)

#Busqueda para definir parámetros en los filtros, atributos y datasets====

#Cuando desconozcamos el nombre exacto del atributo, filtro o dataset; el paquete BiomaRt tiene funciones específicas para hacerlo como las funciones searchDatasets(), searchAtributes() y searchFilters(). Por ejemplo si queremos buscar una base de datos de Homo sapiens podemos usar la función searchDatasets():
biomaRt::searchDatasets(mart = ensembl, pattern = "hsapiens")
biomaRt::searchDatasets(mart = ensembl, pattern = "mmusculus")
biomaRt::searchDatasets(mart = ensembl, pattern = "drerio")
biomaRt::searchDatasets(mart = ensembl, pattern = "rnorvegicus")
biomaRt::searchDatasets(mart = ensembl, pattern = "ggallus")
biomaRt::searchDatasets(mart = ensembl, pattern = "btaurus")


#Para buscar en los atributos aquellos que estén relacionados con el patrón "hgnc", usamos la siguiente línea de comando:
biomaRt::searchAttributes(mart = ensembl, pattern = "hgnc")
#Si queremos encontrar los atributos relacionados a "affymetrics" podemos buscar aquellos que tengan el patrón "affy"
ATRI <- biomaRt::searchAttributes(mart = ensembl, pattern = "affy")
#Por ejemplo, si deseamos hacer una busqueda más fina usando dos parámetros, por ejemplo que el filtro contenga las palabras "ensembl" y la palabra "id" en el nombre lo podemos ejecutar la siguiente linea de código:
biomaRt::searchFilters(mart = ensembl, pattern = "ensembl.*id")
biomaRt::searchAttributes(mart = ensembl, pattern = "ensembl.*id")
biomaRt::searchFilters(mart = ensembl, pattern = "translation")

#Para definir que valores son válidos en el parámetro "chromosome_name" de los filtros podemos usar la función listFilterValues()
biomaRt::listFilterValues(mart = ensembl, filter = "chromosome_name")[20:40]
biomaRt::searchFilterValues(mart = ensembl, filter = "phenotype_description", pattern = "Crohn")

#Fijando resultados====

#Para ahorrar el tiempo en las busquedas que ya realizamos anteriormente, el paquete Biomart guarda dicha busqueda en nuestra computadora. Para saber la localización de dichos archivos usamos la fucnión biomartCacheInfo()
biomaRt::biomartCacheInfo()

#Asociación de terminos GO usando base de datos Biomart====

#Para hacer esta asociación debemos seguir tres sencillos pasos: 
  #Definir el dataset de Biomart correspondiente, para el ejemplo será el de ser humano. 
  #Cargar la lista de genes (de nuestro experimento) 
  #Hacere la asociación con la función getBM()

# Paso 1 definimos el dataset
ensembl <- biomaRt::useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL")
ensembl <- biomaRt::useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
# Paso 2 cargamos la base de datos
immunome <- read.csv("C:/Users/51990/OneDrive/Desktop/Bioinfo_avanzada/Clase4y5/human_gene_set.txt", sep="")
head(immunome)
# Paso 3 asociamos anotaciones
goannot <- biomaRt::getBM(attributes = c("entrezgene_id", "go_id", "name_1006"),
                          filters = "entrezgene_id", 
                          values = immunome,
                          mart = ensembl)

#La función getBM() genera un objeto lista, dentro del cual encontramos un data.frame con el número de columnas que nosotros indicamos en el parámetro attributes de la función getBM

typeof(goannot)
str(goannot)
nrow(goannot)
head(goannot)

#Si generamos una tabla de frecuencias para observar las 20 funciones GO más frecuentes, observamos que
cbind(head(sort(table(goannot$name_1006), decreasing = TRUE),20)/nrow(goannot))
cbind(head(sort(table(goannot$name_1006), decreasing = TRUE),20)/nrow(goannot)*100)

#Asociacion de terminos GO usando base de datos externa====

#Tambien podemos utilizar bases de datos externas para anotar los terminos GO a nuestra lista de genes, con la base de datos org.Hs.eg.db que está incorporada en el paquete AnnotationDbi. Para ello necesitamos cargar la librería:
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
AnnotationDbi::columns(org.Hs.eg.db)
#La base de edatos org.Hs.eg.db tiene las siguientes columnas que puede ser asociadas con nuesetra lista de genes. Si deseamos alguna definción específica de alguna columna podememos usar la función help() y el nombre de la columna.
help("SYMBOL")
#Para imprimir los valores de una columna podemeos usar la función keys() y especificar el nombre de la columna en el parámetro keytype =
keys(org.Hs.eg.db, keytype="GO")
cbind(name=head(keys(org.Hs.eg.db, keytype="GENENAME")),
      uniprot=head(keys(org.Hs.eg.db, keytype = "UNIPROT")))

cbind(name = head(keys(org.Hs.eg.db, keytype="GENENAME"),10),
      symbol = head(keys(org.Hs.eg.db, keytype="SYMBOL"),10),
      uniprot = head(keys(org.Hs.eg.db, keytype="UNIPROT"),10),
      path = head(keys(org.Hs.eg.db, keytype="PATH"),10),
      GO = head(keys(org.Hs.eg.db, keytype="GO"),10),
      entrezID = head(keys(org.Hs.eg.db, keytype="ENTREZID"),10)
      )
#Para poder asociar los téminos GO de la base de datos org.Hs.eg.db con la lista de nuestros genes, que está en el objeto inmunome, podemos hacerlo con la función select() del paquete Biomart. La función select() de Biomart tiene 4 parámetros:
  #x = que define la base de datos externa
  #keys = es el vector caracter que define los ID de los genes
  #columns = La columna de la base de datos externa que se va a asociar a los IDs
  #keytype = define el tipo de datos del parámetro key
head(immunome)
#La función que generaría una base de datos asociación el termino GO con los códigos EntrezID sería la siguiente.
goannot2 <- biomaRt::select(x = org.Hs.eg.db,
                            keys = as.character(immunome$EntrezGeneID), 
                            columns = "GO", 
                            keytype = "ENTREZID")
#Si nosotros quisieramos asociar el código de las rutas KEGG se utilizaría las siguientes líneas de código
pathannot <- biomaRt::select(x = org.Hs.eg.db,
                             keys = as.character(immunome$EntrezGeneID), 
                             columns = "PATH", 
                             keytype = "ENTREZID")

#Analisis de sobreexpresion====

#Para hacer estos análisis debemos cargar nuestra lista de genes de referencia, que para nuestro caso será "gen.ref" y tendrá una longitud de 176 registros
#gen.ref <- read.csv(paste0(dirname(getwd()),"/data/clase_3/human_gene_selected_set.txt"))
gen.ref <- read.csv("C:/Users/51990/OneDrive/Desktop/Bioinfo_avanzada/Clase4y5/human_gene_selected_set.txt", sep="")
head(gen.ref)
dim(immunome)
dim(gen.ref)

#Vamos a utilizar el paquete GOstat para lo cual necesitamos instalarlo por Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GOstats")

#Una vez instalado el paquete "GOstats", cargamos la librería a nuestro entorno de trabajo con la función library().
library("GOstats")
library("org.Hs.eg.db")

#Luego creamos un objeto GOHyperGParams para realizar los cálculos hipergeométricos.
paramsGO <- new("GOHyperGParams", 
                geneIds = gen.ref$EntrezGeneID,
                universeGeneIds = immunome$EntrezGeneID,
                annotation = "org.Hs.eg.db", 
                ontology = "BP", 
                pvalueCutoff = 0.01,
                conditional = FALSE, 
                testDirection = "over")

#Imprimir términos GO. Los primeros (head) serán los más representados
GOresult <- GOstats::hyperGTest(paramsGO)
head(summary(GOresult))

#Biocmanager::install("KEGG")
BiocManager::install("KEGG.db")
library(KEGG.db)
paramsKEGG <- new("KEGGHyperGParams",
                  geneIds = gen.ref$EntrezGeneID,
                  universeGeneIds = immunome$EntrezGeneID,
                  annotation = "org.Hs.eg.db", 
                  pvalueCutoff = 0.01, 
                  testDirection = "over")

over.kegg <- GOstats::hyperGTest(paramsKEGG)
summary(over.kegg)

#Pruebas estadísticas

BiocManager::install("topGO")
library(topGO)

contador <- 1
gene2GO <- c()
genevec <- c()
immunome<-inmunome

for(geneid in immunome$EntrezGeneID){
  genevec[contador] <- geneid
  gene2GO[contador] <- list(goannot2[goannot2$ENTREZID==geneid, "GO"])
  contador <- contador + 1
}

names(gene2GO) <- genevec
genelist <- rep(0,length(immunome$EntrezGeneID))
genelist[immunome$EntrezGeneID %in% gen.ref$EntrezGeneID] <- 1
names(genelist) <- genevec #Tiene el codigo del gen y pone si es 1 o 0
genelist <- factor(genelist)

GOdata.MF <- new("topGOdata", 
                 ontology = "MF", 
                 description = "MF on innate immunity genes", 
                 allGenes = genelist, 
                 annot = annFUN.gene2GO,
                 gene2GO = gene2GO)



resultFisher.MF.classic <- runTest(GOdata.MF, 
                                   algorithm = "classic", 
                                   statistic = "fisher")
allRes.MF.classic <- GenTable(GOdata.MF, 
                              classicFisher = resultFisher.MF.classic,
                              topNodes = 20)
head(allRes.MF.classic)


