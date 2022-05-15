---
title: "Practica_01_calificada"
author: "Aylas Nicole, Cajahuanca Tito, Campos Briggite y Dongo Nathaly"
date: "10/05/2021"
output:
  html_document: default
  word_document: default
  pdf_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown
R version 4.0
+ Cargando librerías necesarias
```{r, message=FALSE, results='hide', warning=FALSE, echo=FALSE}
library(msa)
library(seqinr)
library(ape)
library(Biostrings)
library(ggplot2)
library(grid)
library(gridExtra)
```


### 1. Análisis de la composición de diferentes proteínas presentes en los coronavirus:          SARS-CoV-2, SARS-CoV y MERS-CoV

<div style="text-align: justify"> Debido al contexto mundial de la pandemia, se está investango mucho acerca de la familia del coronavirus (Coronaviridae). Empecemos por el coronavirus del síndrome respiratorio de Oriente Medio (MERS-CoV), es una especie de coronavirus que infecta a humanos, murciélagos y camellos. Se informó por primera vez en 2012 y es una posible causa de epidemias futuras, tiene un receptor de entrada, el cual es DPP4. Luego tenemos al coronavirus del síndrome respiratorio agudo severo (SARS-CoV), surgió en 2003 en los países del sudeste asiático e infectó a 8422 personas en 30 de países y causó 916 muertes, tenía como huésped intermedio a civetas, perros y mapaches, y cuenta como receptor de entrada al virus ACE2. Y por último tenemos al coronavirus de tipo 2 causante del síndrome respiratorio agudo severo (SARS-CoV-2), se descubrió y se aisló por primera vez en Wuhan, China. Tiene un origen zoonótico, es decir, se transmitió de un huésped animal a uno humano y usa también el receptor ACE2 para ingresar a nuestras células. Esta última tuvo una expansión mundial provocó la pandemia de COVID-19 hasta la actualidad. Se mostrará a continuación un estudio acerca de la presencia nucleótidos en las diferentes proteínas de los virus mencionados anteriormente.</div>

<div style="text-align: justify"> A continuación se muestra con detalle las diferentes proteínas para cada tipo de virus. </div>

<div style="text-align: center">
SARS-CoV  | MERS-CoV | SARS-CoV-2  
----------|----------| ----------
ORF1ab    |  ORF1ab  | ORF1ab 
Spike     |  Spike   | Spike 
3a        |  3       | 3 
3b        |  4a      | E
E         |  5       | M
M         |  E       | 6 
6         |  M       | 7a
7a        |  N       | 8a 
7b        |  -       | N 
8a        |  -       | -
8b        |  -       | -
N         |  -       | -  
</div>

<div style="text-align: justify">El análisis se llevó a cabo construyendo una data frame a partir del genoma que codificaba para cada proteína presente en los tres virus. Se declararon las secuencias y se utilizó el paquete "ape" para leer las secuencias:  </div>

```{r, results='hide'}
SARS_COV_2_id<-c("NC_045512.2 (266..21555)", "NC_045512.2 (21563..25384)", "NC_045512.2 (25393..26220)", "NC_045512.2 (26245..26472)", "NC_045512.2 (26523..27191)", "NC_045512.2 (27202..27387)", "NC_045512.2 (27394..27759)",
                 "NC_045512.2 (27894..28259)",  "NC_045512.2 (28274..29533)")

SARS_COV_id<-c("NC_004718.3 (265..21485)", "NC_004718.3 (21492..25259)", "NC_004718.3 (25268..26092)", "NC_004718.3 (25689..26153)", "NC_004718.3 (26117..26347)", "NC_004718.3 (26398..27063)", "NC_004718.3 (26913..27265)" , "NC_004718.3 (27273..27641)", "NC_004718.3 (27638..27772)", "NC_004718.3 (27779..27898)" , "NC_004718.3 (27864..28118)" , "NC_004718.3 (28120..29388)")

MERS_COV_id<-c("NC_019843.3 (279..21514)", "NC_019843.3 (21456..25517)", "NC_019843.3 (25532..25843)","NC_019843.3 (25852..26181)", "NC_019843.3 (26093..26833)", "NC_019843.3 (26840..27514)", "NC_019843.3 (27590..27838)", "NC_019843.3 (27853..28512)", "NC_019843.3 (28566..29807)")

mSARS_COV_2<-ape::read.GenBank(SARS_COV_2_id, as.character = T)
#mSARS_COV_2 descomentar para mostrar variables

mSARS_COV<-ape::read.GenBank(SARS_COV_id, as.character = T)
#mSARS_COV descomentar para mostrar variables

mMERS_COV<-ape::read.GenBank(MERS_COV_id, as.character = T)
#mMERS_COV descomentar para mostrar variables
```

```{r, echo=FALSE}
mSARS_COV_2$ORF1ab<-mSARS_COV_2$ORF1ab[266:21555]
mSARS_COV_2$Spike<-mSARS_COV_2$Spike[21563:25384]
mSARS_COV_2$`3`<-mSARS_COV_2$`3`[25393:26220]
mSARS_COV_2$E<-mSARS_COV_2$E[26245:26472]
mSARS_COV_2$M<-mSARS_COV_2$M[26523:27191]
mSARS_COV_2$`6`<-mSARS_COV_2$`6`[27202:27387]
mSARS_COV_2$`7a`<-mSARS_COV_2$`7a`[27394:27759]
mSARS_COV_2$`8a`<-mSARS_COV_2$`8a`[27894:28259]
mSARS_COV_2$N<-mSARS_COV_2$N[28274:29533]

mSARS_COV$ORF1ab<-mSARS_COV$ORF1ab [265:21485]
mSARS_COV$Spike<-mSARS_COV$Spike [21492:25259]
mSARS_COV$`3a`<-mSARS_COV$`3a` [25268:26092]
mSARS_COV$`3b`<-mSARS_COV$`3b` [25689:26153]
mSARS_COV$E<-mSARS_COV$E[26117:26347]
mSARS_COV$M<-mSARS_COV$M [26398:27063]
mSARS_COV$`6`<-mSARS_COV$`6`[26913:27265]
mSARS_COV$`7a`<-mSARS_COV$`7a` [27273:27641]
mSARS_COV$`7b`<-mSARS_COV$`7b` [27638:27772]
mSARS_COV$`8a`<-mSARS_COV$`8a` [27779:27898]
mSARS_COV$`8b`<-mSARS_COV$`8b`[27864:28118]
mSARS_COV$N<-mSARS_COV$N[28120:29388]

mMERS_COV$ORF1ab<-mMERS_COV$ORF1ab[279:21514]
mMERS_COV$Spike<-mMERS_COV$Spike[21456:25517]
mMERS_COV$`3`<-mMERS_COV$`3`[25532:25843]
mMERS_COV$`4a`<-mMERS_COV$`4a`[25852:26181]
mMERS_COV$`4b`<-mMERS_COV$`4b`[26093:26833]
mMERS_COV$`5`<-mMERS_COV$`5`[26840:27514]
mMERS_COV$E<-mMERS_COV$E[27590:27838]
mMERS_COV$M<-mMERS_COV$M[27853:28512]
mMERS_COV$N<-mMERS_COV$N[28566:29807]

names(mSARS_COV_2)<-c("ORF1ab","Spike","3","E","M","6","7a","8a", "N")
names(mSARS_COV)<-c("ORF1ab","Spike","3a","3b","E","M","6","7a","7b","8a","8b","N")
names(mMERS_COV)<- c("ORF1ab", "Spike", "3", "4a", "4b", "5", "E", "M", "N")
```

<div style="text-align: justify">El primer paso realizado para obtener el data frame fue obtener las frecuencias absolutas de cada nucleótido por proteína con respecto al virus al que pertenece.</div>
```{r}
aSARS_COV_2_ORF1ab<-sort(table(mSARS_COV_2$ORF1ab))
```

```{r, echo=FALSE}
#Frecuencias absolutas 
aSARS_COV_2_ORF1ab<-sort(table(mSARS_COV_2$ORF1ab))
aSARS_COV_2_Spike<-sort(table(mSARS_COV_2$Spike))
aSARS_COV_2_3<-sort(table(mSARS_COV_2$`3`))
aSARS_COV_2_E<-sort(table(mSARS_COV_2$E))
aSARS_COV_2_M<-sort(table(mSARS_COV_2$M))
aSARS_COV_2_6<-sort(table(mSARS_COV_2$`6`))
aSARS_COV_2_7a<-sort(table(mSARS_COV_2$`7a`))
aSARS_COV_2_8a<-sort(table(mSARS_COV_2$`8a`))
aSARS_COV_2_N<-sort(table(mSARS_COV_2$N))

aSARS_COV_ORF1ab<-sort(table(mSARS_COV$ORF1ab))
aSARS_COV_Spike<-sort(table(mSARS_COV$Spike))
aSARS_COV_3a<-sort(table(mSARS_COV$`3a`))
aSARS_COV_3b<-sort(table(mSARS_COV$`3b`))
aSARS_COV_E<-sort(table(mSARS_COV$E))
aSARS_COV_M<-sort(table(mSARS_COV$M))
aSARS_COV_6<-sort(table(mSARS_COV$`6`))
aSARS_COV_7a<-sort(table(mSARS_COV$`7a`))
aSARS_COV_7b<-sort(table(mSARS_COV$`7b`))
aSARS_COV_8a<-sort(table(mSARS_COV$`8a`))
aSARS_COV_8b<-sort(table(mSARS_COV$`8b`))
aSARS_COV_N<-sort(table(mSARS_COV$N))

aMERS_COV_ORF1ab<-sort(table(mMERS_COV$ORF1ab))
aMERS_COV_Spike<-sort(table(mMERS_COV$Spike))
aMERS_COV_3<-sort(table(mMERS_COV$`3`))
aMERS_COV_4a<-sort(table(mMERS_COV$`4a`))
aMERS_COV_4b<-sort(table(mMERS_COV$`4b`))
aMERS_COV_5<-sort(table(mMERS_COV$`5`))
aMERS_COV_E<-sort(table(mMERS_COV$E))
aMERS_COV_M<-sort(table(mMERS_COV$M))
aMERS_COV_N<-sort(table(mMERS_COV$N))
```

<div style="text-align: justify">De la misma manera se procedió a realizar el mismo procedimiento, pero en este caso para hallar las frecuencias relativas.</div>
```{r}
mSARS_COV_2_ORF1ab<-sort(table(mSARS_COV_2$ORF1ab))/length(mSARS_COV_2$ORF1ab)
```

```{r, echo=FALSE}
#Frecuencias relativas
mSARS_COV_2_ORF1ab<-sort(table(mSARS_COV_2$ORF1ab))/length(mSARS_COV_2$ORF1ab)
mSARS_COV_2_Spike<-sort(table(mSARS_COV_2$Spike))/length(mSARS_COV_2$Spike)
mSARS_COV_2_3<-sort(table(mSARS_COV_2$`3`))/length(mSARS_COV_2$`3`)
mSARS_COV_2_E<-sort(table(mSARS_COV_2$E))/length(mSARS_COV_2$E)
mSARS_COV_2_M<-sort(table(mSARS_COV_2$M))/length(mSARS_COV_2$M)
mSARS_COV_2_6<-sort(table(mSARS_COV_2$`6`))/length(mSARS_COV_2$`6`)
mSARS_COV_2_7a<-sort(table(mSARS_COV_2$`7a`))/length(mSARS_COV_2$`7a`)
mSARS_COV_2_8a<-sort(table(mSARS_COV_2$`8a`))/length(mSARS_COV_2$`8a`)
mSARS_COV_2_N<-sort(table(mSARS_COV_2$N))/length(mSARS_COV_2$N)

mSARS_COV_ORF1ab<-sort(table(mSARS_COV$ORF1ab))/length(mSARS_COV$ORF1ab)
mSARS_COV_Spike<-sort(table(mSARS_COV$Spike))/length(mSARS_COV$Spike)
mSARS_COV_3a<-sort(table(mSARS_COV$`3a`))/length(mSARS_COV$`3a`)
mSARS_COV_3b<-sort(table(mSARS_COV$`3b`))/length(mSARS_COV$`3b`)
mSARS_COV_E<-sort(table(mSARS_COV$E))/length(mSARS_COV$E)
mSARS_COV_M<-sort(table(mSARS_COV$M))/length(mSARS_COV$M)
mSARS_COV_6<-sort(table(mSARS_COV$`6`))/length(mSARS_COV$`6`)
mSARS_COV_7a<-sort(table(mSARS_COV$`7a`))/length(mSARS_COV$`7a`)
mSARS_COV_7b<-sort(table(mSARS_COV$`7b`))/length(mSARS_COV$`7b`)
mSARS_COV_8a<-sort(table(mSARS_COV$`8a`))/length(mSARS_COV$`8a`)
mSARS_COV_8b<-sort(table(mSARS_COV$`8b`))/length(mSARS_COV$`8b`)
mSARS_COV_N<-sort(table(mSARS_COV$N))/length(mSARS_COV$N)

mMERS_COV_ORF1ab<-sort(table(mMERS_COV$ORF1ab))/length(mMERS_COV$ORF1ab)
mMERS_COV_Spike<-sort(table(mMERS_COV$Spike))/length(mMERS_COV$Spike)
mMERS_COV_3<-sort(table(mMERS_COV$`3`))/length(mMERS_COV$`3`)
mMERS_COV_4a<-sort(table(mMERS_COV$`4a`))/length(mMERS_COV$`4a`)
mMERS_COV_4b<-sort(table(mMERS_COV$`4b`))/length(mMERS_COV$`4b`)
mMERS_COV_5<-sort(table(mMERS_COV$`5`))/length(mMERS_COV$`5`)
mMERS_COV_E<-sort(table(mMERS_COV$E))/length(mMERS_COV$E)
mMERS_COV_M<-sort(table(mMERS_COV$M))/length(mMERS_COV$M)
mMERS_COV_N<-sort(table(mMERS_COV$N))/length(mMERS_COV$N)
```

<div style="text-align: justify">El segundo paso fue agrupar los parámetros: nucleótido, frecuencia relativa y frecuencia absoluta dentro de un data frame por cada proteína con respecto a su virus.</div>
```{r}
mSARS_COV_2_ORF1ab<-data.frame(n=names(mSARS_COV_2_ORF1ab),valor=as.vector(mSARS_COV_2_ORF1ab),valorR=as.vector(aSARS_COV_2_ORF1ab))
```

```{r, echo=FALSE}
#Creamos una tabla de doble entrada
mSARS_COV_2_Spike<-data.frame(n=names(mSARS_COV_2_Spike),valor=as.vector(mSARS_COV_2_Spike),valorR=as.vector(aSARS_COV_2_Spike))
mSARS_COV_2_3<-data.frame(n=names(mSARS_COV_2_3),valor=as.vector(mSARS_COV_2_3),valorR=as.vector(aSARS_COV_2_3))
mSARS_COV_2_E<-data.frame(n=names(mSARS_COV_2_E),valor=as.vector(mSARS_COV_2_E),valorR=as.vector(aSARS_COV_2_E))
mSARS_COV_2_M<-data.frame(n=names(mSARS_COV_2_M),valor=as.vector(mSARS_COV_2_M),valorR=as.vector(aSARS_COV_2_M))
mSARS_COV_2_6<-data.frame(n=names(mSARS_COV_2_6),valor=as.vector(mSARS_COV_2_6),valorR=as.vector(aSARS_COV_2_6))
mSARS_COV_2_7a<-data.frame(n=names(mSARS_COV_2_7a),valor=as.vector(mSARS_COV_2_7a),valorR=as.vector(aSARS_COV_2_7a))
mSARS_COV_2_8a<-data.frame(n=names(mSARS_COV_2_8a),valor=as.vector(mSARS_COV_2_8a),valorR=as.vector(aSARS_COV_2_8a))
mSARS_COV_2_N<-data.frame(n=names(mSARS_COV_2_N),valor=as.vector(mSARS_COV_2_N),valorR=as.vector(aSARS_COV_2_N))

mSARS_COV_ORF1ab<-data.frame(n=names(mSARS_COV_ORF1ab),valor=as.vector(mSARS_COV_ORF1ab),valorR=as.vector(aSARS_COV_ORF1ab))
mSARS_COV_Spike<-data.frame(n=names(mSARS_COV_Spike),valor=as.vector(mSARS_COV_Spike),valorR=as.vector(aSARS_COV_Spike))
mSARS_COV_3a<-data.frame(n=names(mSARS_COV_3a),valor=as.vector(mSARS_COV_3a),valorR=as.vector(aSARS_COV_3a))
mSARS_COV_3b<-data.frame(n=names(mSARS_COV_3b),valor=as.vector(mSARS_COV_3b),valorR=as.vector(aSARS_COV_3b))
mSARS_COV_E<-data.frame(n=names(mSARS_COV_E),valor=as.vector(mSARS_COV_E),valorR=as.vector(aSARS_COV_E))
mSARS_COV_M<-data.frame(n=names(mSARS_COV_M),valor=as.vector(mSARS_COV_M),valorR=as.vector(aSARS_COV_M))
mSARS_COV_6<-data.frame(n=names(mSARS_COV_6),valor=as.vector(mSARS_COV_6),valorR=as.vector(aSARS_COV_6))
mSARS_COV_7a<-data.frame(n=names(mSARS_COV_7a),valor=as.vector(mSARS_COV_7a),valorR=as.vector(aSARS_COV_7a))
mSARS_COV_7b<-data.frame(n=names(mSARS_COV_7b),valor=as.vector(mSARS_COV_7b),valorR=as.vector(aSARS_COV_7b))
mSARS_COV_8a<-data.frame(n=names(mSARS_COV_8a),valor=as.vector(mSARS_COV_8a),valorR=as.vector(aSARS_COV_8a))
mSARS_COV_8b<-data.frame(n=names(mSARS_COV_8b),valor=as.vector(mSARS_COV_8b),valorR=as.vector(aSARS_COV_8b))
mSARS_COV_N<-data.frame(n=names(mSARS_COV_N),valor=as.vector(mSARS_COV_N),valorR=as.vector(aSARS_COV_N))

mMERS_COV_ORF1ab<-data.frame(n=names(mMERS_COV_ORF1ab),valor=as.vector(mMERS_COV_ORF1ab),valorR=as.vector(aMERS_COV_ORF1ab))
mMERS_COV_Spike<-data.frame(n=names(mMERS_COV_Spike),valor=as.vector(mMERS_COV_Spike),valorR=as.vector(aMERS_COV_Spike))
mMERS_COV_3<-data.frame(n=names(mMERS_COV_3),valor=as.vector(mMERS_COV_3),valorR=as.vector(aMERS_COV_3))
mMERS_COV_4a<-data.frame(n=names(mMERS_COV_4a),valor=as.vector(mMERS_COV_4a),valorR=as.vector(aMERS_COV_4a))
mMERS_COV_4b<-data.frame(n=names(mMERS_COV_4b),valor=as.vector(mMERS_COV_4b),valorR=as.vector(aMERS_COV_4b))
mMERS_COV_5<-data.frame(n=names(mMERS_COV_5),valor=as.vector(mMERS_COV_5),valorR=as.vector(aMERS_COV_5))
mMERS_COV_E<-data.frame(n=names(mMERS_COV_E),valor=as.vector(mMERS_COV_E),valorR=as.vector(aMERS_COV_E))
mMERS_COV_M<-data.frame(n=names(mMERS_COV_M),valor=as.vector(mMERS_COV_M),valorR=as.vector(aMERS_COV_M))
mMERS_COV_N<-data.frame(n=names(mMERS_COV_N),valor=as.vector(mMERS_COV_N),valorR=as.vector(aMERS_COV_N))
```

<div style="text-align: justify">Así mismo, se concatenaron los nucleótidos de cada proteína, los valores de la frecuencia relativa y de la frecuencia absoluta.</div>

```{r}
#Nucleotidos
n_mSARS_COV_2=c(as.vector(mSARS_COV_2_ORF1ab$n),as.vector(mSARS_COV_2_Spike$n),as.vector(mSARS_COV_2_3$n),as.vector(mSARS_COV_2_E$n),as.vector(mSARS_COV_2_M$n),as.vector(mSARS_COV_2_6$n),as.vector(mSARS_COV_2_7a$n),as.vector(mSARS_COV_2_8a$n),as.vector(mSARS_COV_2_N$n))
n_mSARS_COV=c(as.vector(mSARS_COV_ORF1ab$n),as.vector(mSARS_COV_Spike$n),as.vector(mSARS_COV_3a$n),as.vector(mSARS_COV_3b$n),as.vector(mSARS_COV_E$n),as.vector(mSARS_COV_M$n),as.vector(mSARS_COV_6$n),as.vector(mSARS_COV_7a$n),as.vector(mSARS_COV_7b$n),as.vector(mSARS_COV_8a$n),as.vector(mSARS_COV_8b$n),as.vector(mSARS_COV_N$n))
n_mMERS_COV=c(as.vector(mMERS_COV_ORF1ab$n),as.vector(mMERS_COV_Spike$n),as.vector(mMERS_COV_3$n),as.vector(mMERS_COV_4a$n),as.vector(mMERS_COV_4b$n),as.vector(mMERS_COV_5$n),as.vector(mMERS_COV_E$n),as.vector(mMERS_COV_M$n),as.vector(mMERS_COV_N$n))

#Valores frecuencia absoluta
valor_mSARS_COV_2=c(mSARS_COV_2_ORF1ab$valor,mSARS_COV_2_Spike$valor,mSARS_COV_2_3$valor,mSARS_COV_2_E$valor,mSARS_COV_2_M$valor,mSARS_COV_2_6$valor,mSARS_COV_2_7a$valor,mSARS_COV_2_8a$valor,mSARS_COV_2_N$valor)
valor_mSARS_COV=c(mSARS_COV_ORF1ab$valor,mSARS_COV_Spike$valor,mSARS_COV_3a$valor,mSARS_COV_3b$valor,mSARS_COV_E$valor,mSARS_COV_M$valor,mSARS_COV_6$valor,mSARS_COV_7a$valor,mSARS_COV_7b$valor,mSARS_COV_8a$valor,mSARS_COV_8b$valor,mSARS_COV_N$valor)
valor_mMERS_COV=c(mMERS_COV_ORF1ab$valor,mMERS_COV_Spike$valor,mMERS_COV_3$valor,mMERS_COV_4a$valor,mMERS_COV_4b$valor,mMERS_COV_5$valor,mMERS_COV_E$valor,mMERS_COV_M$valor,mMERS_COV_N$valor)

#Valores frecuencia relativa
valorR_aSARS_COV_2=c(aSARS_COV_2_ORF1ab,aSARS_COV_2_Spike,aSARS_COV_2_3,aSARS_COV_2_E,aSARS_COV_2_M,aSARS_COV_2_6,aSARS_COV_2_7a,aSARS_COV_2_8a,aSARS_COV_2_N)
valorR_aSARS_COV=c(aSARS_COV_ORF1ab,aSARS_COV_Spike,aSARS_COV_3a,aSARS_COV_3b,aSARS_COV_E,aSARS_COV_M,aSARS_COV_6,aSARS_COV_7a,aSARS_COV_7b,aSARS_COV_8a,aSARS_COV_8b,aSARS_COV_N)
valorR_mMERS_COV=c(aMERS_COV_ORF1ab,aMERS_COV_Spike,aMERS_COV_3,aMERS_COV_4a,aMERS_COV_4b,aMERS_COV_5,aMERS_COV_E,aMERS_COV_M,aMERS_COV_N)
```

```{r, echo=FALSE}
# Virus
virus_mSARS_COV_2=c(rep("SARS_COV_2", each=length(names(mSARS_COV_2))*4))
virus_mSARS_COV=c(rep("SARS_COV", each=length(names(mSARS_COV))*4))
virus_mMERS_COV=c(rep("MERS_COV", each=length(names(mMERS_COV))*4))
```

<div style="text-align: justify">De esta manera se creo una base de datos donde se pueden observar los tres virus y sus respectivas proteínas, nucleótidos, frecuencias absolutas y relativas.</div>

```{r}
dt <- data.frame(virus=c(virus_mSARS_COV_2,virus_mSARS_COV,virus_mMERS_COV),
                 proteina= rep(c(names(mSARS_COV_2),names(mSARS_COV),names(mMERS_COV)),each=4),
                 n=c(n_mSARS_COV_2,n_mSARS_COV,n_mMERS_COV),
                 FA= c(valorR_aSARS_COV_2,valorR_aSARS_COV,valorR_mMERS_COV),
                 FR=c(valor_mSARS_COV_2,valor_mSARS_COV,valor_mMERS_COV))
head(dt)
```

### Análisis Estadístico Descriptivo 
+ Cargando paquetes necesarios para el análisis:

<div style="text-align: justify">A continuación se graficaron los nucleótidos de cada proteína por cada virus utilizando el paquete ggplot2.</div>

```{r}
df1_filtered<-dt[dt$virus=="MERS_COV",]
p <- ggplot(df1_filtered, aes(x=proteina, y=FR, fill=n)) + 
  geom_bar(stat="identity", position=position_dodge()) 

df2_filtered<-dt[dt$virus=="SARS_COV",]
q <- ggplot(df2_filtered, aes(x=proteina, y=FR, fill=n)) + 
  geom_bar(stat="identity", position=position_dodge()) 

df3_filtered<-dt[dt$virus=="SARS_COV_2",]
r <- ggplot(df3_filtered, aes(x=proteina, y=FR, fill=n)) + 
  geom_bar(stat="identity", position=position_dodge()) 

grid.arrange(
  
  p + scale_fill_brewer(palette="Paired") + theme(axis.text.x = element_text(angle = 90)) + theme_minimal() + ggtitle("Frecuencia relativa de nucleotidos en MERS-CoV") +
    xlab("") + ylab(""), 
  
  q + scale_fill_brewer(palette="Paired") + theme(axis.text.x = element_text(angle = 90)) + theme_minimal() + ggtitle("Frecuencia relativa de nucleótidos en SARS-CoV") +
    xlab("") + ylab("Frecuencia relativa"),
  
  r + scale_fill_brewer(palette="Paired") + theme(axis.text.x = element_text(angle = 90)) + theme_minimal() + ggtitle("Frecuencia relativa de nucleótidos en SARS-CoV-2") +
    xlab("Proteina") + ylab("")
  
  , nrow=3)

```

### Análisis Estadístico Inferencial
+ Mediante Anova de dos vías y test Tukey
```{r}
fit<-aov(dt$FR~dt$n*dt$proteina)
tukeyobj<-TukeyHSD(fit)
View(tukeyobj[[3]])
#TukeyHSD(fit)
#summary(fit)
```

<div style="text-align: justify">Una vez realizado se eligen los valores significativos para poder analizarlo. Estos valores están asociado bajo la misma proteína y núcleótidos complementarios, pues bajo esta premisa se plantea la hipótesis nula:  </div>

<div style="text-align: justify"> H0: Los nucleotidos de una misma proteína es similar su base de nucleotido complementario </div>

<div style="text-align: justify"> A continuación se muestran los datos elegidos del test Tukey </div>

<div style="text-align: center">

N°  | Comparación         | p-Value 
----|---------------------| -------
1   | t:ORF1ab-a:ORF1ab   | 7.64E-01
2   | g:ORF1ab-c:ORF1ab   | 1.00E+00
3   |      g:N-c:N        | 6.58E-01
4   |      t:N-a:N        | 2.50E-07
5   |  t:Spike-a:Spike    | 8.33E-03
6   |  g:Spike-c:Spike    | 1.00E+00
7   |      g:M-c:M        | 1.00E+00
8   |      t:M-a:M        | 5.79E-03
9   |      g:E-c:E        | 1.00E+00
10  |      t:E-a:E        | 9.12E-12
11  |      t:3-a:3        | 1.92E-03

</div>

<div style="text-align: justify"> En el caso 1, se observa a la proteían ORF1ab con su par de base complementario, pero que poseen una media significativamente distinta, dado probablemente por ser de coronavirus diferentes, SARS-COV, SARS-Cov-2 y MERS-CoV. A diferencia de ello, para las bases complementarias G-C la media es similar, dado que el p-value es mayor a 0.05, esto puede ser porque coincide la secuencia de uno de los coronavirus </div>

<div style="text-align: justify"> Por otro lado, en el caso 3 y 4, se observa la proteína N, que a pesar de ser en ambos casos, nucleotidos de bases complementarias poseen una media significativamente distinta, posiblemente porque pertenecen a coronavirus diferentes, SARS-COV, SARS-Cov-2 y MERS-CoV </div>

<div style="text-align: justify"> Bajo este mismo analisis, en el caso 5 se observa la proteína Spike para el par de base T-A, una media significativamente distinta al igual que el caso 1. Sin embargo, para el caso 6 que compara la misma proteína pero el par de base G-C, no se puede rechazar la hiótesis, por lo tanto se asume que existe una similitud </div>

<div style="text-align: justify"> Para la proteína M, sucede lo mismo pues el caso 7 muestra el par de base G-C, en el cual se observa que existe una similitud. MInetras que en el caso 8 que tiene el par de base T-A,  presenta una media significativamente distinta, posiblemente porque pertenecen a coronavirus diferentes, SARS-Cov-2 y MERS. </div>

<div style="text-align: justify"> Ahora bien para el caso de la proteína E, el comportaiento es igual, para el par de base G-C, existe una similitud. Mientras que para el par de base T-A,  presenta una media significativamente distinta, posiblemente porque pertenecen a coronavirus diferentes, SARS-Cov-2 y MERS. </div>

<div style="text-align: justify"> En el caso 11, se observa que a pesar de ser nucleotidos de bases complementarias poseen una media significativamente distinta, posiblemente porque pertenecen a coronavirus diferentes, SARS-Cov-2 y MERS </div>

<div style="text-align: justify"> Finalmente del analisis se puede concluir que los pares de bases GC para las proteínas si tienen similitud; es decir que si cumplen bajo la especulación que al ser complementarias estás van a ser similares. No obstante, para el par de base TA difiere, pues aquí a pesar de ser complementarios tienen una media significativamente distinta, porque son diferentes virs. Es importante resaltar que no se cumple para la proteína N, pues aquí en ambos pares de bases son distintos. Estos cambios, probablemente son porque los virus difieren en las secuencias para sus proteínas. </div>

### 2. Analisis de SARS-CoV-2 en distintas especies

<div style="text-align: justify"> Se sospecha que el coronavirus 2 del síndrome respiratorio agudo severo (SARS-CoV-2, anteriormente 2019-nCoV) se originó en 2019 en China a partir de un murciélago infectado por coronavirus del género Rhinolophus. Tras la aparición inicial, posiblemente facilitada por un huésped puente de mamíferos, el SARS-CoV-2 se transmite actualmente en todo el mundo a través de una transmisión eficiente de persona a persona. Los resultados obtenidos de estudios experimentales indican que especies animales como gatos, hurones, perros mapache, macacos cynomolgus, macacos rhesus, venado cola blanca, conejos, murciélagos frugívoros egipcios y hámsteres sirios son susceptibles a la infección por SARS-CoV-2, y que La transmisión de gato a gato y de hurón a hurón puede tener lugar a través del contacto y el aire. Sin embargo, se han informado infecciones naturales de SARS-CoV-2 solo en perros y gatos domésticos, tigres, leones, leopardos de las nieves, pumas, visones y hurones de granja. Como se mostrará a continuación se realizó un estudios sobre los aminoácidos de la proteínas Spike en distintos organismos. </div>

<div style="text-align: center">
Especie                  | Nombre Común  
------------------------ | -------------
Hommo sapiens            | Humano
Chlorocebus sabaeus      | Mono 
Canis lupus familiaris   | Perro
Felis catus              | Gato
Panthera tigris          | Tigre
Panthera leo             | León
Neovison vison           | Vison Americano
Mustela lutreola         | Vison Europeo
Mustela putorius furo    | Hurón
Mesocricetus auratus     | Hámster Dorado
</div>

<div style="text-align: justify">El análisis se llevo a cabo construyendo una data frame a partir de los número de acceso del genoma que codificaban para una proteina en cada organismo. Se declararon las secuencias y se utilizó el paquete "ape" para leer las secuencias:  </div>
```{r, results='hide'}

secuencias_org_id<-c("NC_045512.2","MW477798", "MT270814" , "MW064259", 
                     "MT704314.1" , "MT704310" , "MW626378.1" , "MW562252" , 
                     "MZ099821.1" , "MT835139")
secuencias_org<-ape::read.GenBank(secuencias_org_id, as.character = T)

```

```{r, echo=FALSE}
secuencias_org$NC_045512.2<-secuencias_org$NC_045512.2[21563:25384]
secuencias_org$MW477798<-secuencias_org$MW477798[21534:25355]
secuencias_org$MT270814<-secuencias_org$MT270814[21524:25345]
secuencias_org$MW064259<-secuencias_org$MW064259[21298:25119]
secuencias_org$MT704314.1<-secuencias_org$MT704314.1[21563:25384] 
secuencias_org$MT704310<-secuencias_org$MT704310[21563:25384] 
secuencias_org$MW626378.1<-secuencias_org$MW626378.1[21525:25346]
secuencias_org$MW562252<-secuencias_org$MW562252[21575:25393]
secuencias_org$MZ099821.1<-secuencias_org$MZ099821.1[21538:25359]
secuencias_org$MT835139<-secuencias_org$MT835139[21509:25330]


n_especies<-c("Hommo sapiens","Chlorocebus sabaeus","Canis lupus familiaris",
                         "Felis catus","Panthera tigris","Panthera leo","Neovison vison",
                         "Mustela lutreola","Mustela putorius furo","Mesocricetus auratus")

names(secuencias_org)<-c("Humano","Mono verde","Perro","Gato","Tigre","Leon",
                         "Vision A.","Vision E.", "Hurón","Hamster")
```

<div style="text-align: justify">Así mismo, se crearon dos funciones que reciben como parámetros de entrada los nombres de las especies, los nombres comúnes, y las secuencias que se desean introducir. Posteriormente devuelven una data frame que incluyen también el nombre común, la especie, y la frecuencia relativa y absoluta que en el caso de la función para ADN será de los nucleótidos, y de la función para proteínas de aminoácidos.

Como se puede ver en la frecuencia de proteínas, se obtuvo a partir de la traducción de la cadena de ADN, por lo que se utilizaron la función DNAString() y translate() del paquete BioStrings.</div>

```{r, message=FALSE}

frecuencia_proteinas <- function(secuencias, especies, nombres ) {
  nucleotids_all<-c()
  pFA<-c()
  pFR<-c()
  nombres_all<-c()
  especies_all<-c()
  AA_all<-c()
  adn_FA<-c()
  adn_FR<-c()
  proteinas<-list()
  
    for (i in 1:length(secuencias)){
      #Calculando las frecuencias relativas y absolutas -- ADN
      dna <- Biostrings::DNAString(paste(secuencias[[i]], collapse=""))
      proteina<-Biostrings::translate(dna)
      proteina<-as.character(proteina)
      proteinas[i]<-proteina
      proteina <- strsplit(proteina, "")[[1]]
      pFA<-c(pFA,as.vector(table(proteina)))
      pFR<-c(pFR,as.vector(table(proteina)/length(proteina)))
      # Nombres 
      nombres_all<-c(nombres_all,rep(names(secuencias)[i],length(as.vector(table(proteina)))))  
      # Especies
      especies_all<-c(especies_all,rep(especies[i],length(as.vector(table(proteina)))))
      AA_all<-c(AA_all,names(table(proteina)))
    } 
    
  data.frame(nombre_cientifico=especies_all,
               nombre_comun=nombres_all,
               AA=AA_all, 
               FA=pFA, 
               FR=pFR, secuencia_prot=proteinas)
}

frecuencia_ADN<- function(secuencias, especies, nombres ) {
  nucleotids_all<-c()
  pFA<-c()
  pFR<-c()
  nombres_all<-c()
  especies_all<-c()
  AA_all<-c()
  adn_FA<-c()
  adn_FR<-c()

    for (i in 1:length(secuencias)){
      #Calculando las frecuencias relativas y absolutas -- ADN
      adn_FA<-c(adn_FA,as.vector(table(secuencias[[i]])))
      adn_FR<-c(adn_FR,as.vector(table(secuencias[[i]])/length(secuencias[[i]])))
      # Nombres 
      nombres_all<-c(nombres_all,rep(names(secuencias)[i],length(table(secuencias[[i]]))))  
      # Especies
      especies_all<-c(especies_all,rep(especies[i],length(table(secuencias[[i]]))))
      nucleotids_all<-c(nucleotids_all,names(table(secuencias[[i]])))
    } 
    
  data.frame(nombre_cientifico=especies_all,
               nombre_comun=nombres_all,
               nucleotido=nucleotids_all, 
               FA=adn_FA, 
               FR=adn_FR)

}

```


<div style="text-align: justify">De esta manera creamos dos bases de datos para las secuencias en formato ADN y en formato proteínas (AA).</div>

```{r, echo=FALSE}
df1<-frecuencia_ADN(secuencias_org,n_especies, names(secuencias_org))

df2<-frecuencia_proteinas(secuencias_org[-2],n_especies[-2], names(secuencias_org[-2]))
#La razón por la cuál se excluye la secuencia N° 2 es porque falta algunas bases, por lo que no pueden ser traducidad por el paquete
```

+ Base de datos para frecuencias de nucleótidos (ADN).
```{r, echo=FALSE}
head(df1)
```

+ Base de datos para frecuencias de aminoácidos (proteínas).
```{r, echo=FALSE}
head(df2[,1:5])
```

### *Análisis Estadístico Descriptivo*
  
+ Cargando paquetes necesarios para el análisis:

```{r, message=FALSE}
p <- ggplot(df1, aes(x=nombre_comun, y=FR, fill=nucleotido)) + 
   geom_bar(stat="identity", position=position_dodge()) 

df1_filtered<-df1[df1$nombre_comun!="Mono verde",]
q <- ggplot(df1_filtered, aes(x=nombre_comun, y=FR, fill=nucleotido)) + 
   geom_bar(stat="identity", position=position_dodge()) 

grid.arrange(
  
  p + scale_fill_brewer(palette="Paired") + theme(axis.text.x = element_text(angle = 90)) + theme_minimal() + ggtitle("                            Frecuencia relativa de nucleotidos") +
  xlab("Organismos") + ylab("Frecuencia relativa"), 
  
  q + scale_fill_brewer(palette="Paired") + theme(axis.text.x = element_text(angle = 90)) + theme_minimal() + ggtitle("                            Frecuencia relativa de nucleótidos filtrada") +
  xlab("Organismos") + ylab("Frecuencia relativa")
  
  , nrow=2)
```

```{r, message=FALSE}
ggplot(df2, aes(x=FR, y=AA, fill=nombre_comun))+  
  geom_boxplot(show.legend = F, aes(fill=nombre_comun))+
  geom_point(position=position_jitterdodge(jitter.width=2, dodge.width = 0), 
             pch=21, aes(fill=factor(nombre_comun)), show.legend = T)+ 
  ggtitle("                            Frecuencia relativa de Aminoácidos") +
  xlab("Frecuencia relativa") + ylab("Aminoácidos")

```

### *Análisis Estadístico Inferencial*
+ *Anova de una vía*
```{r, message=FALSE}
fit<-aov(df2$FR~df2$AA)
summary(fit)
```
```{r}
longitud_vector<-c()
for (i in secuencias_org){
  longitud_vector<-c(longitud_vector,(length(i)))
}
data.frame(nombres=names(secuencias_org),Longitud=longitud_vector)
```

*Matriz de distancia y arbol filogenetico de todo el genoma*
Para esta parte se utilizó el paquete Biostrings, específicamente la función DNAStringSet, devuelve un objeto que tiene una clase particular, DNA. Podemos crear fácilmente a DNAStringSet de nuestro vector el cual se extrajo como se muestra en el bucle for del siguiente codigo. 

También se utilizó la función msa del paquete msa, que por defecto escoge a msaClustalW, esta función llama el alineamiento multiple con el algoritmo ClustalW, tiene como entrada la secuencia (objeto definido con DNAStringSet), el tipo de data "DNA", los demás parámetros se utilizaron como default a excepcion de la matriz, que era clustalw 

Tambien se hizo un arbol filogenético obtenido a partir de la matriz de distancia (con la función dist.alignment()). Por otro lado, la función que se uso con para el arbol fue nj(). Esta función realiza la estimación del árbol de unión de vecinos de Saitou y Nei (1987). 

```{r, message=FALSE, results=FALSE}
nombres_org<-c("Hommo_sapiens","Chlorocebus_sabaeus","Canis_lupus_familiaris",
                      "Felis_catus","Panthera_tigris","Panthera_leo","Neovison_vison",
               "Mustela_lutreola", "Mustela_putorius_furo", "Mesocricetus_auratus") #redefinir nombres ya que 
lista_secuencias<-c()
for (i in secuencias_org){
  lista_secuencias<-c(lista_secuencias,(paste(i, collapse="")))
}
lista_secuencias<-Biostrings::DNAStringSet(lista_secuencias)
names(lista_secuencias)<-nombres_org

matrix_Spike<-msa::msa(lista_secuencias,  type="DNA", substitutionMatrix="clustalw")
# msaConvert() that allows for converting multiple sequence alignment objects to other types/classes.
Aln_Spike <- msaConvert(matrix_Spike, type="seqinr::alignment")

library(seqinr) #Carga la libreria seqinr
d_Spike <- dist.alignment(Aln_Spike, "identity")
matrix_Spike<-as.matrix(d_Spike)[2:10,1, drop=FALSE]
colnames(matrix_Spike)<-c("Distancia Homo sapiens (humano)")
#Construcción del arbol filogenetico
Spike_tree <- ape::nj(d_Spike)
plot(Spike_tree, main="Arbol filogenetico de proteinas Spike en distintos organismos")

```

La matriz de distancia es: 
```{r, echo=FALSE}
print(matrix_Spike)
```

Se realiza el mismo procedimiento con las scuencias de todo el genoma obtenido con la función  de ape::read.GenBank(). 

```{r}
secuencias_org_all_genome<-ape::read.GenBank(secuencias_org_id, as.character = T)
                                
```

Para correr obtener el árbol se debe descomentar el siguiente codigo, por ahora no se correrá ya que toma mucho tiempo.
```{r, echo=TRUE}
# names(secuencias_org_all_genome)<- nombres_org 
# lista_secuencias_all<-c()
# for (i in secuencias_org_all_genome){
#   lista_secuencias_all<-c(lista_secuencias_all,(paste(i, collapse="")))
# }
# lista_secuencias_all<-Biostrings::DNAStringSet(lista_secuencias_all)
# names(lista_secuencias_all)<-nombres_org
# 
# #Produciendo alineamiento
# matrix_Aln_all<-msa::msa(lista_secuencias_all,  type="DNA", substitutionMatrix="clustalw")
# # msaConvert() that allows for converting multiple sequence alignment objects to other types/classes.
# Aln_all <- msaConvert(matrix_Aln_all, type="seqinr::alignment")
# 
# library(seqinr) #Carga la libreria seqinr
# d_all <- dist.alignment(Aln_all, "identity")
# matrix_all_dist<-as.matrix(d_all)[2:10,1, drop=FALSE]
# colnames(matrix_all_dist)<-c("Distancia Homo sapiens (humano)")
# 
# #Construcción del arbol filogenetico
# All_genome_Tree <- ape::nj(d_all)
# plot(All_genome_Tree, main="Arbol filogenetico del genoma entero en distintos organismos")

```

La matriz de distancia es: 
```{r, echo=TRUE}
# print(matrix_all_dist)
```

Graficando ambos árboles filogenéticos.Vemos que la distribución es la misma que la anterior, significa que los árboles que se generar a partir de la proteína Spike son los mismos que analizando el genoma completo. Sería util realizar el mismo análisis para genomas recolectados de distintos países. 

```{r, message=FALSE}
# par(mfrow=c(1,2))
# plot(Spike_tree, main="Arbol filogenetico de todo el  \n genoma en distintos organismos")
# plot(All_genome_Tree, main="Arbol filogenetico de proteinas \n Spike en distintos organismos")
```
+ *Análisis de la composición de aminoacidos y el contenido GC*
A continuacion se declaró una funcion que tiene recibe como entrada un objeto de tipo lista, que en este caso es secuencias_org_all_genome, el objeto que contiene los genomas completos de cada organismo leídos anteriormente por la función ape::read.GenBank. La función también recibe la posición, dado que sabemos que tenemos 10 organismos. Dentro de la función se define un tamaño de ventana, entonces se genera un vector que contiene el porcentage de GC por cada segmento que pertenece al paso, cuyo ancho es igual al del tamaño de la ventana. 

```{r}
#Analisis sobre el contedido GC funcion
GC_porcentaje <- function (secuencias_org_input,pos){
  a<-list()
  contador_2<-1
      for (tamano.ventana in c(200,500,1000)){
        porcentaje_GC<-vector()
        contador<-1
        for (i in seq(from=1,to=length(secuencias_org_input[[pos]])-tamano.ventana,by=tamano.ventana/10)){
          j<-i+tamano.ventana
          porcentaje_GC[contador]<- seqinr::GC(secuencias_org_input[[pos]][i:j])
          contador<-contador+1
        }
        porcentaje_GC <- porcentaje_GC[!is.na(porcentaje_GC)]
        a[[contador_2]]<-porcentaje_GC
        contador_2<-contador_2+1
      }
  
  return(a)
}

#Aplicando la funcion que extrae el contenido GC
Humano_GC<-GC_porcentaje(secuencias_org_all_genome,1)
```

Y así lo definimos para el resto de secuencias, por otro lado, graficamos los resultados para poder ver los picos de GC para ventanas de 200, 500 y 1000 nucleótidos, como se muestra a continuación. 

```{r, echo=FALSE}
Mono_GC<-GC_porcentaje(secuencias_org_all_genome,2)
Perro_GC<-GC_porcentaje(secuencias_org_all_genome,3)
Gato_GC<-GC_porcentaje(secuencias_org_all_genome,4)
Tigre_GC<-GC_porcentaje(secuencias_org_all_genome,5)
Leon_GC<-GC_porcentaje(secuencias_org_all_genome,6)
Vison_Americano_GC<-GC_porcentaje(secuencias_org_all_genome,7)
Vison_europeo_GC<-GC_porcentaje(secuencias_org_all_genome,8)
Huron_GC<-GC_porcentaje(secuencias_org_all_genome,9)
Hamster_GC<-GC_porcentaje(secuencias_org_all_genome,10)
```


```{r, message=FALSE, results='hide', warning=FALSE}
#Para 200 
x<-1
dat<-list(humano=Humano_GC[[x]], mono=Mono_GC[[x]], perro=Perro_GC[[x]], gato=Gato_GC[[x]],
          tigre=Tigre_GC[[x]], leon=Leon_GC[[x]], vison_am=Vison_Americano_GC[[x]], 
          vison_eu=Vison_europeo_GC[[x]], huron=Huron_GC[[x]] ,hamster=Hamster_GC[[x]])

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(unlist(dat),type="n",xlim=c(1,max(sapply(dat,length))),xlab="Posicion/200", ylab="Porcentaje de GC", main="Porcentaje de GC con una ventana de 200 nucleótidos")
mapply(lines,dat,col=seq_along(dat),lty=2)
legend("topright", inset=c(-0.35,0), names(dat),lty=2,col=seq_along(dat),title="Organismo")
```

```{r, echo=FALSE, message=FALSE, results='hide', warning=FALSE}
#Para 50
x<-2
dat<-list(humano=Humano_GC[[x]], mono=Mono_GC[[x]], perro=Perro_GC[[x]], gato=Gato_GC[[x]],
          tigre=Tigre_GC[[x]], leon=Leon_GC[[x]], vison_am=Vison_Americano_GC[[x]], 
          vison_eu=Vison_europeo_GC[[x]], huron=Huron_GC[[x]] ,hamster=Hamster_GC[[x]])

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(unlist(dat),type="n",xlim=c(1,max(sapply(dat,length))),xlab="Posicion/500",ylab="Porcentaje de GC", main="Porcentaje de GC con una ventana de 500 nucleótidos")
mapply(lines,dat,col=seq_along(dat),lty=2)
legend("topright", inset=c(-0.35,0), names(dat),lty=2,col=seq_along(dat),title="Organismo")
```

```{r,  echo=FALSE, echo=FALSE, message=FALSE, results='hide', warning=FALSE}
#Para 1000
x<-3
dat<-list(humano=Humano_GC[[x]], mono=Mono_GC[[x]], perro=Perro_GC[[x]], gato=Gato_GC[[x]],
          tigre=Tigre_GC[[x]], leon=Leon_GC[[x]], vison_am=Vison_Americano_GC[[x]], 
          vison_eu=Vison_europeo_GC[[x]], huron=Huron_GC[[x]] ,hamster=Hamster_GC[[x]])

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(unlist(dat),type="n",xlim=c(1,max(sapply(dat,length))),xlab="Posicion/1000", ylab="Porcentaje de GC", main="Porcentaje de GC con una ventana de 1000 nucleótidos")
mapply(lines,dat,col=seq_along(dat),lty=2)
legend("topright", inset=c(-0.35,0), legend=names(dat), lty=2,col=seq_along(dat),title="Organismo")
```

##3. *Análisis extras*

Para esta parte se utilizó un paquete llamado CovidMutaciones(). El cual se descargó del repositorio de GitHub (devtools::install_github("MSQ-123/CovidMutations")).



```{r}
library(CovidMutations)
#The example data:
data("nucmer")
head(nucmer)
# Fix IUPAC codes
nucmer<-nucmer[!nucmer$qvar%in%c("B","D","H","K","M","N","R","S","V","W","Y"),]
#nucmer<- mergeEvents(nucmer = nucmer)## This will update the nucmer object
```

+ Anotación completa de todas las mutaciones identificadas por este estudio. Las columnas se describen aquí. Muestra: ID de muestra de GISAID; refpos: posición en el genoma de referencia NC_045512.2 ; refvar: composición de nucleótidos de la referencia en la coordenada refpos (un "." indica una inserción); qvar: variante en la muestra de la consulta (un "." indica una eliminación); qlength: longitud del genoma de consulta (el genoma de referencia siempre tiene una longitud de 29.903 nucleótidos)

Para proporcionar efectos de cada SNP, inserción y deleción en el genoma del virus.
```{r}
data("refseq")
data("gff3")

head(gff3)
annot <- setNames(gff3[,10],gff3[,9])  #annot: subset the gene and its product.
outdir <- tempdir()
head(indelSNP(nucmer = nucmer,
         saveRda = FALSE,
         refseq = refseq,
         gff3 = gff3,
         annot = annot,
         outdir = outdir))
```
+ región: región anotada en la posición del evento (secuencia codificante, intergénica o UTR); variant: un cambio de proteína (mostrado como código de aminoácido) o la posición genómica (si el evento afecta a una región no codificante); varclass: clase variante; annotation: nombre completo de la proteína codificada por la región afectada (si codifica); varname: nombre completo de la variante de proteína; varclade: nombre completo de la variante de nucleótidos.

El paquete también grafica las estadísticas de mutación::
```{r}
data("covid_annot")
head(covid_annot)
covid_annot <- as.data.frame(covid_annot)
outdir <- tempdir()

print("Para ver las graficas revisar en el directorio de su computadora:")
print(tempdir)


plotMutAnno(results = covid_annot,figureType = "MostMut", outdir = outdir)
plotMutAnno(results = covid_annot,figureType = "MutPerSample", outdir = outdir)
plotMutAnno(results = covid_annot,figureType = "VarClasses", outdir = outdir)
plotMutAnno(results = covid_annot,figureType = "VarType", outdir = outdir)
plotMutAnno(results = covid_annot,figureType = "NucleoEvents", outdir = outdir)
plotMutAnno(results = covid_annot,figureType = "ProEvents", outdir = outdir)

```

Ahora le agregaremos la ubicacion geografica de los datos. 
```{r}
data("nucmer")
data("chinalist")
outdir <- tempdir()
head(nucmerRMD(nucmer = nucmer, outdir = outdir, chinalist = chinalist),3)
```


A continuación utilizaremos los datos de aqui para hacer una data frame que tenga los valores 

País | SNP | Frecuencia

```{r}
covid_gisaid<-nucmerRMD(nucmer = nucmer, outdir = outdir, chinalist = chinalist)

sort(table(toString(covid_gisaid$country), decreasing = T))
newdata<-as.data.frame(table(covid_gisaid$rpos,covid_gisaid$country),col.names=c("posicion","pais","frecuencia"))
colnames(newdata)<-c("posicion","pais","frecuencia")
head(newdata,3)
```
Graficando:
```{r, message=FALSE, echo=FALSE, warning=FALSE, message=FALSE}
ggplot(newdata, aes(x=posicion, y=frecuencia, fill=pais))+
  geom_point(position=position_jitterdodge(jitter.width=2, dodge.width = 0),
             pch=21, aes(fill=factor(pais)), show.legend = T)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_discrete(breaks = newdata$posicion[c(T,rep(F,1400))])+
  ggtitle("                            Frecuencia relativa de Mutaciones") +
  xlab("Posición") + ylab("Frecuencia")

```

Realizaremos una analisis anova para ver que son significativamente distintas las frecuencias con respecto a las posiciones en los distintos países. Donde nuestra hipotesis nula significa que las medias de las frecuencias por posicion no son significativamente distintas, mientras que la hipotesis alternativa significa que las medias de las frecuencias de mutaciones por posicion son significativamente distintas. 

```{r}
summary(aov(newdata$frecuencia~newdata$posicion))
```

Se rechaza la hipotesis nula, las medias son significativamente distintas. 
```{r}
summary(aov(newdata$frecuencia~newdata$pais))
```
Se rechaza la hipotesis nula, que es análoga a la anterior solo que y no analiza las posiciones, sino por países. 

Por último se realiza una anova de dos vías. Recortado ya que si se considera todas las secuencias el análisis es demasiado grande. 
```{r}
#summary(aov(newdata$frecuencia[1:20000]~newdata$pais[1:20000]*newdata$p#osicion[1:20000]))
```

