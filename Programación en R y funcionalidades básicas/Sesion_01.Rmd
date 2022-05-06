---
title: "Session 2 – Using R, Installing Packages and Importing/Exporting Data"
output: html_document
date: '2022-05-06'
---
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Tu primera sesión

Partes de R Studio 
Hacer click en View > Panes 
Empezar usando la consola 


### 1.1 Realización de algunas funciones básicas y obtención de ayuda

```{r}
x <- 5
y <- 3
z <- x*y
A*x
A*y
A*(x - y)
res <- (x^2)*y + A 
range (res)
mean (res)
median(res)
min(res)
max(res)
```

Obteniendo ayuda 

```{r}
# get help about function 'mean'
?mean
# another way to ask help
help(mean)
```

Las funciones ?y help()proporcionan un enlace a una nueva ventana con el manual correspondiente de esa función y cómo se usa y (en la mayoría de los casos) ejemplos.

Algunas funciones básicas para vectores numéricos.

```{r}
# Arithmetic Mean
mean(my_vector) 
# Standard Deviation
sd (my_vector) 
# Variance
var (my_vector) 
# Median Value
median (my_vector) 
# Maximum Value
max (my_vector)
# Minimum Value 
min (my_vector)
# Sums values in vector
sum(my_vector)
# Provides the number of elements in the vector 
length(my_vector) 
# Round numbers to of elements in the vector to number of digits
round(3.1415, digits = 2)
```

### 1.2 The working directory 

- Enumerar todos los objetos "en la memoria", es decir, mostrar los nombres de los objetos en su espacio de trabajo usando la función ls()

- También puede enumerar objetos y su estructura con la función ls.str().
- También puede eliminar objetos de su espacio de trabajo con la función rm()
- En el caso de que sólo quieras guardar algunos objetos de tu wokspace puedes usar la función save de R.
- Sin embargo, cuando se quiere guardar un solo objeto es conveniente utilizar la función saveRDS que guardará los datos en el formato RDS.

```{r}
# Borramos todo lo que haya en la memoria
rm(list = ls())

# Creamos unas variables
x <- 20
y <- 34
z <- "casa"

# Guardamos los objetos en mis_Objetos.RData
save.image(file = "mis_Objetos.RData")

# Guardando los objetos 'x' e 'y'
save(x, y, file = "mis_dos_Objetos.RData") 

saveRDS(x, file = "mi_objeto.rds")
```

## 2. Instalación de paquetes

CRAN es el repositorio oficial de paquetes de R, que cuenta con miles de paquetes de R gratuitos. La mayoría de ellos han sido desarrollados por científicos de datos, estadísticos, profesores universitarios e investigadores.

```{r}
install.packages("ggplot2")
library(ggplot2)
install.packages("devtools")
library(devtools)

#require(devtools)
#install_version("ggplot2", version = "0.9.1", repos = "http://cran.us.r-project.org")
```

La mayoría de los paquetes R genéticos, genómicos y transcriptómicos son bioconductores . Hasta la fecha (26/01/22), este archivo tiene 2083 paquetes. Para instalar los paquetes presentes en este archivo, debe seguir las siguientes instrucciones.

```{r}
## test if 'BiocManager' is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
BiocManager::install("DECIPHER")
```

```{r}
# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("tidyverse/ggplot2")
```


## 3. Estructuras de datos 

### 3.1 Tipos de datos

Numeric 
```{r}
my_numeric_vector <- c(1,3,45,56,1)
my_numeric_vector
my_numeric_vector <- 1:10
my_numeric_vector
```


Strings
```{r}
my_mix_vector <- c("my", "bioinformatics", "class", 3.141593)
my_mix_vector
#[1] "my"             "bioinformatics" "class"          "3.141593"      
str(my_mix_vector)
```

Factors
```{r}
my_factor_vector <- c(1,0,1,1,1,0)
my_factor_vector <- factor(my_factor_vector)
my_factor_vector
```

Logical 
```{r}
my_logical_vector <- ifelse(my_numeric_vector > 5, TRUE, FALSE)
my_logical_vector
```


### 3.2 Importar datos
```{r}
library(readxl)
getwd()

# XLSX
Data_procesada <- read_excel("Data procesada.xlsx", sheet = "Hoja 1")
#as.data.frame(Data_procesada)

# CSV
write.csv(Data_procesada,"Mi_primer_archivo_guardado.csv", row.names = FALSE)
Mi_primer_archivo_guardado <- read.csv("Mi_primer_archivo_guardado.csv")

#TXT
write.table(Data_procesada,"Mi_primer_archivo_guardado.txt", append = FALSE, sep = " ", dec = ".",  row.names = FALSE, col.names = TRUE)
Mi_primer_archivo_guardado <- read.csv("Mi_primer_archivo_guardado.txt", sep="")
```

Funciones basicas con data frames 
```{r}
head(Data_procesada)
str(Data_procesada)
```
```{r}
subset(Data_procesada, select = `Tasa de infeccion`)
subset(Data_procesada, Vacunados  >= 3)
```

```{r}
my_dataframe_2 <- cbind(Data_procesada,Data_procesada$FECHA, Data_procesada$Positivos)
my_dataframe_2
table(Data_procesada$FECHA, Data_procesada$Positivos)
```

### 4 Secuencias 
```{r}
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
```


```{r}
## Instalando paquetes con dependencias
install.packages("seqinr", dependencies = T)
install.packages("ape", dependencies = T)
library(seqinr)
library(ape)
```

```{r}
## 1. Cargando una secuencia manualmente ====

path <- "sequence.fasta"
actin.necator <- seqinr::read.fasta(path, seqtype = "AA")

attr(actin.necator$XP_013298257.1,"Annot")
tabla <- table(actin.necator$XP_013298257.1)
sort(table(actin.necator$XP_013298257.1))
sort(table(actin.necator$XP_013298257.1),decreasing = T)
plot(sort(table(actin.necator$XP_013298257.1),decreasing = T), col=colores1)
barplot(sort(table(actin.necator$XP_013298257.1),decreasing = T), col=colores1)
longitud <- length(actin.necator$XP_013298257.1)
freqre <- tabla/longitud
barplot(freqre)

#2. Cargando multiples secuencias ====

path <- "app_Hs_multiple_seq.txt"
app.Hs <- seqinr::read.fasta(path, as.string = T, seqtype = "AA")
```

```{r}

#3. Importando secuencias desde R ====

seqinr::choosebank(infobank = T)
seqinr::choosebank("swissprot")
app2 <- query("app2", "K=Amyloid AND sp=Homo sapiens")
summary(app2)
app2$call
app2$nelem
app2$req
app2$req[[2]]
seqinr::getSequence(app2$req[[2]])
```

```{r}
#4. Usando ape ====

acs.num <- c("MT419843.1","MT450923.1","MT419850.1",
             "MT419838.1","MT419849.1","MT419840.1")
library(ape)
cov2 <- ape::read.GenBank(acs.num, as.character = T)
cov2
```

```{r}
#Crear pdf con los graficos
pdf(file = paste0(getwd(),"/comparacion.pdf"), paper = "US")
par(mfrow = c(3,2))

for (gen in 1:length(cov2)) {
  print(table(cov2[gen]))
  barplot(round(table(cov2[[gen]])*100/length(cov2[[gen]]), 1), main = names(cov2)[gen], las = 1, col=colores1)
}
dev.off()
```


