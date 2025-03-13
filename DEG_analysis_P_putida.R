# RNA-seq analysis in R
# 5 November 2024
# Author: Gonzalo Ariza

#Data files

#Libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("curl")
install.packages(c("httr", "UCSC.utils", "GenomeInfoDb", "Biostrings", "KEGGREST", "AnnotationDbi"))

BiocManager::install("edgeR")
install.packages("ggplot2")
install.packages("statmod")
install.packages(c("dplyr", "tibble"))
install.packages("ggrepel")

library(edgeR)
library(ggplot2)
library(statmod)
library(dplyr)
library(tibble)
library(ggrepel)

#Establecer directorio de trabajo
setwd(dir = "E:/tfm/r_analysis")

# 1. Importar metadatos, puede ser usando un archivo preexistente o definiendolos como en este caso 
# Definir los mutantes y las condiciones
group <- factor(rep(c("fleN", "fleN_lapA_pleD", "fleN_lapA_yhjH", "fleN_lapA", 
                      "fleQ", "lapA_pleD", "lapA_yhjH", "lapA", "WT"), each = 3), 
                levels = c("WT", "fleN", "fleN_lapA_pleD", "fleN_lapA_yhjH", 
                           "fleN_lapA", "fleQ", "lapA_pleD", "lapA_yhjH", "lapA"))
print(table(group))

#Establecemos el directorio de los conteos con los archivos .tsv y guardamos los archivos en una lista
directorio_tsv = "/media/gonzalo/EXTERNAL_USB/tfm/datos/processed/6_Gene_counts/"
archivos <- list.files(path = directorio_tsv, pattern = "\\.tsv$", full.names = TRUE)
head(archivos)

#Leemos y procesamos los archivos, elimando las 2 primeras columnas y estableciendo el nombre de cada fila según su Geneid (PP_XXXX)
seqdata <- lapply(archivos, function(file) {
  datos <- read.table(file, row.names = "Geneid", header = TRUE, sep = "\t")
  datos <- datos[, 2:ncol(datos)]
  colnames(datos)[ncol(datos)] <- "Count"  # Renombrar la última columna a "Count"
  return(datos)
})
head(seqdata[[1]])

#2. Convertir los counts en el objeto DGEList
# Combinar los conteos en una matriz
conteos_matriz <- do.call(cbind, lapply(seqdata, `[[`, "Count"))

# Asegurarse de que los nombres de las filas sean los nombres de los genes
rownames(conteos_matriz) <- rownames(seqdata[[1]])
head(conteos_matriz)

#Creamos una variable que contenga el nombre de cada condición
condition <- sub("_counts\\.tsv$", "", basename(archivos))  # Limpiar sufijo '_counts.tsv'
colnames(conteos_matriz) <- condition

# Verifica que el número de columnas de conteos_matriz sea igual al número de elementos en condition
if (ncol(conteos_matriz) != length(group)) {
  stop("El número de condiciones no coincide con el número de columnas en la matriz de conteos.")
}

# Visualización de verificación
head(conteos_matriz)

# Crear el objeto DGEList
y <- DGEList(counts = conteos_matriz)

# Verifica el objeto DGEList
print(y)
print(colnames(y))
y$samples$group <- group
print(y$samples)

#3. Filtrado para quitar los conteos bajos
keep <- filterByExpr(y)
head(keep)
summary(keep)
y <- y[keep, keep.lib.sizes=FALSE]

# Ver los primeros genes que no pasan el filtro
genes_no_pasan <- rownames(y)[!keep]
print(genes_no_pasan)

length(y$samples$group)  # Número de muestras en DGEList
length(keep) 

#Porcentaje de genes que han pasado el filtro de baja expresión
total_genes <- nrow(conteos_matriz)
genes_filtrados <- sum(keep)
porcentaje_genes_filtrados <- (genes_filtrados / total_genes) * 100
print(porcentaje_genes_filtrados)

#boxplot logCPM antes de la normalización
logCPM <- cpm(y, log = TRUE)

#Representamos solo algunas condiciones para que no sea muy sobrecargado
par(mar = c(8, 4, 4, 2))
boxplot(logCPM[, c(1,2,3,13,14,15,25,26,27)], 
        main = "Recuento de lecturas antes de la normalización",
        ylab = "logCPM",
        las = 2,
        outline = FALSE,           
        cex.axis = 0.8,            
        cex.lab = 1.2,             
        cex.main = 1.4)
mtext("Muestras", side = 1, line = 5)

#4. Normalizacion
y <- calcNormFactors(y)
head(y)

logCPM_norm <- cpm(y, log = TRUE)

boxplot(logCPM_norm[, c(1,2,3,13,14,15,25,26,27)], 
        main = "Recuento de lecturas después de la normalización",
        ylab = "logCPM",
        las = 2,                  
        outline = FALSE,           
        cex.axis = 0.8,            
        cex.lab = 1.2,             
        cex.main = 1.4)
mtext("Muestras", side = 1, line = 5)

#Muestras que tienen un factor de normalización menor que 1
factores_norm <- y$samples$norm.factors
muestras_menor_1 <- rownames(y$samples)[factores_norm < 1]
muestras_menor_1

#PlotMDS (para observar la distribución de las muestras y la variabilidad biológica)
plotMDS(y)
pch <- c(16,16,16,15,17,16,16,15,17)
colors <- c("darkgreen","orange", rep("red", 3), "blue", rep("magenta", 3))
plotMDS(y, col=colors[group], pch=pch[group])
legend("topright", legend=levels(group), pch=pch, col=colors, ncol=1, cex = 0.6)
title(main = "Plot MDS muestras RNA-seq")

#5. Diseño de la matriz modelo
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
head(design)

# Verifica los nuevos nombres
print(colnames(design))

#6. Estimación de las dispersiones
y <- estimateDisp(y, design, robust = TRUE)
y$common.dispersion
y$trended.dispersion
y$tagwise.dispersion
plotBCV(y)
y$samples

#7. Ajuste con GLM Model
fit <- glmQLFit(y, design, robust = TRUE)
head(fit)
plotQLDisp(fit)

#8. Probando expresión diferencial WT vs fleQ, el primer grupo es el testado y el segundo hace de control
WTvsfleQ <- makeContrasts(fleQ - WT, levels = design)
WTvsfleQ

res <- glmQLFTest(fit, contrast = WTvsfleQ)

class(res)
names(res)
head(res$table)

#Para saber el numero de genes (o tests estadisticos realizados):
dim(res$table)

res_corrected <- topTags(res, n = Inf) #La función topTags por defecto reporta los 10 top genes, ponemos n = Inf para que devuelva todos)
dim(res_corrected)
head(res_corrected)

#Filtrado de genes (FDR)
#Filtrado de genes (FDR+logFC)
nrow(res_corrected$table[res_corrected$table$FDR <= 0.05 & abs(res_corrected$table$logFC) >= 1,])
nrow(res_corrected$table[res_corrected$table$FDR <= 0.05 & res_corrected$table$logFC >= 1,])
nrow(res_corrected$table[res_corrected$table$FDR <= 0.05 & res_corrected$table$logFC <= (-1),])

#Otra forma:
is.del <- decideTests(res, adjust.method = "BH", p.value = 0.05, lfc = 1)
head(is.del)
table(is.del)
summary(is.del)

#Visualización
par(mfrow =c(1,2))
plotMD(res, status = is.del)

#Tabla de datos
data <- res_corrected$table

#Añadir nombre Geneid a la primera columna
data <- data %>% 
  rownames_to_column(var = "Geneid")

#Ordenar los genes de mayor a menor según su logFC
data <- data %>% 
  arrange(desc(logFC))

#Añadiendo columna de DE (expresion diferencial)
data$DE <- "NO"
head(data)
data$DE[data$logFC > 1 & data$FDR <= 0.05] <- "UP"
data$DE[data$logFC < -1 & data$FDR <= 0.05] <- "DOWN"
head(data)
nrow(data)  

sum(data$FDR <= 0.05 & abs(data$logFC) > 1)
head(data[data$FDR <= 0.05 & abs(data$logFC) > 1, ])

#Guardamos esos genes para centrarnos en ellos
DEG <- data$Geneid[data$FDR <= 0.05 & abs(data$logFC) > 1]
head(DEG)
length(DEG)

#Guardar resultados en excel
write.table(as.data.frame(data), file = "Resultados_WTvsfleQ.csv", sep = ";", dec = ",", row.names = FALSE)
write.table(as.data.frame(data[data$FDR <= 0.05 & abs(data$logFC) > 1, ]), file = "DEG_WTvsfleQ.csv", sep = ";", dec = ",", row.names = FALSE)

# Seleccionar los primeros 7 genes con logFC > 1
top_genes_pos <- data[data$logFC > 1 & data$FDR <= 0.05, ]
top_genes_pos <- top_genes_pos[order(-top_genes_pos$logFC), ][1:7, ]

# Seleccionar los primeros 7 genes con logFC < -1
top_genes_neg <- data[data$logFC < -1 & data$FDR <= 0.05, ]
top_genes_neg <- top_genes_neg[order(top_genes_neg$logFC), ][1:7, ]

# Combinar ambos resultados
top_genes <- rbind(top_genes_pos, top_genes_neg)

#Creamos el gŕafico de Volcán donde se representen dichos genes
ggplot(data, aes(x = logFC, y = -log10(FDR), col = DE)) + 
  geom_point(size = 2) + 
  theme_minimal() + 
  geom_vline(xintercept = c(-1, 1), col = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "green", linetype = "dashed") +
  geom_text_repel(data = top_genes, aes(label = Geneid), size = 5, 
                  fontface = "bold", box.padding = 0.3, 
                  point.padding = 0.5, segment.color = "gray", max.overlaps = 20) + 
  ggtitle("lapA_yhjH vs lapA_pleD") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "blue"))

#Otra comparacion
WTvsfleN <- makeContrasts(fleN - WT, levels = design)
WTvsfleN

res_2 <- glmQLFTest(fit, contrast = WTvsfleN)
res_corrected_2 <- topTags(res_2, n = Inf)

nrow(res_corrected_2$table[res_corrected_2$table$FDR <= 0.05 & abs(res_corrected_2$table$logFC) >= 1,])
nrow(res_corrected_2$table[res_corrected_2$table$FDR <= 0.05 & res_corrected_2$table$logFC >= 1,])
nrow(res_corrected_2$table[res_corrected_2$table$FDR <= 0.05 & res_corrected_2$table$logFC <= (-1),])
is.del_2 <- decideTests(res_2, adjust.method = "BH", p.value = 0.05, lfc = 1)
summary(is.del_2)

#Tabla de datos
data_2 <- res_corrected_2$table

#Añadir nombre Geneid a la primera columna
data_2 <- data_2 %>% 
  rownames_to_column(var = "Geneid")

#Ordenar los genes de mayor a menor según su logFC
data_2 <- data_2 %>% 
  arrange(desc(logFC))

#Añadiendo columna de DE
data_2$DE <- "NO"
head(data_2)

data_2$DE[data_2$logFC > 1 & data_2$FDR <= 0.05] <- "UP"
data_2$DE[data_2$logFC < -1 & data_2$FDR <= 0.05] <- "DOWN"
head(data_2)
tail(data_2)
nrow(data_2)
print(data)

sum(data_2$FDR <= 0.05 & abs(data_2$logFC) > 1)
head(data_2[data_2$FDR <= 0.05 & abs(data_2$logFC) > 1, ])

#Guardar resultados en excel
write.table(as.data.frame(data_2), file = "Resultados_WTvsfleN.csv", sep = "\t", dec = ",", row.names = FALSE)
write.table(as.data.frame(data_2[data_2$FDR <= 0.05 & abs(data_2$logFC) > 1, ]), file = "Resultados_DEG_WTvsfleN.csv", sep = "\t", dec = ",", row.names = FALSE)

#Otra comparacion
WTvslapA <- makeContrasts(lapA_pleD - lapA, levels = design)
WTvslapA

res_3 <- glmQLFTest(fit, contrast = WTvslapA)
res_corrected_3 <- topTags(res_3, n = Inf)

nrow(res_corrected_3$table[res_corrected_3$table$FDR <= 0.05 & abs(res_corrected_3$table$logFC) >= 1,])
nrow(res_corrected_3$table[res_corrected_3$table$FDR <= 0.05 & res_corrected_3$table$logFC >= 1,])
nrow(res_corrected_3$table[res_corrected_3$table$FDR <= 0.05 & res_corrected_3$table$logFC <= (-1),])
is.del_3 <- decideTests(res_3, adjust.method = "BH", p.value = 0.05, lfc = 1)
summary(is.del_3)

#Tabla de datos
data_3 <- res_corrected_3$table

#Añadir nombre Geneid a la primera columna
data_3 <- data_3 %>% 
  rownames_to_column(var = "Geneid")

#Ordenar los genes de mayor a menor según su logFC
data_3 <- data_3 %>% 
  arrange(desc(logFC))

#Añadiendo columna de DE
data_3$DE <- "NO"
head(data_3)
data_3$DE[data_3$logFC > 1 & data_3$FDR <= 0.05] <- "UP"
data_3$DE[data_3$logFC < -1 & data_3$FDR <= 0.05] <- "DOWN"
head(data_3)
tail(data_3)
nrow(data_3)
print(data)

sum(data_3$FDR <= 0.05 & abs(data_3$logFC) > 1)
head(data_3[data_3$FDR <= 0.05 & abs(data_3$logFC) > 1, ])

DEG <- data$Geneid[data$FDR <= 0.05 & abs(data$logFC) > 1]
DEG_2 <- data_2$Geneid[data_2$FDR <= 0.05 & abs(data_2$logFC) > 1]
DEG_3 <- data_3$Geneid[data_3$FDR <= 0.05 & abs(data_3$logFC) > 1]
genes_comunes <- Reduce(intersect, list(DEG, DEG_2, DEG_3))
print(genes_comunes)

#Variables para agrupar los genes según su DE
gene_fleq_up <- data$Geneid[data$DE == "UP"]
gene_fleq_down <- data$Geneid[data$DE == "DOWN"]
gene_fleq_no <- data$Geneid[data$DE == "NO"]

gene_flen_up <- data_2$Geneid[data_2$DE == "UP"]
gene_flen_down <- data_2$Geneid[data_2$DE == "DOWN"]
gene_flen_no <- data_2$Geneid[data_2$DE == "NO"]

gene_lapa_up <- data_3$Geneid[data_3$DE == "UP"]
gene_lapa_down <- data_3$Geneid[data_3$DE == "DOWN"]
gene_lapa_no <- data_3$Geneid[data_3$DE == "NO"]


install.packages("VennDiagram")
library(VennDiagram)

# Creamos una lista para el diagrama de Venn
venn_genes <- list(
  "Reprimidos por fleQ" = gene_fleq_up,
  "Activados por fleQ" = gene_fleq_down,
  "No significativos fleQ" = gene_fleq_no,
  "Reprimidos por fleN" = gene_flen_up,
  "Activados por fleN" = gene_flen_down,
  "No significativos fleN" = gene_flen_no,
  "Activados por diGMPc" = gene_lapa_up,
  "Reprimidos por diGMPc" = gene_lapa_down,
  "No significativos diGMPc" = gene_lapa_no
)

# Crear el diagrama de Venn
venn.plot <- venn.diagram(
  x = venn_genes,
  filename = NULL,  #Para mostrar el gráfico en la ventana gráfica
  fill = c("lightblue", "pink", "lightgreen", "lightyellow"),
  alpha = 0.7,  #Transparencia
  cex = 1.5,  #Tamaño del texto
  cat.cex = 1.5,  # Tamaño de las etiquetas
  main = "Diagrama de Venn de genes con expresión diferencial significativa"
)
grid.draw(venn.plot)

if (!requireNamespace("ComplexUpset", quietly = TRUE)) {
  install.packages("ComplexUpset")
}

#Creacion de Upset plot
library(ComplexUpset)
library(ggplot2)

genes <- unique(unlist(venn_genes)) # Todos los genes únicos
binary_matrix <- as.data.frame(sapply(venn_genes, function(x) genes %in% x))
colnames(binary_matrix) <- names(venn_genes)
binary_matrix <- cbind(Gene = genes, binary_matrix)
sets <- colnames(binary_matrix)[-1]
print(binary_matrix)

#Podemos guardar los resultados en un archivo
write.csv(binary_matrix, file = "patrones_expresion.csv", row.names = FALSE, quote = TRUE)

#Representamos el Upset plot
upset(binary_matrix, 
      intersect = colnames(binary_matrix)[-1],  # Especifica las columnas que representan los conjuntos
      base_annotations = list(
        'Intersection Size' = intersection_size(counts = TRUE)),  # Muestra el tamaño de la intersección 
      ggtitle("UpSet Plot de genes con expresión diferencial en dfleQ, dfleN y dlapA frente a WT"))


#Pasos para crear un heatmap

#Transformamos los datos
logCPM <- cpm(y, log = TRUE)
head(logCPM)
head(DEG)
colnames(logCPM)

#Filtramos la matriz total con los genes de DE en las condiciones de interes
DEG_selection <- logCPM[rownames(logCPM) %in% DEG, c(1,2,3,13,14,15,22,23,24)]
head(DEG_selection)
nrow(DEG_selection)

#Creación de heatmap
install.packages("pheatmap")
library(pheatmap)
pheatmap(DEG_selection, scale = "row", cluster_rows = F, cluster_cols = T, 
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         clustering_method = "average",
         display_numbers = F,
         fontsize_number = 6,
         fontsize_row = 7,
         border_color = NA,
         show_rownames = FALSE,
         main = "Heatmap usando el aglomeramiento average con DEG WT vs fleQ")
