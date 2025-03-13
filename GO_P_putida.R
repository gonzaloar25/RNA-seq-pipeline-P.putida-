#Gene Ontology
# Instalar topGO si no está disponible
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("topGO")

# Cargar bibliotecas necesarias
library(topGO)

# Cargar el archivo genes_ontology.csv, descargado de la base de datos Pseudomonas genome DB
archivo_genes <- "/media/gonzalo/EXTERNAL_USB/gene_ontology_PP.csv"  # Ajusta la ruta
genes_ontology <- read.csv(archivo_genes, header = TRUE, stringsAsFactors = FALSE)

# Crear el mapeo de genes a términos GO
geneID2GO <- split(genes_ontology$Accession, genes_ontology$Locus.Tag)

# Ver el resultado del mapeo
head(geneID2GO)

# Inspeccionar las primeras filas
head(genes_ontology$Locus.Tag)

# Crear un factor para los genes de interés
gene_list <- as.integer(genes_ontology$Locus.Tag %in% DEG)  # Asegúrate de que 'DEG' esté definido, del script DEG_analysis
names(gene_list) <- genes_ontology$Locus.Tag
print(gene_list)

# Filtrar los genes con logFC > 1 o logFC < -1 y FDR <= 0.05
genes_filtrados <- res_corrected$table[res_corrected$table$FDR <= 0.05 & abs(res_corrected$table$logFC) >= 1, ]
genes_filtrados$Locus_tag <- rownames(genes_filtrados)
head(genes_filtrados)

# Extraer los nombres de los genes filtrados
gene_names <- rownames(genes_filtrados)

# Crear un vector con los valores de logFC para todos los genes
all_logFC <- res_corrected$table$logFC

# Asignar los nombres de los genes (utilizando las filas de res_corrected$table)
names(all_logFC) <- rownames(res_corrected$table)
head(all_logFC)

# Asignar el FDR y que los nombres de los genes coincidan entre all_logFC y res_corrected$table
FDR_values <- res_corrected$table$FDR
names(FDR_values) <- rownames(res_corrected$table)

selected_genes <- (abs(all_logFC) > 1) & (FDR_values < 0.05)
table(selected_genes)


# Crear la función de selección de genes significativos
geneSelectionFun <- function(allGenes) {
  # Seleccionar los genes significativos: logFC > 1 o logFC < -1, y FDR < 0.05
  selectedGenes <- (abs(allGenes) > 1) & (res_corrected$table$FDR < 0.05)
  return(selectedGenes)
}

# Crear el objeto topGOdata
GOdata <- new("topGOdata",
              ontology = "BP",  # Puedes cambiar a "MF" o "CC" según el tipo de análisis
              allGenes = all_logFC,
              annot = annFUN.gene2GO,
              gene2GO = geneID2GO,  # Asegúrate de tener este mapeo correctamente definido
              nodeSize = 10,  # Número mínimo de genes en un nodo
              geneSelectionFun = geneSelectionFun)

# Verificar el objeto creado
GOdata


# Crear un dataframe para almacenar los resultados
genes_filtrados_df <- data.frame(
  Locus_tag = rownames(genes_filtrados),  # Asegúrate de que genes_filtrados tenga la columna 'Locus_tag'
  name = NA,
  description = NA,
  GO.MF = NA,
  GO.BP = NA,
  GO.CC = NA,
  # Nueva columna para las descripciones
  stringsAsFactors = FALSE
)

# Asegúrate de que `genes_ontology` tiene las columnas correctas. Aquí usas "Locus.Tag", "GO.Term", "Namespace", "Gene.Name", y "Product.description"
genes_filtrados_df <- genes_ontology %>%
  filter(Locus.Tag %in% genes_filtrados$Locus_tag) %>%
  group_by(Locus.Tag) %>%
  summarise(
    name = paste(unique(Gene.Name[!is.na(Gene.Name)]), collapse = "; "),
    description = paste(unique(Product.Description[!is.na(Product.Description)]), collapse = "; "),
    GO.BP = paste(unique(GO.Term[Namespace == "biological_process"]), collapse = "; "),
    GO.MF = paste(unique(GO.Term[Namespace == "molecular_function"]), collapse = "; "),
    GO.CC = paste(unique(GO.Term[Namespace == "cellular_component"]), collapse = "; ")
  ) %>%
  ungroup()
# Guardar los resultados en un archivo Excel y CSV
#install.packages("writexl")  # Solo si no lo has instalado previamente
library(writexl)

write.csv(genes_filtrados_df, file = "GOs_Assigned_to_Genes.csv", row.names = FALSE, quote = TRUE)

write.xlsx(genes_filtrados_df, file = "GOs_Assigned_to_Genes.xlsx", rowNames = FALSE)

### Anotación automática de términos en archivos csv de DEGs

# Cargamos el archivo con todos los términos Attributes de KT2440 de PseudomonasDB
terms <- read.csv(file = "GOs_Assigned_to_Genes.csv", header = TRUE)
head(terms)

# Cargamos cada archivo de DEGs
DEGs <- read.csv(file = "DEG_WTvsfleQ.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Anotamos los términos en la tabla de DEGs 
DEGs.genes <- as.vector(DEGs$Geneid)
terms[1,1]
i <- 1
PP.temp <- "empty"

for (i in 1:nrow(terms)) {
  # Si el gen está en los genes de DEGs
  if ((terms[i, 1]) %in% DEGs.genes) {
    PP.temp <- as.character(terms[i, 1])  # Locus Tag
    Name.temp <- as.character(terms[i, 2])  # Nombre
    Description.temp <- as.character(terms[i, 3])  # Descripción
    GO.BP.temp <- as.character(terms[i, 4])  # GO Biological Process
    GO.MF.temp <- as.character(terms[i, 5])  # GO Molecular Function
    GO.CC.temp <- as.character(terms[i, 6])  # GO Cellular Component
    
    # Buscamos el número de la fila en DEGs donde el 'Geneid' coincide
    temp.row.number <- which(DEGs$Geneid == PP.temp)
    
    # Añadimos la información en las columnas correspondientes
    if (length(temp.row.number) > 0) {  # Si se encuentra la fila
      DEGs[temp.row.number, 8] <- Name.temp
      DEGs[temp.row.number, 9] <- Description.temp
      DEGs[temp.row.number, 10] <- GO.BP.temp
      DEGs[temp.row.number, 11] <- GO.MF.temp
      DEGs[temp.row.number, 12] <- GO.CC.temp
    }
  }
}

# Renombrar las columnas 8 a 12
colnames(DEGs)[8:12] <- c("Name", "Description", "GO.BP", "GO.MF", "GO.CC")
colnames(DEGs)

# Guardamos los datos
write.csv(DEGs, file = "ResultadosAnotados_WT_vs_fleQ.csv", row.names = FALSE, quote = TRUE)
write.xlsx(DEGs, file = "ResultadosAnotados_WT_vs_fleQ.xlsx", rownames = FALSE)

#Enriquecimiento
DEGs.file <- read.csv ("ResultadosAnotados_WT_vs_fleQ.csv")

# Coger los IDs de los genes de las filas seleccionadas
genes_ID <- DEGs.file[, 1]

# Get background annotation
geneID2GO <- split(genes_ontology$Accession, genes_ontology$Locus.Tag)
background_genes <- names(geneID2GO)

head(GO2geneID)
head(geneID2GO)

# Compare genes vs bg_genes
compared_genes <- factor(as.integer(background_genes %in% genes_ID))
names(compared_genes) <- background_genes

# Create topGO object
GOdata <- new("topGOdata",
              ontology = "BP",  # Puedes cambiar a "MF" o "CC" según el tipo de análisis
              allGenes = all_logFC,
              annot = annFUN.gene2GO,
              gene2GO = geneID2GO,  # Asegúrate de tener este mapeo correctamente definido
              nodeSize = 10,  # Número mínimo de genes en un nodo
              geneSelectionFun = geneSelectionFun)

# Verificar el objeto creado
GOdata

#Fisher test
# El modo "elim", en lugar del "classic", elimina la redundancia en la jerarquía de términos GO.
resultFisher <- runTest(GOdata, algorithm = "elim" , statistic = "fisher")

#Crear e imprimir la tabla del enriquecimiento
Nodes <- 40 #Puede cambiarse
allRes <- GenTable(GOdata, classicFisher = resultFisher, topNodes = Nodes)
allRes

#Añadir el valor de classicFisher
allRes$classicFisher <- as.numeric(allRes$classicFisher)

#Añadir el p-valor ajustado
p_values <- allRes$classicFisher
p_value_adjusted <- p.adjust(p_values, method = "BH")
allRes$p_value_adjusted <- p_value_adjusted

# Calcular el Gene ratio
allRes$GenRatio <- allRes$Significant / allRes$Annotated

#Crear el gráfico de dispersión con los 15 términos más representativos
top_15 <- allRes %>% top_n(15, GenRatio)

ggplot(top_15, aes(x = GenRatio, 
                   y = reorder(Term, GenRatio), 
                   size = Annotated, 
                   color = p_value_adjusted)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  labs(x = "Gen Ratio", y = "Término GO CC", 
       size = "Genes Asociados", color = "p_valor adj", title = "Genes regulador por fleQ") +
  theme_minimal()

write.xlsx(allRes, file = "GO_enrichment_dfleQ_CC.xlsx")