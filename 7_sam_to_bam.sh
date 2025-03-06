#!/bin/bash

# Directorio donde están las carpetas de entrada
directorio="/media/gonzalo/EXTERNAL_USB/rnaseq/datos/processed/3_Alignment"

# Recorre todas las subcarpetas
for carpeta in "$directorio"/*/; do

    # Verifica si ya existe un archivo BAM en la carpeta
    archivo_BAM=$(find "$carpeta" -type f -name "*.bam")
    
    if [ -n "$archivo_BAM" ]; then
        echo "Se encontró un archivo BAM en $carpeta. Omitiendo esta carpeta."
        continue
    fi

    # Encuentra el archivo SAM en la carpeta (contiene .sam en su nombre)
    archivo_SAM=$(find "$carpeta" -type f -name "*.sam")
    
    if [ -n  "$archivo_SAM" ]; then  
        # Extrae el nombre base del archivo sin la ruta y la extensión
        nombre_base=$(basename "$archivo_SAM" ".sam")

        echo "Ejecutando conversión del archivo $archivo_SAM a formato BAM"
        samtools view -@ 10 -Sbh "$archivo_SAM" > "$carpeta/$nombre_base.bam"
        samtools sort -@ 10 "$carpeta/$nombre_base.bam" -o "$carpeta/$nombre_base.sorted.bam"
        samtools index -@ 9 "$carpeta/$nombre_base.sorted.bam"
    else
        echo "No se puedo ejecutar el comando: falta archivo SAM en $carpeta"
    fi
done
