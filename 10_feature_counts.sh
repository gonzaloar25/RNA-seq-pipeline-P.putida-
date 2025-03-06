#!/bin/bash

# Directorio donde están las carpetas de entrada
directorio="/media/gonzalo/EXTERNAL_USB/rnaseq/datos/processed/4_Mark_duplicates"

# Directorio donde se guardarán los archivos de conteo
output_dir="/media/gonzalo/EXTERNAL_USB/rnaseq/datos/processed/5_Gene_counts"

# Verifica si el directorio de salida existe, si no, lo crea
if [ ! -d "$output_dir" ]; then
    echo "Creando directorio de salida: $output_dir"
    mkdir -p "$output_dir"
else
    echo "El directorio de salida ya existe: $output_dir"
fi

# Recorre todas las subcarpetas
for carpeta in "$directorio"/*/; do

    # Encuentra el archivo BAM en la carpeta (contiene .bam en su nombre)
    archivo_BAM=$(find "$carpeta" -type f -name "*.bam")

    if [ -n "$archivo_BAM" ]; then  
        # Extrae el nombre base del archivo sin la ruta y la extensión
        nombre_base=$(basename "$archivo_BAM" "_deduplicated.bam")

        # Archivo de salida en el directorio de conteo
        archivo_salida_tsv="$output_dir/${nombre_base}_counts.tsv"

        # Ejecuta featureCounts si aún no se había hecho anteriormente
        if [ -f "$archivo_salida_tsv" ]; then
            echo "Archivo de salida ya existe para $nombre_base. Omitiendo cuantificación de genes."
        else
            echo "Ejecutando featureCounts en el archivo $archivo_BAM"
            featureCounts -t CDS,tRNA,ncRNA,tmRNA,rRNA -p -g locus -a /media/gonzalo/EXTERNAL_USB/rnaseq/genome/Pseudomonas_putida_KT2440_definitivo.gff -o $archivo_salida_tsv $archivo_BAM 
        fi
    else
        echo "No se puede ejecutar el comando: falta archivo BAM en $carpeta"
    fi
done


