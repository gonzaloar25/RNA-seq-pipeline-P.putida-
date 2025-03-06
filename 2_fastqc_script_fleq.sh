#!/bin/bash

# Directorio donde están las carpetas
directorio="/media/gonzalo/EXTERNAL_USB/tfm/datos/fleq/raw_data"

# Directorio donde redirigimos la salida del comando
output_dir="/media/gonzalo/EXTERNAL_USB/tfm/datos/processed/1_Quality/FastQC"

# Recorre todos los archivos fastq.gz
for archivo in "$directorio"/*_R1*.fastq.gz; do
    # Verifica que haya encontrado al menos un archivo R1
    if [ -f "$archivo" ]; then
        # Extrae el nombre base del archivo sin la ruta y la extensión
        nombre_base=$(basename "$archivo" ".fastq.gz")
        
        # Crea una subcarpeta de salida para este análisis dentro de la carpeta del output
        carpeta_salida="$output_dir/${nombre_base}_fastqc_output"
        mkdir -p "$carpeta_salida"

        echo "Ejecutando FastQC en $archivo"
        fastqc "$archivo" -o "$carpeta_salida"
    else
        echo "No se encontró archivo R1 en $directorio"
    fi
done
