#!/bin/bash

# Directorio donde están las carpetas
directorio="/media/gonzalo/EXTERNAL_USB/rnaseq/datos/processed/2_Trimming"

# Directorio donde redirigimos la salida del comando
output_dir="/media/gonzalo/EXTERNAL_USB/rnaseq/datos/processed/2_Trimming/FastQC"

# Verifica si el directorio de salida existe, si no, lo crea
if [ ! -d "$output_dir" ]; then
    echo "Creando directorio de salida: $output_dir"
    mkdir -p "$output_dir"
else
    echo "El directorio de salida ya existe: $output_dir"
fi

# Recorre todas las subcarpetas
for carpeta in "$directorio"/*/; do
    # Encuentra el archivo R1 en la carpeta (normalmente contiene _R1_ en su nombre)
    archivo_R1=$(find "$carpeta" -type f -name "*_R1*.fq.gz")
    
    # Ejecuta FastQC sobre el archivo R1 encontrado
    if [ -n "$archivo_R1" ]; then
     # Extrae el nombre base del archivo sin la ruta y la extensión
        nombre_base=$(basename "$archivo_R1" ".fastq.gz")
        
        # Crea una subcarpeta de salida para este análisis dentro de la carpeta del output
        carpeta_salida="$output_dir/${nombre_base}_fastqc_output"
        mkdir -p "$carpeta_salida"

        echo "Ejecutando FastQC en $archivo_R1"
        fastqc "$archivo_R1" -o "$carpeta_salida"
    else
        echo "No se encontró archivo R1 en $carpeta"
    fi
done
