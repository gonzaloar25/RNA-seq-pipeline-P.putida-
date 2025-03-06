#!/bin/bash

# Directorio donde están las carpetas
directorio="/media/gonzalo/EXTERNAL_USB/rnaseq/datos/raw_data/"

# Directorio donde redirigimos la salida del comando
output_dir="/media/gonzalo/EXTERNAL_USB/rnaseq/datos/processed/2_Trimming"

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
    archivo_R1=$(find "$carpeta" -type f -name "*_R1*.fastq.gz")
    # Encuentra el archivo R2 en la carpeta (normalmente contiene _R2_ en su nombre)
    archivo_R2=$(find "$carpeta" -type f -name "*_R2*.fastq.gz") 
    
    # Ejecuta trim_galore sobre los archivos encontrados
    if [ -n "$archivo_R1" ] && [ -n "$archivo_R2" ]; then
     # Extrae el nombre base del archivo sin la ruta y la extensión
        nombre_base_R1=$(basename "$archivo_R1" "_R1_001.fastq.gz")
        
        # Crea una subcarpeta de salida para este análisis dentro de la carpeta del output, si ya existe la omite
        carpeta_salida="$output_dir/${nombre_base_R1}_trimmed"

        if [ -d "$carpeta_salida" ]; then
            echo "Carpeta de salida ya existe para $nombre_base_R1. Omitiendo."
            continue
        fi

        mkdir -p "$carpeta_salida"

        echo "Ejecutando TrimGalore! en $archivo_R1 y $archivo_R2"
        trim_galore --paired --quality 25 --length 35 --clip_R2 2 "$archivo_R1" "$archivo_R2" -o "$carpeta_salida"
    else
        echo "No se puedo ejecutar el comando: falta archivo R1 o R2 en $carpeta"
    fi
done
