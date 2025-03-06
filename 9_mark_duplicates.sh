#!/bin/bash

# Directorio donde están las carpetas de entrada
directorio="/media/gonzalo/EXTERNAL_USB/rnaseq/datos/processed/3_Alignment"

# Directorio donde se guardarán los archivos fusionados
output_dir="/media/gonzalo/EXTERNAL_USB/rnaseq/datos/processed/4_Mark_duplicates"

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
    archivo_BAM=$(find "$carpeta" -type f -name "*sorted.bam")

    if [ -n "$archivo_BAM" ]; then  
        # Extrae el nombre base del archivo sin la ruta y la extensión
        nombre_base=$(basename "$archivo_BAM" "_sorted.bam")

        # Crea una subcarpeta de salida para este análisis dentro de la carpeta del output, si ya existe la omite
        carpeta_salida="$output_dir/${nombre_base}_deduplicated"

        if [ -d "$carpeta_salida" ]; then
            echo "Carpeta de salida ya existe para $nombre_base. Omitiendo creación de directorio."
        else
            mkdir -p "$carpeta_salida"
        fi

        # Archivos de salida: BAM deduplicado y métricas
        archivo_salida_bam="$carpeta_salida/${nombre_base}_deduplicated.bam"
        archivo_salida_metrics="$carpeta_salida/${nombre_base}_dup_metrics.txt"

        # Ejecuta MarkDuplicate si aún no se había hecho anteriormente
        if [ -f "$archivo_salida_bam" ]; then
            echo "Archivo de salida ya existe para $nombre_base. Omitiendo marcaje de duplicados."
        else
            echo "Ejecutando MarkDuplicates en el archivo $archivo_BAM"
            picard MarkDuplicates I="$archivo_BAM" O="$archivo_salida_bam" M="$archivo_salida_metrics"
        fi

        # Indexación del archivo BAM deduplicado
        echo "Indexando el archivo $archivo_salida_bam con samtools"
        samtools index -@ 8 "$archivo_salida_bam"

    else
        echo "No se puede ejecutar el comando: falta archivo BAM en $carpeta"
    fi
done


