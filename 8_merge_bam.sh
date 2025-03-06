#!/bin/bash

# Directorio donde están las carpetas de entrada
directorio="/media/gonzalo/EXTERNAL_USB/tfm/datos/processed/3_Alignment"

# Directorio donde se guardarán los archivos fusionados
output_dir="/media/gonzalo/EXTERNAL_USB/tfm/datos/processed/4_Merged_BAM"

# Array para almacenar los sorted.bam
archivos_BAM=()

# Verifica si el directorio de salida existe, si no, lo crea
if [ ! -d "$output_dir" ]; then
    echo "Creando directorio de salida: $output_dir"
    mkdir -p "$output_dir"
else
    echo "El directorio de salida ya existe: $output_dir"
fi

# Recorre todas las subcarpetas
for carpeta in "$directorio"/*/; do
    # Encuentra el archivo sorted.bam en la carpeta
    archivo_BAM=$(find "$carpeta" -type f -name "*.sorted.bam")
    
    if [ -n "$archivo_BAM" ]; then
        archivos_BAM+=("$archivo_BAM")  
    else
        echo "No se encontró archivo .sorted.bam en $carpeta"
    fi
done

# Fusiona los archivos BAM en grupos de 4 o 2 si el nombre contiene 'afleQ'
i=0
while [ $i -lt ${#archivos_BAM[@]} ]; do
    # Obtiene el nombre base del primer archivo
    nombre_base=$(basename "${archivos_BAM[i]}" | sed -E 's/(.*)RNA.*$/\1/; s/_S[0-9]*_L[0-9]*\.sorted\.bam$//')

    # Verifica si el nombre de la carpeta contiene 'AfleQ'
    if [[ "$nombre_base" == AfleQ* ]]; then
        step=2
    else
        step=4
    fi

    # Crea una subcarpeta para los archivos fusionados
    carpeta_merged="$output_dir/$nombre_base"

    # Verifica si la carpeta ya existe
    if [ -d "$carpeta_merged" ]; then
        echo "La carpeta $carpeta_merged ya existe. Omitiendo la fusión para este conjunto."
        i=$((i + step))
        continue
    fi

    mkdir -p "$carpeta_merged"

    # Define el archivo de salida para la fusión
    archivo_merged="$carpeta_merged/${nombre_base}_merged.bam"

    # Fusiona archivos BAM en grupos de 2 o 4, dependiendo del prefijo del nombre
    echo "Fusionando archivos BAM: ${archivos_BAM[@]:i:step}"
    samtools merge -@ 9 "$archivo_merged" "${archivos_BAM[@]:i:step}"
    samtools sort -@ 9 "$archivo_merged" 

    # Indexar el archivo BAM resultante
    echo "Indexando el archivo BAM fusionado: $archivo_merged"
    samtools index -@ 9 "$archivo_merged"

    # Incrementa el índice en función de si se usó 2 o 4 archivos
    i=$((i + step))
done

