#!/bin/bash

# Directorio donde están los resultados de FastQC
directorio_fastqc="/media/gonzalo/EXTERNAL_USB/rnaseq/datos/processed/1_Quality/FastQC"

# Directorio donde guardaremos los informes de MultiQC (tenemos que crearlo previamente o usar mkdir)
directorio_multiqc="/media/gonzalo/EXTERNAL_USB/rnaseq/datos/processed/1_Quality/MultiQC"

# Verifica si el directorio de salida existe, si no, lo crea
if [ ! -d "$output_dir" ]; then
    echo "Creando directorio de salida: $output_dir"
    mkdir -p "$output_dir"
else
    echo "El directorio de salida ya existe: $output_dir"
fi

# Contador para los lotes
contador_otros=0

# Arrays para almacenar los directorios de los resultados
declare -a carpetas_otros=()

# Recorre todas las carpetas en el directorio de resultados
for carpeta in "$directorio_fastqc"/*/; do
        # Agrega las demás carpetas al array
        carpetas_otros+=("$carpeta")
        contador_otros=$((contador_otros + 1))

        # Si hemos acumulado 4 carpetas de otros, ejecutamos MultiQC
        if (( contador_otros % 4 == 0 )); then
            # Extrae el nombre de la muestra del primer archivo FastQC en el lote
            nombre_muestra=$(basename "${carpetas_otros[0]}" | cut -d'_' -f1)  
            informe_multiqc="$directorio_multiqc/${nombre_muestra}_multiqc_report"

            echo "Ejecutando MultiQC en ${carpetas_otros[@]}"
            multiqc "${carpetas_otros[@]}" -o "$directorio_multiqc" -n "${nombre_muestra}_multiqc_report"

            # Vacía el array y reinicia el contador
            carpetas_otros=()
        fi
    fi
done


