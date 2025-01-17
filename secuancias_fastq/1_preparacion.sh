#!/bin/bash

# Este scrip  realiza demultiplexeo

# Creamos los param-files. Nota: los archivos param  deberan ser editados para especificar el path donde estan los .fastq de las muestras 
#de cada placas y los barcodes.

#ipyrad -n GBS_PLACA1
#ipyrad -n GBS_PLACA2
#ipyrad -n GBS_PLACA3
#ipyrad -n GBS_PLACA4

# significado de los argmentos 
# -n crea nuevos archivo param 

# Realizaremos el demultiplexeo
#ipyrad -p params-GBS_PLACA1.txt -s 1 -f
#ipyrad -p params-GBS_PLACA2.txt -s 1 -f 
#ipyrad -p params-GBS_PLACA3.txt -s 1 -f
ipyrad -p params-GBS_PLACA4.txt -s 1 -f
 
# significado de los argmentos 
# -s selecionar el paso de la pipeline de ipyrad que se quiere realizar
# -f sobrescribir los datos


#process_radtags -p /mnt/e/P1_P2/P1/ -P -b /mnt/e/P1_P2/P1/barcodes_P1.txt -o /mnt/e/sorted_fastq/ --renz-2 pstI mspI -q -r -D 

# significado de los argmentos 
# -p ruta a un directorio con las secuancias 
# -P los archivos del directorios son de doble digestion
# -b ruta a los barcodes
# -o ruta a un directorio de salida 
# -renz-2 se usa para especificar las dos enzimas
# -threads 
# -q descartar lecturas con baja calidad 
# -r códigos de barras de rescate y etiquetas RAD
# -D guarda las lecturas descartadas de baja calidad en un archivo
# -t cortar el final de las lecturas a valor que se especifica 




