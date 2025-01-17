#!/bin/bash

# Este script realiza pasos extra para la limpieza de datos, en este paso se usa Trimmomatic:
# cortar o eliminar las secuancias con calidad menor a 20 o 25 
# eliminar adaptadores de illumina
# elimar las secuencias de poli G (errores de illumina)
# eliminar las secuancias de 30pb o menos
# cortar los extremos 3 y 5, ademas, hacer el recorte por ventanas deslizantes
# hacer una correccion de la calidad al empalmar las dos secuancias foward y reverse


for files in /mnt/f/GBS_PLACA4_v2_fastqs/*_R1_.fastq.gz; do 

base=$(basename "$files" "_R1_.fastq.gz")
read1=/mnt/f/GBS_PLACA4_v2_fastqs/${base}_R1_.fastq.gz
read2=/mnt/f/GBS_PLACA4_v2_fastqs/${base}_R2_.fastq.gz
 
OutputForwardPaired=/mnt/f/trim/${base}.1.fastq.gz

OutputReversePaired=/mnt/f/trim/${base}.2.fastq.gz

fastp --in1 $read1 --in2 $read2 --out1 $OutputForwardPaired --out2 $OutputReversePaired --detect_adapter_for_pe --trim_poly_g -l 140 --cut_front --cut_tail --cut_right_window_size 4 --cut_right_mean_quality 15 --correction
done

# $ se usa para referirse a lo que esta dentro de una variable, similar a cuando se le da un valor a una varaibles en R y la varaible se unsa en una funcion (e.g anova <- True, y anova se usa en analis).
# en el caso del scritp $ se usa para almacenar rutas absolutas a los archivos fastq, ademas sirve para especificar los nombres de los archivos dentro de un directorio y poderlos identificar para asi analisarlos por 
# separada
