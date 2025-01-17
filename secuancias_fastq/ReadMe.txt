# Explicacion de las carpetas y scripts
# secuancias_fastq = carpeta donde estan todos los scripts para realizar un ensamble de novo:
 
# 1_preparacion.sh realiza el demultiplexeo, primero correr para generar los param-files, 
# despues, recuerda ocultar las líneas de los param-files para volver a correr el script
# 2_preparacion.sh genera análisis de calidad con fastqc
# 3_preparacion.sh realiza filtros de calidas extra

# 1_ensamble_*.sh realiza un ensamble con stacks con una muestra de 9 indiviuos con repeticiones interplaca e intraplaca, las salidas de estos scrips se editaron 
# y se movieron a la carpeta pasos_en_R/ensamble_parametros
 
# el ensamble de todas las placas se ralizo con el script: 1_ensamble_m5_M2.sh, el cual se ejecuto en un cluster y los datos de salida fueron depositados en 
# en la carpeta salida_cluster, donde se separaron las muestras del bosque y las usadas en las camaras de ozono, el siguiente paso se realizara en la carpeta:
# pasos_en_R/filtros_post_ensamble

# paquetes o software que se usaron 
# ipyrad 0.9.102
# stacks 
# FastQC v0.12.1 
# fastp 0.24.0
# PLINK v1.90b6.21 

intalacion de software 
# conda install bioconda::ipyrad 
# conda install bioconda::stacks 
# conda install bioconda::fastqc 
# conda install bioconda::plink