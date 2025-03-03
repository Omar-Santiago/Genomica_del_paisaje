#!/bin/bash

# este scrip genera reportes qc para muestras de las 4 placas

src=/mnt/f/trim

#files=$(echo -e "A10_N07P06\nB113_N05P34\nD70_N09P19\nE89_N16P29\nK118_NNAP45\nA122_N03P40\nB84_NNAP24\nC36_N06P12\nD163_N18P43\nE166NNAP38\nA131_N10P46\nC17_N5P6\nD168_N19P47\nE66_N14P20\nF100_N10P32\nA27_N8P9\nB103_N11P34\nD142_N14P38\nE47_N11P11\nF108_N8P32")

files=$(echo -e "B103r_N11P34\nB103_N11P34\nK11r_N4P1\nK11_N4P1\nPP4r_NNAP1\nPP4_NNAP1\nM2r_N5P9\nM2_N5P9\nAMP08r_NNAPNA\nAMP05r_NNAPNA\nAMP05_NNAPNA\nCA10_N04TAF\nF101_N9P32\nAMP08_NNAPNA\nCC26r_N04TAF\nCC26_N04TAF\nCA10r_N04TAF\nF101r_N9P32")

for samples in $files; do 

fastqc $src/${samples}.1.fastq.gz -o $src/reportes_placas
fastqc $src/${samples}.2.fastq.gz -o $src/reportes_placas

done


