#!/bin/bash

# Este scrip  realiza análisis de calidad con fastqc, los análisis se hacen para pocas muestras y de forma aleatoria

fastqc /mnt/e/P1_P2/GBS_PLACA1_fastqs/A131_N10P46_R1_.fastq.gz -o /mnt/e/P1_P2/P1/P1_FASTQC/
fastqc /mnt/e/P1_P2/GBS_PLACA1_fastqs/A131_N10P46_R2_.fastq.gz -o /mnt/e/P1_P2/P1/P1_FASTQC/

fastqc /mnt/e/P1_P2/GBS_PLACA1_fastqs/B105_N5P34_R1_.fastq.gz -o /mnt/e/P1_P2/P1/P1_FASTQC/
fastqc /mnt/e/P1_P2/GBS_PLACA1_fastqs/B105_N5P34_R2_.fastq.gz -o /mnt/e/P1_P2/P1/P1_FASTQC/

fastqc /mnt/e/P1_P2/GBS_PLACA1_fastqs/D168_N19P47_R1_.fastq.gz -o /mnt/e/P1_P2/P1/P1_FASTQC/
fastqc /mnt/e/P1_P2/GBS_PLACA1_fastqs/D168_N19P47_R2_.fastq.gz -o /mnt/e/P1_P2/P1/P1_FASTQC/

fastqc /mnt/e/P1_P2/GBS_PLACA1_fastqs/F66_N8P13_R1_.fastq.gz -o /mnt/e/P1_P2/P1/P1_FASTQC/
fastqc /mnt/e/P1_P2/GBS_PLACA1_fastqs/F66_N8P13_R2_.fastq.gz -o /mnt/e/P1_P2/P1/P1_FASTQC/

fastqc /mnt/e/P1_P2/GBS_PLACA1_fastqs/H22_N6P6_R1_.fastq.gz -o /mnt/e/P1_P2/P1/P1_FASTQC/
fastqc /mnt/e/P1_P2/GBS_PLACA1_fastqs/H22_N6P6_R2_.fastq.gz -o /mnt/e/P1_P2/P1/P1_FASTQC/