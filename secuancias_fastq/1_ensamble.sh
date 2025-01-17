#!/bin/bash

# Este script  realiza un ensable de novo

src=/mnt/e/P1_P2/P1/trim
#for files in $src/*.trimmed.paired_r1.fastq.gz; do
#sample=$(basename "$files" ".trimmed.paired_r1.fastq.gz")
#echo "las muestras forward son: ${sample}.trimmed.paired_r1.fq.gz"

#files="A1_N10P2 A19_N5P A38_N5P12 A131_N10P46 B103r_N11P34 B105_N5P34"

# Lista de muestras, recordar usar \n
files=$(echo -e "G89_N11P30\nG12_NNAP5\nI74_N14P29\nF104_N8P32\nI58_N12P21\nE51_N11P15\nI97_N11P34\nK99_N10P33\nH16_N9P6\nH42_N4P12")

# Inicializar un contador para los IDs únicos
id=1

#
# Build loci de novo in each sample for the single-end reads only. If paired-end reads are available, 
# they will be integrated in a later stage (tsv2bam stage).
# This loop will run ustacks on each sample, e.g.
#   ustacks -f ./samples/sample_01.1.fq.gz -o ./stacks --name sample_01 -M 4 -t 8
#
 
for sample in $files; do
#echo "Procesando muestra: $sample"
#echo "Archivo de entrada: $src/${sample}.trimmed.paired_r1.fastq.gz"
#sample=$(basename "$files" ".trimmed.paired_r1.fastq.gz")
ustacks -f $src/${sample}.1.fastq.gz -o $src/stacks --name $sample -M 4 -t 4 -i $id --force-diff-len
((id++))
done

# 
# Build the catalog of loci available in the metapopulation from the samples contained
# in the population map. To build the catalog from a subset of individuals, supply
# a separate population map only containing those samples.
#
cstacks -n 4 -P $src/stacks/ -M /mnt/e/secuancias_fastq/popmap_abies_2.txt 

#
# Run sstacks. Match all samples supplied in the population map against the catalog.
#
sstacks -P $src/stacks/ -M /mnt/e/secuancias_fastq/popmap_abies_2.txt

#
# Run tsv2bam to transpose the data so it is stored by locus, instead of by sample. We will include
# paired-end reads using tsv2bam. tsv2bam expects the paired read files to be in the samples
# directory and they should be named consistently with the single-end reads,
# e.g. sample_01.1.fq.gz and sample_01.2.fq.gz, which is how process_radtags will output them.
#
tsv2bam -P $src/stacks/ -M /mnt/e/secuancias_fastq/popmap_abies_2.txt --pe-reads-dir $src -t 4
#
# Run gstacks: build a paired-end contig from the metapopulation data (if paired-reads provided),
# align reads per sample, call variant sites in the population, genotypes in each individual.
#
gstacks -P $src/stacks/ -M /mnt/e/secuancias_fastq/popmap_abies_2.txt

#
# Run populations. Calculate Hardy-Weinberg deviation, population statistics, f-statistics and 
# smooth the statistics across the genome. Export several output files.
#
populations -P $src/stacks/ -M /mnt/e/secuancias_fastq/popmap_abies_2.txt -r 0.65 --vcf --genepop --fstats --smooth --hwe