#!/bin/bash

#SBATCH --mem 64000
#SBATCH -n 8

src=~/trim

files=$(echo -e "B103r_N11P34\nB103_N11P34\nK11r_N4P1\nK11_N4P1\nPP4r_NNAP1\nPP4_NNAP1\nM2r_N5P9\nM2_N5P9\nAMP08r_NNAPNA\nAMP05r_NNAPNA\nAMP05_NNAPNA\nCA10_N04TAF\nF101_N9P32\nAMP08_NNAPNA\nCC26r_N04TAF\nCC26_N04TAF\nCA10r_N04TAF\nF101r_N9P32")

id=1

for sample in $files; do
ustacks -f $src/${sample}.1.fastq.gz -o $src/parametros/m5 --name $sample -m5 -M 1 -t 8 -i $id --force-diff-len
((id++))
done

cstacks -n 1 -P $src/parametros/m5 -M /LUSTRE/Genetica/oclemente/secuancias_fastq/popmap_abies_2.txt 

sstacks -P $src/parametros/m5 -M /LUSTRE/Genetica/oclemente/secuancias_fastq/popmap_abies_2.txt 

tsv2bam -P $src/parametros/m5 -M /LUSTRE/Genetica/oclemente/secuancias_fastq/popmap_abies_2.txt  --pe-reads-dir $src -t 8

gstacks -P $src/parametros/m5 -M /LUSTRE/Genetica/oclemente/secuancias_fastq/popmap_abies_2.txt 

populations -P $src/parametros/m5 -M /LUSTRE/Genetica/oclemente/secuancias_fastq/popmap_abies_2.txt  -r 0.65 --min-maf 0.05 --plink --genepop --structure --vcf 
