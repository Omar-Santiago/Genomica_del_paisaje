#!/bin/bash

#SBATCH --mem 64000
#SBATCH -n 8

src=~/trim

files=$(echo -e "B103r_N11P34\nB103_N11P34\nK11r_N4P1\nK11_N4P1\nPP4r_NNAP1\nPP4_NNAP1\nM2r_N5P9\nM2_N5P9\nAMP08r_NNAPNA\nAMP05r_NNAPNA\nAMP05_NNAPNA\nCA10_N04TAF\nF101_N9P32\nAMP08_NNAPNA\nCC26r_N04TAF\nCC26_N04TAF\nCA10r_N04TAF\nF101r_N9P32")

id=1

for sample in $files; do
ustacks -f $src/${sample}.1.fastq.gz -o $src/parametros/m3 --name $sample -m3 -M 1 -t 4 -i $id --force-diff-len
((id++))
done

cstacks -n 1 -P $src/parametros/m3 -M /LUSTRE/Genetica/oclemente/secuancias_fastq/popmap_abies_2.txt 

sstacks -P $src/parametros/m3 -M /LUSTRE/Genetica/oclemente/secuancias_fastq/popmap_abies_2.txt 

tsv2bam -P $src/parametros/m3 -M /LUSTRE/Genetica/oclemente/secuancias_fastq/popmap_abies_2.txt  --pe-reads-dir $src -t 4

gstacks -P $src/parametros/m3 -M /LUSTRE/Genetica/oclemente/secuancias_fastq/popmap_abies_2.txt 

populations -P $src/parametros/m3 -M /LUSTRE/Genetica/oclemente/secuancias_fastq/popmap_abies_2.txt  -r 0.65 --min-maf 0.05 --plink --genepop --structure --vcf 

#plink --file $src/stacks/m3/populations.plink --cluster --mds-plot 3 --allow-extra-chr --out $src/stacks/m3/mds_m3

#plink --file $src/stacks/m3/populations.plink --cluster --pca 3 --allow-extra-chr --out $src/stacks/m3/pca_m3

