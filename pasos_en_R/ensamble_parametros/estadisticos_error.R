#
# script tomado de Mastretta-Yanes et al., 2015 y Giles-Perez et al., 2022
# agregar citas y guthub*****
#
# Este script estima las distancias euclidianas y tasas de error
# 

library(ade4)
library(adegenet)
library(vcfR)
library(lattice)


# leemos los archivos vcf que nos devuelve el modulo populations de stacks 
# este paso lee los nombres de los archivos vcf
temp <- list.files(pattern="*.vcf")
# leemos los vcf, y a cada elemento de la lista agregamos el nombre antes extraido
vcf.files <- lapply(setNames(temp, make.names(gsub("*.vcf$", "", temp))), read.vcfR)
# convertir los objetos vcf a genlight
genlight.objects <- lapply(vcf.files, vcfR2genlight)

# llamamos las funciones previamene obtenidas de Giles-Perez et al., 2022
source("Distancia_euclideana.r")
source("snp_errorT.r") #Alicia Mastretta-Yanes (https://github.com/AliciaMstt/RAD-error-rates).
source("LociAllele_errorT.r")

# extraemos la tasa de errores de los snps
ErrorRates <- vector(mode = "list", length = length(genlight.objects))
for (i in 1:length(genlight.objects)){
  ErrorRates[[i]]<-SNP_error(genlight.objects[[i]], names(genlight.objects)[i])
}
write.table((do.call(rbind, ErrorRates)),"SNPerrorRates.txt", sep ="\t", quote = F, row.names = F)
spn<-read.table("SNPerrorRates.txt", sep ="\t", header = T)


# perform the same for Euclidean distances
EuclDist <- vector(mode = "list", length = length(genlight.objects))
for (i in 1:length(genlight.objects)){
  EuclDist[[i]]<-Distancia_euclideana(genlight.objects[[i]], names(genlight.objects)[i])
}
write.table((do.call(rbind, EuclDist)),"EuclDist.txt", sep ="\t", quote = F, row.names = F)
eucli<-read.table("EuclDist.txt", sep ="\t", header = T)

# and finally for allele error rate
LocAllError <- vector(mode = "list", length = length(genlight.objects))
for (i in 1:length(genlight.objects)){
  LocAllError[[i]]<-LociAllele_error(genlight.objects[[i]], names(genlight.objects)[i])
}
write.table((do.call(rbind, LocAllError)),"LocAllError.txt", sep ="\t", quote = F, row.names = F)
allel<-read.table("LocAllError.txt", sep ="\t",header = T)

### vamos hacer unas graficas



















