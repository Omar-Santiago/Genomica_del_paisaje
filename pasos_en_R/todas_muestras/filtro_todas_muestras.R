library(usethis)
library(devtools)
library("qvalue")
library(SNPfiltR) # filtros
library(vcfR)
library(fsthet) # prueba de outliers
library(dartR) # eliminar loci desde R
library(ade4)
library(adegenet)
library(lattice)

# citas y scritps
# DeRaad, D.A. (2022), SNPfiltR: an R package for interactive and 
# reproducible SNP filtering. Molecular Ecology Resources, 22, 
# 2443–2453, 1–15. https://doi.org/10.1111/1755-0998.13618.

# Flanagan, S. P., & Jones, A. G. (2017). Constraints on the 
# FST–heterozygosity outlier approach. Journal of Heredity, 
# 108(5), 561-573.

# O'Leary, S. J., Puritz, J. B., Willis, S. C., Hollenbeck, C. M., 
# & Portnoy, D. S. (2018). These aren’t the loci you’e looking for: 
# Principles of effective SNP filtering for molecular ecologists.
#  Mol Ecol. 2018; 27: 1–14. https://doi.org/10.1111/mec.14792

# https://devonderaad.github.io/SNPfiltR/ 
# https://github.com/spflanagan/fsthet_analysis 

rm(list = ls())

########## filtros #########

vcfR_0 <- read.vcfR("F:/salida_cluster/populations.snps.vcf")
vcfR_0

# eliminamos los loci con una profundidad menor a 10
#hard_filter(vcfR=vcfR_0)
vcfR_1<-hard_filter(vcfR=vcfR_0, depth = 10)
vcfR_1

# aplicamos un filtro del balance de alelos
vcfR_2<-filter_allele_balance(vcfR_1)
vcfR_2

# eliminar loci con profundidad mayor a 100
max_depth(vcfR_2)
vcfR_3<-max_depth(vcfR_2, maxdepth = 100)
max_depth(vcfR_3)

# eliminamos individuos con missing data mayor al 90%
popmap <- read.csv("popmap_placa.csv", header = T)
missing_by_sample(vcfR=vcfR_3, popmap = popmap)
vcfR_4<-missing_by_sample(vcfR=vcfR_3, cutoff = .9)
popmap<-popmap[popmap$id %in% colnames(vcfR_4@gt),]
vcfR_4

# eliminamos singletons
vcfR_5<-min_mac(vcfR_4, min.mac = 2)
vcfR_5

#eliminamos locus con datos faltantes 
# analisis exploratorio
#miss<-assess_missing_data_pca(vcfR=vcfR_5, popmap = popmap, thresholds = c(.5,.6,.7,.8), clustering = FALSE)
vcfR_miss0.75<-missing_by_snp(vcfR_5, cutoff = .75)
vcfR_miss0.75

# eliminamos locus en desequilibrio de ligamiento
pre_vcf <- distance_thin(vcfR_miss0.75, min.distance = 100)
pre_vcf
vcfR::write.vcf(pre_vcf, "pre_vcf_todos.vcf")

########### eliminar sesgo placa ##############
# convetimosa objeto genlight
placa_genl <- gl.read.vcf("pre_vcf_todos.vcf")
pop(placa_genl) <- popmap$pop

# filtramos los loci con valores altos de hw
hw_placa <- gl.filter.hwe(placa_genl, subset = "each", n.pop.threshold = 1,
                          method_sig = "Exact", multi_comp = FALSE, multi_comp_method = "BY", 
                          alpha_val = 0.05, pvalue_type = "midp", min_sample_size = 5,
                          verbose = 3)

# guardamos los locus en un archivo de texto y los eliminamos en vcftools
all_loci <- as.data.frame(vcfR_0@fix)
loc.names <- hw_placa@loc.names
lista_reducida_1 <- as.data.frame(loc.names)
lista_eliminar_1<-all_loci[!all_loci$ID %in% lista_reducida$loc.names,]
write.csv(lista_eliminar_1$ID, "loci_eliminados_hw.csv")

# creamos una lista con los snps eliminados:
all_names <- as.data.frame(colnames(vcfR_0@gt))
pop_placa <- read.csv("popmap_placa.csv", header = T)
popmap<-pop_placa[pop_placa$id %in% colnames(pre_vcf@gt),]
ind_eliminar<-all_names[!all_names$`colnames(vcfR_0@gt)` %in% popmap$id,]
write.table(ind_eliminar, "todas_placa_ind_eliminados.txt")

# este paso elimina el individuo PP2_NNAP1
# vcftools --vcf /mnt/f/salida_cluster/populations.snps.vcf --remove todas_placa_ind_eliminados.txt --recode --recode-INFO-all --out bosque_limpios

# vcftools --vcf bosque_limpios.recode.vcf --exclude loci_eliminados_hw.txt --recode --recode-INFO-all --out vcfR_p

# creamos una matriz con formato similar a genepop
vcfR_p <- read.vcfR("vcfR_p.recode.vcf")
geno <- extract.gt(vcfR_p) # Character matrix containing the genotypes
geno <- t(geno)
G <- matrix(geno, nrow = nrow(geno), ncol = ncol(geno))

G[geno %in% c("0/0", "0|0", NA)] <- "0101"
G[geno  %in% c("0/1", "0|1")] <- "0102"
G[geno  %in% c("1/0", "1|0")] <- "0201"
G[geno %in% c("1/1", "1|1")] <- "0202"

colnames(G) <- vcfR_p@fix[,3]
pop.info <- popmap$pop
ind.names <- popmap$id
G <- cbind(pop.info, ind.names, G)

# calculamos los Fst y Ht para los loci (Flanagan & Jones, (2017)
fsts<-calc.actual.fst(G,"fst")
head(fsts)

par(mar=c(4,4,1,1))
plot(fsts$Ht, fsts$Fst,xlab="Ht",ylab="Fst",pch=19)

# estimamos los outlayers y los guardamos para eliminarlos 
quant.out<-fst.boot(G, bootstrap = F)
outliers<-find.outliers(fsts,boot.out=quant.out)
snps_eliminar <- as.vector(outliers$Locus)
write.csv(snps_eliminar, "outliers_efecto_placa_todos.csv")
head(outliers)
out.dat<-fsthet(G)
head(out.dat)

# este paso elimina los loci atipicos 
# vcftools --vcf vcfR_p.recode.vcf --exclude outliers_efecto_placa_todos.txt --recode --recode-INFO-all --out todos_snps_limpios

# creamos un genlight y eliminamos los SNPs con fst altos que causan 
# el efecto placa. Es necesario instalar dart
vcfR_limpio <- read.vcfR("todos_snps_limpios.recode.vcf")
gen_limpio <- gl.read.vcf("todos_snps_limpios.recode.vcf")
pop(gen_limpio) <- popmap$pop

# Realizamos un DAPC para ver el efecto de la placa 
pra_limpio<-xvalDapc(tab(gen_limpio,NA.method="mean"), pop(gen_limpio))
plotdapc <- scatter(pra_limpio$DAPC, mstree = T, clabel = T, lwd = 2, grid=T, cex=3,
                    col = colores, cleg = 0.80, posi.da = "bottomleft" )



############## error inerplaca e intrapalca  ############
#leemos los archivos vcf que nos devuelve el modulo populations de stacks 
# este paso lee los nombres de los archivos vcf
temp <- list.files(pattern="*_snps_limpios.recode.vcf")
# leemos los vcf, y a cada elemento de la lista agregamos el nombre antes extraido
vcf.files <- lapply(setNames(temp, make.names(gsub("*.snps.vcf$", "", temp))), read.vcfR)
# convertir los objetos vcf a genlight
genlight.objects <- lapply(vcf.files, vcfR2genlight)

# llamamos las funciones previamene obtenidas de Giles-Perez et al., 2022
source("F:/pasos_en_R/ensamble_parametros/Distancia_euclideana.r")
source("F:/pasos_en_R/ensamble_parametros/snp_errorT.r") #Alicia Mastretta-Yanes (https://github.com/AliciaMstt/RAD-error-rates).
source("F:/pasos_en_R/ensamble_parametros/LociAllele_errorT.r")

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

################ dapc ########

# muestras paisaje
# vcftools --vcf todos_snps_limpios.recode.vcf --keep muestras_bosque.txt --recode --recode-INFO-all --out muestras_paisaje

# mustras camaras ozono
# vcftools --vcf todos_snps_limpios.recode.vcf --remove muestras_bosque.txt --recode --recode-INFO-all --out muestras_camaras

# cargamos los dos archivos y los transformamos a genlight
gen_paisaje<- gl.read.vcf("muestras_paisaje.recode.vcf")
gen_camaras <- gl.read.vcf("muestras_camaras.recode.vcf")

# introducimos un popmap con las poblaciones artificiales
pop_placa <- read.csv("popmap_placa.csv", header = T)
popmap_p<-pop_placa[pop_placa$id %in% gen_paisaje@ind.names,]
popmap_c<-pop_placa[pop_placa$id %in% gen_camaras@ind.names,]
pop(gen_paisaje) <- popmap_p$pop
pop(gen_camaras) <- popmap_c$pop

# Realizamos un DAPC de las muestras del paisaje y despues de las camaras 
dapc_paisaje<-xvalDapc(tab(gen_paisaje,NA.method="mean"), pop(gen_paisaje))
dapc_camaras<-xvalDapc(tab(gen_camaras,NA.method="mean"), pop(gen_camaras))

plotdapc <- scatter(dapc_paisaje$DAPC, mstree = T, clabel = T, lwd = 2, grid=T, cex=3,
                    col = colores, cleg = 0.80, posi.da = "bottomleft" )

plotdapc <- scatter(dapc_camaras$DAPC, mstree = T, clabel = T, lwd = 2, grid=T, cex=3,
                    col = colores, cleg = 0.80, posi.da = "bottomleft" )
