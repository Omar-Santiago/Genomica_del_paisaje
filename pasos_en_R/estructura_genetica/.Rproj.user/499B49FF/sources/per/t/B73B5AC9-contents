######## Paquetes ###########
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("LEA")
library(LEA)
library(vcfR)

setwd("F:/pasos_en_R/estructura_genetica/")


############# snmf ##########
# Se requiere un archivo en formato genepop, asi que primero usaremos PGDSpider para convertir el vcf a ped para 
# despues usar la funcion ped2geno

ind_fil_snmf <- ped2geno("filtrado.ped")
snmf_estructura = snmf("filtrado.geno", 
                       iterations=100, K=1:10, rep=30, 
                       entropy=T, CPU=4, ploidy=2, 
                       project="new")
plot(snmf_estructura,lwd=8,col="#9215D0",pch=1, cex=2)
snmf_ind_estructura = load.snmfProject("filtrado.snmfProject")
best=which.min(cross.entropy(snmf_ind_estructura, K=2))
my.colors <- c("tomato", "lightblue",
               "olivedrab", "gold")
barchart(snmf_ind_estructura, K = 2, run = best,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuos",
         ylab = "Proporciones ancestrales",
         main = "Ancestria") -> bp
axis(1, at = 1:length(bp$order),
     labels = bp$order, las=1,
     cex.axis = .4)

# ahora que estimamos la estructura neutra vamos a imputar el LFMM para despues hacer un pca
ind_pca <- ped2lfmm("filtrado.ped")
impute(snmf_estructura, "filtrado.lfmm",  method = 'mode', K = 1, run = best)

# ahora corremos un pca 
pc <- pca("filtrado.lfmm_imputed.lfmm")
project = load.pcaProject("filtrado.lfmm_imputed.pcaProject") 

# con esto estimaremos la estructura gentica neutra 
estructura <- read.table("filtrado.lfmm_imputed.pca/filtrado.lfmm_imputed.projections")
estructura1 <- estructura[1:4]
write.csv(estructura1, "Structura4PCAs.csv")

plot(estructura1$V1, estructura1$V2)






