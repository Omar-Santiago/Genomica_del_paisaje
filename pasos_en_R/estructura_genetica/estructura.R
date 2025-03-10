######## Paquetes ###########
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("adegenet")
#install.packages("adegenet")
#install.packages("pcadapt")
library(usethis)
library("devtools")
#install_github("Andrew-Shirk/sGD")
# library("sGD") # para obtener medidad de diversida segun coordenadas

library(LEA) # para correr snmf
library(vcfR)
library(ade4)
library(adegenet) # para dapc, find.cluster, estimar alelos compartidos  
library(ggplot2)
library(dplyr)
library(dartR) # para

library(terra)
library(scatterpie) # para hacer mapas

library(ggplot2)
library(reshape2) # para generar graficos de la distribucion de daño

library(gdsfmt)
library(SNPRelate) # para analisis de identidad por estado

library(fsthet) # estimar loci candidatos y separarlos de los neutrales
library(pcadapt) # estimar loci candidatos y separarlos de los neutrales
library(qvalue)


library(bedr)

#
########## find.cluster + dapc #######
# creamos un genlight
vcfR_paisaje <- read.vcfR("F:/Gih/Genomica_del_paisaje/pasos_en_R/filtros_post_ensamble/muestras_paisaje.recode.vcf")
vcf_paisaje<-vcfR2genlight(vcfR_paisaje)

# corremos find.cluster, se escoje un numero maximo de clusters por si existe un efecto placa,
grp <- find.clusters(vcf_paisaje, max.n.clust = 5, n.pca = 100, choose = FALSE, stat = "BIC", method = "kmeans")
plot(grp$Kstat, type = "o", xlab = "numero de grupos (K)",
     ylab = "BIC",
     main = "find.clusters")

grp$grp

# guardamos los resultados
grupos_formados <- as.data.frame(grp$grp)
write.csv(grupos_formados, "grupos_formados_find.clusters_2.csv")

# corremos dapc con los grupos formados
dapc1 <- dapc(vcf_paisaje, grp$grp, n.pca = 50, n.da = 4) 
scatter(dapc1) # graficamos, los resultados indican que el efecto placa sigue
#
########## snmf ##########
# Se requiere un archivo en formato genepop, asi que primero usaremos PGDSpider para convertir el vcf a ped para 
# despues usar la funcion ped2geno, el archivo vcf de partida es muestras_paisaje.recode.vcf

ped2geno("paisaje.ped")
ped_data <- read.table("paisaje.ped", header = FALSE, sep = " ")
geno <- read.geno("paisaje.geno") # esto lo hacemos para ver el orden de los individuos

# corremos un snmf con pocas iteraciones
snmf_estructura = snmf("paisaje.geno", 
                       iterations=1000, K=1:10, rep=30, 
                       entropy=T, CPU=4, ploidy=2, 
                       project="new")
plot(snmf_estructura,lwd=8,col="#9215D0",pch=1, cex=2)
snmf_ind_estructura = load.snmfProject("paisaje.snmfProject")
best=which.min(cross.entropy(snmf_ind_estructura, K=1))
my.colors <- c("tomato", "lightblue",
               "olivedrab", "gold")

barchart(snmf_ind_estructura, K = 1, run = best,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuos",
         ylab = "Proporciones ancestrales",
         main = "Ancestria") -> bp
axis(1, at = 1:length(bp$order),
     labels = bp$order, las=1,
     cex.axis = .4)

########## colinearidad ########

# cargamos los datos de las coordenadas y datos geneticos 
coord <- read.csv("F:/fenotipo_ambiente/ambiente/Qgis (1)/ind_gbs_coord.csv")
puntos <- st_as_sf(coord, coords = c("longitude", "latitude"), crs = 6369)
print(puntos)

IBD_dist <- geo_dist(puntos, type = "Euclidean")

proj <- CRS("+init=epsg:6369")
xy_points <- SpatialPoints(coord[, c(2, 3)], proj4string = proj)
IBD_dist <- distmat(xy_points,method="ed")
IBD_dist [1:10,1:10]

vcfR_paisaje <- read.vcfR("F:/Gih/Genomica_del_paisaje/pasos_en_R/todas_muestras/muestras_paisaje.recode.vcf")

paisaje_gen <- vcfR2genind(vcfR_paisaje)
pop(paisaje_gen) <- coord$pop

radius <- 500 # for this demo, units are in meters 6000
IBD_GD <- sGD(paisaje_gen,
              coord,
              IBD_dist,
              radius,
              min_N=5,
              metrics="GD")

write.csv(IBD_GD, "IBD_GD_resultado.csv")

# esta aproximacion es menos eficiente y una mas actual es la propuesta por Bishop: 
# https://doi.org/10.1111/2041-210X.14090  
# pero esta aproximacion se realizara en la parte de diversidad_genetica


########## dapc categorias de saño #####
vcfR_paisaje <- read.vcfR("F:/Gih/Genomica_del_paisaje/pasos_en_R/filtros_post_ensamble/muestras_paisaje.recode.vcf")
vcf_paisaje<-vcfR2genlight(vcfR_paisaje)

# cargamos los datos de las categorias de daño
meta <- read.csv("popmap_paisaje_dapc.csv")
pop(vcf_paisaje) <- meta$ozono

# xvalDAPC realiza un 
paj<-xvalDapc(tab(vcf_paisaje,NA.method="mean"), pop(vcf_paisaje))
plotdapc <- scatter(paj$DAPC, mstree = T, clabel = T, lwd = 2, grid=T, cex=3,
                    cleg = 0.80, posi.da = "bottomleft", col = c("darkgreen","gold2", "chocolate1", "red4")  )


#
########## mapa con la distrubucion de las categorias de daño ######
# ahora vamos a generar un mapa con la distrubucion de las categorias de daño
# scrip tomado de: https://github.com/AliciaMstt/monitoreo-oyameles/blob/main/scripts/3_figuras_propuesta.html

meta <- read.csv("popmap_paisaje_dapc.csv")
meta$sano <- ifelse(meta$ozono == "0%", 1, 0)
meta$menos10 <- ifelse(meta$ozono == "10%", 1, 0)
meta$d10_40 <- ifelse(meta$ozono == "10_40%", 1, 0)
meta$d40_70 <- ifelse(meta$ozono == "40_70%", 1, 0)

#cat <- c("sano", "menos10", "d10_40", "d40_70")

ggplot() +
  geom_scatterpie(data = meta, aes(x = p_log, y = p_lat), 
                  cols = c("sano", "menos10", "d10_40", "d40_70")) +

  scale_fill_manual(values = c("darkgreen", "gold2", "chocolate1", "red4"),
                    name = "Categorias de daño") +
  coord_fixed() +  # Mantiene la proporción correcta
  labs(x = "Longitud", y = "Latitud") +
  theme_minimal()

#my_cols2<-c("gold2", "chocolate1", "orangered", "red4", "darkorchid4")
#desired_order_percentage<-c("less than 10%", "10 to 40%", "40 to 50%", "50 to 70%", "more than 70%")

imagen <- rast("paisaje.tif")
#plotRGB(imagen, r=1, g=2, b=3, stretch="lin")
#imagen_rgb <- as.raster(imagen)

imagen_norm <- imagen / max(values(imagen), na.rm = TRUE) * 255

# Convertir a data frame para ggplot2
imagen_df <- as.data.frame(imagen_norm, xy = TRUE)


# Renombrar bandas para facilidad de uso en ggplot
colnames(imagen_df) <- c("lon", "lat", "red", "green", "blue")

# plot sampled plots
p_satmap <- ggplot() +
  geom_raster(data = imagen_df, aes(x = lon, y = lat, fill = rgb(red, green, blue, maxColorValue = 255))) +
  scale_fill_identity() +  # Usar colores tal cual en la imagen
  geom_scatterpie(data = meta,
                  aes(x = longitude, y = latitude, group = ozono),
                  cols=cat,
                  pie_scale = 2, color = NA, alpha = 1) +
  ggtitle("a)") +
  scale_fill_manual(values = c("darkgreen", "gold2", "chocolate1", "red4"),
                    name = "Estado de salud") +
  theme_minimal() +
  theme(text = element_text(size = 20), legend.position = "none")

print(p_satmap)

#

########## estructura genetica espacial a resolucion fina (FSSR) ########
# vamos a evaluar la autocorrelacion espacial, usaremos el coeficiente de parentesco de
# Loiselle para ello lo implementaremos en Spagedi (terminal). En esta seccion crearemos apartir de un vcf 
# el archivo que necesita Spagedi 

vcf_data <- vcfR_paisaje <- read.vcfR("F:/Gih/Genomica_del_paisaje/pasos_en_R/filtros_post_ensamble/muestras_paisaje.recode.vcf")
geno <- extract.gt(vcf_data) # Character matrix containing the genotypes
geno <- t(geno)
G <- matrix(geno, nrow = nrow(geno), ncol = ncol(geno))

G[geno %in% c(NA)] <- "0"
G[geno %in% c("0/0", "0|0")] <- "11/11"
G[geno  %in% c("0/1", "0|1")] <- "11/22"
G[geno  %in% c("1/0", "1|0")] <- "22/11"
G[geno %in% c("1/1", "1|1")] <- "22/22"

meta3 <- read.csv("popmap_paisaje_dapc.csv")

colnames(G) <- vcf_data@fix[,3]
Lat <- meta3$latitude
Long <- meta3$longitude
ind <- meta3$ID
G <- cbind(ind,Lat,Long,G)

write.csv(G, "relatness.csv")

# corremos SpaGeDi 1.5
# las opciones utilixadas fueron estimar los siguientes coeficientes de relacion:
# kinship coefficient estimated according to J. Nason (described in Loiselle et al. 1995).
# kinship coefficient estimated according to Ritland (1996).
# relationship coefficient estimated according to Queller and Goodnight (1989).
# relationship coefficient estimated according to Lynch and Ritland (1999)
# relationship coefficient estimated according to Wang (2002)
# relationship coefficient estimated according to Li et al. (1993).
# fraternity coefficient, Lynch and Ritland (1999)
# fraternity coefficient, Wang (2002)
# distance measure described in Rousset (2000)(a by Rousset)

# se realizaron 100 permutaciones y se optuvieron las matrices de distancia genetica

# ahora vamos a crear hetmaps 

Loiselle <- read.csv("kinship_Loiselle_1995.csv", row.names = 1)
kiship<-as.matrix(Loiselle)
gl.plot.heatmap(kiship)
heatmap(kiship, col = colorRampPalette(c("blue", "white", "red"))(100), symm = TRUE)
# agrupa a los de la placa 3, no sirve

Ritlan <-  read.csv("k_Ritlan_1996.csv", row.names = 1)
K_Ritlan <- as.matrix(Ritlan)
gl.plot.heatmap(K_Ritlan)
heatmap(K_Ritlan, col = colorRampPalette(c("blue", "white", "red"))(100), symm = TRUE)
# no agrupa a los de la placa 3, sive, todos iguales

Moran <- read.csv("Moran_I_1999.csv", row.names = 1)
r_M <- as.matrix(Moran)
gl.plot.heatmap(r_M) 
heatmap(r_M, col = colorRampPalette(c("blue", "white", "red"))(100), symm = TRUE)
# no agrupa a los de la placa 3, sive, todos iguales

Queller <- read.csv("r_Queller_Goodnight_1989.csv", row.names = 1)
r_QG <- as.matrix(Queller)
gl.plot.heatmap(r_QG) # se ve interesante
heatmap(r_QG, col = colorRampPalette(c("blue", "white", "red"))(100), symm = TRUE)
# no agrupa a los de la placa 3, sive, se ve interesante

Lynch <- read.csv("r_Linch_Ritland_1999.csv", row.names = 1)
r_LR <-  as.matrix(Lynch)
gl.plot.heatmap(r_LR)
heatmap(r_LR, col = colorRampPalette(c("blue", "white", "red"))(100), symm = TRUE)
# no agrupa a los de la placa 3, sive, todos iguales

Wand <- read.csv("r_Wang_2002.csv", row.names = 1)
r_W <- as.matrix(Wand)
gl.plot.heatmap(r_W) # se ve interesante
heatmap(r_W, col = colorRampPalette(c("blue", "white", "red"))(100), symm = TRUE)
# no agrupa a los de la placa 3, sive, interesante

Li <- read.csv("r_Li_1993.csv", row.names = 1)
r_L <- as.matrix(Li)
gl.plot.heatmap(r_L)
heatmap(r_L, col = colorRampPalette(c("blue", "white", "red"))(100), symm = TRUE)
# no agrupa a los de la placa 3, sive, todos iguales

f_LR <- read.csv("f_LR_1999.csv", row.names = 1)
fLR <- as.matrix(f_LR)
gl.plot.heatmap(fLR)
heatmap(fLR, col = colorRampPalette(c("blue", "white", "red"))(100), symm = TRUE)
# no agrupa a los de la placa 3, sive, todos iguales

f_Wang <- read.csv("f_Wang_2002.csv", row.names = 1)
f_W <- as.matrix(f_Wang)
gl.plot.heatmap(f_W)
heatmap(f_W, col = colorRampPalette(c("blue", "white", "red"))(100), symm = TRUE)
# no agrupa a los de la placa 3, sive, interesante
# interesante

a <- read.csv("a_Rousset_2000.csv", row.names = 1)
a_R <- as.matrix(a)
gl.plot.heatmap(a_R)
heatmap(a_R, col = colorRampPalette(c("blue", "white", "red"))(100), symm = TRUE)
# no agrupa a los de la placa 3, sive, iguales

vcf_data <- vcfR_paisaje <- read.vcfR("F:/Gih/Genomica_del_paisaje/pasos_en_R/filtros_post_ensamble/muestras_paisaje.recode.vcf")
genidata <- vcfR2genind(vcf_data)
sd <- propShared(genidata)
gl.plot.heatmap(sd)
heatmap(sd, col = colorRampPalette(c("blue", "white", "red"))(100), symm = TRUE)
# agrupa a los de la placa 3, no sirve


f_w_2 <- read.csv("f_wang_2.csv", row.names = 1)
fw2 <- as.matrix(f_w_2)
gl.plot.heatmap(fw2)
heatmap(fw2, col = colorRampPalette(c("blue", "white", "red"))(100), symm = TRUE)


r_w_2 <- read.csv("r_wang_2.csv", row.names = 1)
rw2 <- as.matrix(r_w_2)
gl.plot.heatmap(rw2)
rw2_long <- melt(rw2)# Convertir a formato largo
ggplot(rw2_long, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-0.5, 0.5)) +
  theme_minimal()

dgeo <- read.csv("distancia_geografica.csv", row.names = 1)
distgeo <- as.matrix(dgeo)
gl.plot.heatmap(distgeo)
heatmap(distgeo, col = colorRampPalette(c("blue", "white", "red"))(100), symm = TRUE)


########## dapc paisaje y camaras ######

pyc <- gl.read.vcf("F:/Gih/Genomica_del_paisaje/pasos_en_R/todas_muestras/todos_snps_limpios.recode.vcf")
pyc

meta_pyc <- read.csv("pyc_meta.csv")

########### todos menos los AE
ae <- c("AP", "AX", "CA", "CC")
keep_AE <- meta_pyc[!meta_pyc$categoria %in% ae, ]
pyc_ae <- gl.keep.ind(pyc, keep_AE$Column1, recalc = FALSE, mono.rm = FALSE, verbose = 3)
pop(pyc_ae) <- keep_AE$categoria

ae_dapc<-xvalDapc(tab(pyc_ae,NA.method="mean"), pop(pyc_ae))
plot_dapc_ae <- scatter(ae_dapc$DAPC, mstree = T, clabel = T, lwd = 2, grid=T, cex=3,
                    cleg = 0.80, posi.da = "bottomleft", col = c("#8867A1", "darkgreen","gold2", "chocolate1", "#5F90BB"))
resul_ae <- ae_dapc$DAPC$posterior

resul_ae <- as.data.frame(ae_dapc$DAPC$ind.coord)
rgl::plot3d(resul_ae[,1:3],type="n",col=cols)
with(resul_ae,
     rgl::text3d(resul_ae[,1],
                 resul_ae[,2],
                 resul_ae[,3],
                 rownames(resul_ae),
                 cex=1))

########### todos menos los AP
ap <- c("AE", "AX", "CA", "CC")
keep_AP <- meta_pyc[!meta_pyc$categoria %in% ap, ]
pyc_ap <- gl.keep.ind(pyc, keep_AP$Column1, recalc = FALSE, mono.rm = FALSE, verbose = 3)
pop(pyc_ap) <- keep_AP$categoria

ap_dapc<-xvalDapc(tab(pyc_ap,NA.method="mean"), pop(pyc_ap))
plot_dapc_ap <- scatter(ap_dapc$DAPC, mstree = T, clabel = T, lwd = 2, grid=T, cex=3,
                        cleg = 0.80, posi.da = "bottomleft", col = c("#8867A1", "darkgreen","gold2", "chocolate1", "#5F90BB"))
resul_ap <- ap_dapc$DAPC$posterior

########### todos menos los AX
ax <- c("AE", "AP", "CA", "CC")
keep_AX <- meta_pyc[!meta_pyc$categoria %in% ax, ]
pyc_ax <- gl.keep.ind(pyc, keep_AX$Column1, recalc = FALSE, mono.rm = FALSE, verbose = 3)
pop(pyc_ax) <- keep_AX$categoria

ax_dapc<-xvalDapc(tab(pyc_ax,NA.method="mean"), pop(pyc_ax))
plot_dapc_ax <- scatter(ax_dapc$DAPC, mstree = T, clabel = T, lwd = 2, grid=T, cex=3,
                        cleg = 0.80, posi.da = "bottomleft", col = c("#8867A1", "darkgreen","gold2", "chocolate1", "#5F90BB"))
resul_ax <- ax_dapc$DAPC$posterior


########### todos menos los CC
cc <- c("AE", "AP", "AX", "CA")
keep_CC <- meta_pyc[!meta_pyc$categoria %in% cc, ]
pyc_cc <- gl.keep.ind(pyc, keep_CC$Column1, recalc = FALSE, mono.rm = FALSE, verbose = 3)
pop(pyc_cc) <- keep_CC$categoria

cc_dapc<-xvalDapc(tab(pyc_cc,NA.method="mean"), pop(pyc_cc))
plot_dapc_cc <- scatter(cc_dapc$DAPC, mstree = T, clabel = T, lwd = 2, grid=T, cex=3,
                        cleg = 0.80, posi.da = "bottomleft", col = c("#8867A1", "darkgreen","gold2", "chocolate1", "#5F90BB"))
resul_cc <- cc_dapc$DAPC$posterior

########### todos menos los CA
ca <- c("AE", "AP", "AX", "CC")
keep_CA <- meta_pyc[!meta_pyc$categoria %in% ca, ]
pyc_ca <- gl.keep.ind(pyc, keep_CA$Column1, recalc = FALSE, mono.rm = FALSE, verbose = 3)
pop(pyc_ca) <- keep_CA$categoria

ca_dapc<-xvalDapc(tab(pyc_ca,NA.method="mean"), pop(pyc_ca))
plot_dapc_ca <- scatter(ca_dapc$DAPC, mstree = T, clabel = T, lwd = 2, grid=T, cex=3,
                        cleg = 0.80, posi.da = "bottomleft", col = c("#8867A1", "darkgreen","gold2", "chocolate1", "#5F90BB"))
resul_ca <- ca_dapc$DAPC$posterior
#


###### identity by state #####

snpgdsVCF2GDS("F:/Gih/Genomica_del_paisaje/pasos_en_R/filtros_post_ensamble/muestras_paisaje.recode.vcf", "test.gds", method="biallelic.only")
snpgdsSummary("test.gds")
genofile <-  snpgdsOpen("test.gds")

sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

meta <- read.csv("popmap_paisaje_dapc.csv")

pop_code <- meta$ozono
table(pop_code)
pop_code
head(cbind(sample.id, pop_code))

# se estima y grafica el promedio de identidades de ibs
ibs <- snpgdsIBS(genofile, num.thread=2)
pop.idx <- order(pop_code)
image(ibs$ibs[pop.idx, pop.idx], col=terrain.colors(16), cex.axis = 1)
gl.plot.heatmap(ibs$ibs)
heatmap(ibs$ibs, col = colorRampPalette(c("blue", "white", "red"))(100), symm = TRUE)

# para realizar un analisis de escalamiento multidimensional en la matris nxn del ibs
loc <- cmdscale(1 - ibs$ibs, k = 2)
x <- loc[,1]; y <- loc[,2]
race <- as.factor(pop_code)
plot(x,y,col=race,xlab="", ylab="",
     main="Multidimensional scaling analysis (ibs)")
legend("bottomright", legend=levels(race), pch="0",text.col=1:nlevels(race))

# se realiza una agrupacion por ibs y se determina los grupos mediante permutacion
set.seed(100)
ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile))

rv2 <- snpgdsCutTree(ibs.hc, samp.group = as.factor(pop_code))
plot(rv2$dendrogram, leaflab="none", main = "heatmap", pch=2)
legend("bottomleft", legend = levels(race), col=1:nlevels(race), pch=19, ncol=4)

#
###### separar loci candidatos de neutrales ####

# creamos una matriz con formato similar a genepop
vcfR_paisaje <- read.vcfR("F:/Gih/Genomica_del_paisaje/pasos_en_R/filtros_post_ensamble/muestras_paisaje.recode.vcf")
geno <- extract.gt(vcfR_paisaje) # Character matrix containing the genotypes
geno <- t(geno)
G <- matrix(geno, nrow = nrow(geno), ncol = ncol(geno))

G[geno %in% c("0/0", "0|0", NA)] <- "0101"
G[geno  %in% c("0/1", "0|1")] <- "0102"
G[geno  %in% c("1/0", "1|0")] <- "0201"
G[geno %in% c("1/1", "1|1")] <- "0202"

colnames(G) <- vcfR_paisaje@fix[,3]
pop.info <- meta$ozono
ind.names <- meta$ID
G <- cbind(pop.info, ind.names, G)

# calculamos los Fst y Ht para los loci (Flanagan & Jones, (2017)
fsts<-calc.actual.fst(G,"fst")
head(fsts)

par(mar=c(4,4,1,1))
plot(fsts$Ht, fsts$Fst,xlab="Ht",ylab="Fst",pch=19)

# estimamos los outlayers y los guardamos para eliminarlos 
quant.out<-fst.boot(G, bootstrap = T)
outliers<-find.outliers(fsts,boot.out=quant.out)
snps_outliers <- as.vector(outliers$Locus)
write.csv(snps_outliers, "outliers_pop_daño.csv")
head(outliers)
out.dat<-fsthet(G)
head(out.dat)

# ahora estimaremos los outliers con pcadapt
path_to_file <- "F:/Gih/Genomica_del_paisaje/pasos_en_R/filtros_post_ensamble/muestras_paisaje.recode.vcf"
#vcfR_paisaje <- read.vcf("F:/Gih/Genomica_del_paisaje/pasos_en_R/filtros_post_ensamble/muestras_paisaje.recode.vcf")
#x.bed <- vcf2bed(vcfR_paisaje)
path_to_file <- "snps_limpios.ped"

vcfR_paisaje <- read.vcfR("F:/Gih/Genomica_del_paisaje/pasos_en_R/filtros_post_ensamble/muestras_paisaje.recode.vcf")
geno <- extract.gt(vcfR_paisaje) # Character matrix containing the genotypes
G <- matrix(geno, nrow = nrow(geno), ncol = ncol(geno))

G[geno %in% c(NA)] <- "0"
G[geno  %in% c("0/0")] <- "0"
G[geno  %in% c("1/0", "0/1")] <- "1"
G[geno %in% c("1/1", "1|1")] <- "2"

rownames(G) <- vcfR_paisaje@fix[,3]
colnames(G) <- meta$ID

filename <- read.pcadapt(G, type = "pcadapt")
filename <- read.pcadapt(path_to_file, type = "vcf")

x <- pcadapt(filename, K=50)

# usaremos el metodo de q y fdr
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1
outpcadapt <- which(qval < alpha)
outpcadapt

opcion1 <- outpcadapt
opcion_vcf <- outpcadapt



write.csv(outpcadapt, "outpcadapt.csv")

# especificamos las categorias de daño
summary(x$scores)
print(pop_code)
plot(x, option = "score", pop=pop_code)
#

##### estructura sin outliers #####
# vcftools --vcf /mnt/f/Gih/Genomica_del_paisaje/pasos_en_R/filtros_post_ensamble/muestras_paisaje.recode.vcf --exclude out_pop_daño.txt --recode --recode-INFO-all --out vcf_sin_outliers
# vcftools --vcf /mnt/f/Gih/Genomica_del_paisaje/pasos_en_R/filtros_post_ensamble/muestras_paisaje.recode.vcf --snps out_pop_daño.txt --recode --recode-INFO-all --out vcf_con_outliers

# muestras sin outliers
vcf_sin <- gl.read.vcf("vcf_sin_outliers.recode.vcf")
# cargamos los datos de las categorias de daño
meta <- read.csv("popmap_paisaje_dapc.csv")
pop(vcf_sin) <- meta$ozono
paj_sin<-xvalDapc(tab(vcf_sin,NA.method="mean"), pop(vcf_sin))
plotdapc <- scatter(paj_sin$DAPC, mstree = T, clabel = T, lwd = 2, grid=T, cex=3,
                    cleg = 0.80, posi.da = "bottomleft", col = c("darkgreen","gold2", "chocolate1", "red4")  )


# muestras solo outliers
vcf_con <- gl.read.vcf("vcf_con_outliers.recode.vcf")
pop(vcf_con) <- meta$ozono
paj_con<-xvalDapc(tab(vcf_con,NA.method="mean"), pop(vcf_con))
plotdapc <- scatter(paj_con$DAPC, mstree = T, clabel = T, lwd = 2, grid=T, cex=3,
                    cleg = 0.80, posi.da = "bottomleft", col = c("darkgreen","gold2", "chocolate1", "red4")  )
#

# muestras con outliers
snpgdsVCF2GDS("vcf_con_outliers.recode.vcf", "con_out.gds", method="biallelic.only")
snpgdsSummary("con_out.gds")
genofile <-  snpgdsOpen("con_out.gds")
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

pop_code <- meta$ozono
table(pop_code)
pop_code
head(cbind(sample.id, pop_code))

# se estima y grafica el promedio de identidades de ibs
ibs <- snpgdsIBS(genofile, num.thread=2)
pop.idx <- order(pop_code)
image(ibs$ibs[pop.idx, pop.idx], col=terrain.colors(16), cex.axis = 1)
gl.plot.heatmap(ibs$ibs)
heatmap(ibs$ibs, col = colorRampPalette(c("blue", "white", "red"))(100), symm = TRUE)


# para realizar un analisis de escalamiento multidimensional en la matris nxn del ibs
loc <- cmdscale(1 - ibs$ibs, k = 2)
x <- loc[,1]; y <- loc[,2]
race <- as.factor(pop_code)
plot(x,y,col=race,xlab="", ylab="",
     main="Multidimensional scaling analysis (ibs)")
legend("bottomright", legend=levels(race), pch="0",text.col=1:nlevels(race))

# se realiza una agrupacion por ibs y se determina los grupos mediante permutacion
set.seed(100)
ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile))
rv2 <- snpgdsCutTree(ibs.hc, samp.group = as.factor(pop_code))
plot(rv2$dendrogram, leaflab="none", main = "heatmap", pch=2)
legend("bottomleft", legend = levels(race), col=1:nlevels(race), pch=19, ncol=4)


# muestras sin outliers
snpgdsVCF2GDS("vcf_sin_outliers.recode.vcf", "sin_out.gds", method="biallelic.only")
snpgdsSummary("sin_out.gds")
genofile <-  snpgdsOpen("sin_out.gds")
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

pop_code <- meta$ozono
table(pop_code)
pop_code
head(cbind(sample.id, pop_code))

# se estima y grafica el promedio de identidades de ibs
ibs <- snpgdsIBS(genofile, num.thread=2)
pop.idx <- order(pop_code)
image(ibs$ibs[pop.idx, pop.idx], col=terrain.colors(16), cex.axis = 1)
gl.plot.heatmap(ibs$ibs)
heatmap(ibs$ibs, col = colorRampPalette(c("blue", "white", "red"))(100), symm = TRUE)
