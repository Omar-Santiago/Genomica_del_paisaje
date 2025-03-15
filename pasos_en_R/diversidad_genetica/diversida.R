##### Paquetes ####
#remotes::install_github("hemstrow/snpR")
library(snpR)
#vamos a usar algtr
library(algatr)
# Install required packages
#wingen_packages() 
library(wingen) # no se que estimadores usa
library(sp)
library(raster)
library(terra)
library(ggplot2)
library(sf)
library(vcfR)
library(ade4)
library(adegenet) # versiones mas recientes usan pegas para diversidad
library(hierfstat) # estimaciones de diversidad genetica
library(dartR)

#
##### preparar datos ######
vcfR_paisajeisaje <- "F:/Gih/Genomica_del_paisaje/pasos_en_R/todas_muestras/muestras_paisaje.recode.vcf"
genl <- gl.read.vcf(vcfR_paisajeisaje)

coord <- read.csv("F:/fenotipo_ambiente/ambiente/Qgis (1)/ind_gbs_coord.csv")
pop(genl) <- coord$ozono 

gn <- gl2gi(genl)
pai_hier <- genind2hierfstat(gn)

tascols <- c("skyblue", "#74c476", "#FDD835", "orange")
#
##### calculate allelic richness #### 
ar <- allelic.richness(pai_hier)
summary(ar$Ar)
ar <- as.data.frame(ar$Ar)
mean.ar <- colMeans(ar)
mean.ar
# Calculate a measure of variance for allelic richness, e.g. standard error. 
# Standard error = standard deviation / sqrt(n), where n is the number of genotyped loci
sd(ar$`0%`, na.rm=TRUE)/sqrt(nrow(ar) - length(which(is.na(ar$`0%`))))

# Boxplot of allelic richness per population:
#first, extend the margins of the graphing window to fit long axis labels
boxplot(ar, ylab="Allelic richness", las=2, col=tascols)

##### estadisticos de heterocigosis ####
# Calculates basic statistics for each loci (Hs, Ho, Fis etc.) dartR y hierfstat
# The function gl.report.heterozygosity reports the observed, expected, and unbiased heterozygosities 
# and Fis (inbreeding coefficient) by population or the observed heterozygosity for each individual 
hets <- gl.report.heterozygosity(genl, method="pop", plot_colors_pop = tascols)

# ahora lo realizamos a nivel individual
ind.hets <- gl.report.heterozygosity(genl, method="ind")
# 
##### pi ####



##### diversidad genetica en el espacio #####
# cargar datos geograficos
coord <- read.csv("F:/fenotipo_ambiente/ambiente/Qgis (1)/ind_gbs_coord.csv")
puntos <- st_as_sf(coord, coords = c("longitude", "latitude"), crs = "+proj=longlat")
puntos_proj <- st_transform(puntos, crs = 6369)
print(puntos_proj)

# cargamos un shape del desierto de los leones
map_dir <- "F:/fenotipo_ambiente/ambiente/Qgis (1)/Temp_suelo/mayo_2018.tif"
map_vec_sf <-  raster(map_dir)
envlayer <- projectRaster(map_vec_sf,  crs = 6369)

# raster para las ventanas deslizantes
liz_lyr <- coords_to_raster(puntos_proj, res = 90, buffer = 1, plot = TRUE)

# vemos el ancho de la venta que se usara para los calculos 
sample_count <- preview_gd(liz_lyr, puntos_proj, 
                           method = "window", wdim = 3, fact = 0)


# Visualize the sample count layer
ggplot_count(sample_count)

# realizamos el analis (especificar las medidad de diversidad genetica)
liz_vcf <- read.vcfR("F:/Gih/Genomica_del_paisaje/pasos_en_R/todas_muestras/muestras_paisaje.recode.vcf")

# realizamos las estimacion de la diversidad 
# se deben tener minimo dos individuos para estimar la diversidad, de los contrario
# se producen NAs, en nuesstro caso asi pasa.  
wgd <- window_gd(liz_vcf,
                      puntos_proj,
                 liz_lyr,
                 stat = c("allelic_richness"),
                 wdim = 3,
                 fact = 0, 
                 min_n = 1
)

# Plot map of allelic_richness
ggplot_gd(wgd$allelic_richness) + ggtitle("allelic richness")
# hacemos interpolacion
kgd <- krig_gd(wgd$allelic_richness, index = 1:2, liz_lyr, disagg_grd = 5)
summary(kgd)
ggplot_gd(kgd) + ggtitle("Kriged Ar")


# Plot map of Ho_hierfstat
ggplot_gd(wgd$Ho_hierfstat) + ggtitle("Ho_hierfstat")
# hacemos interpolacion
kgd2 <- krig_gd(wgd$Ho_hierfstat, index = 1:2, liz_lyr, disagg_grd = 5)
summary(kgd2)
ggplot_gd(kgd2) + ggtitle("Kriged Ho_hierfstat")


# Plot map of Ho_hierfstat
ggplot_gd(wgd$Hs_hierfstat) + ggtitle("Hs_hierfstat")
# hacemos interpolacion
kgd3 <- krig_gd(wgd$Hs_hierfstat, index = 1:2, liz_lyr, disagg_grd = 5)
summary(kgd3)
ggplot_gd(kgd3) + ggtitle("Kriged Hs_hierfstat")


# Plot map of Ho_hierfstat
ggplot_gd(wgd$Ht_hierfstat) + ggtitle("Ht_hierfstat")
# hacemos interpolacion
kgd4 <- krig_gd(wgd$Ht_hierfstat, index = 1:2, liz_lyr, disagg_grd = 5)
summary(kgd4)
ggplot_gd(kgd4) + ggtitle("Kriged Ht_hierfstat")





