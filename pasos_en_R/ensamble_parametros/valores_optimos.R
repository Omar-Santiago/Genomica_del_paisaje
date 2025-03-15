###### paquetes  ####
#install.packages("adegenet", dependencies = T)

library(ade4)
library(adegenet)
library(vcfR)
library(dartR)

##

###### PCA y distancia ####
# graficas para el mds con el valor de m3 en el ensamble de novo
m3 <- read.table("mds_m3.mds", header=TRUE)

colors <- c("#6B990F", "#6B990F", "#1E8E99", "#1E8E99", "#920000", "#920000", "#B66DFF",
                    "#B66DFF", "#662700", "#662700", "#264DFF", "#264DFF", "#00FF00", "#00FF00", 
            "#FF6E00", "#FF6E00", "black", "black")
m3$color <- colors

plot(m3$C1, m3$C2)
text(m3$C1, m3$C2, labels = m3$IID, col=m3$color, pos = 4)

plot(m3$C2, m3$C3)
text(m3$C2, m3$C3, labels = m3$IID, col=m3$color, pos = 4)


# graficas para el pca con el valor de m3 en el ensamble de novo
eigen_m3 <-  read.table("pca_m3.eigenvec", header = FALSE)
eigen_m3$color <- colors
plot(eigen_m3$V3, eigen_m3$V4)
text(eigen_m3$V3, eigen_m3$V4, labels = eigen_m3$V2, col=m3$color,  pos = 4)

eigen_m5 <- read.table("pca_m5.eigenvec")
eigen_m5$color <- colors
plot(eigen_m5$V3, eigen_m5$V4)
text(eigen_m5$V3, eigen_m5$V4, labels = eigen_m5$V2, col=eigen_m5$color,  pos = 4)

eigen_m7 <- read.table("pca_m7.eigenvec")
eigen_m7$color <- colors
plot(eigen_m7$V3, eigen_m7$V4)
text(eigen_m7$V3, eigen_m7$V4, labels = eigen_m7$V2, col=eigen_m7$color,  pos = 4)

eigen_m3M2 <- read.table("pca_m3M2.eigenvec")
eigen_m3M2$color <- colors
plot(eigen_m3M2$V3, eigen_m3M2$V4)
text(eigen_m3M2$V3, eigen_m3M2$V4, labels = eigen_m3M2$V2, col=eigen_m3M2$color,  pos = 4)

eigen_m5M2 <- read.table("pca_m5M2.eigenvec")
eigen_m5M2$color <- colors
plot(eigen_m5M2$V3, eigen_m5M2$V4)
text(eigen_m5M2$V3, eigen_m5M2$V4, labels = eigen_m5M2$V2, col=eigen_m5M2$color,  pos = 4)

# graficas para el mds con el valor de m3 en el ensamble de novo
m3 <- read.table("mds_m3.mds", header=TRUE)
m3$color <- colors
plot(m3$C1, m3$C2)
text(m3$C1, m3$C2, labels = m3$IID, col=m3$color, pos = 4)

m5 <- read.table("mds_m5.mds", header=TRUE)
m5$color <- colors
plot(m5$C1, m5$C2)
text(m5$C1, m5$C2, labels = m5$IID, pos = 4, col = m5$color )

# graficas para el mds con el valor de m7 en el ensamble de novo
m7 <- read.table("mds_m7.mds", header=TRUE)
m7$color <- colors
plot(m7$C1, m7$C2)
text(m7$C1, m7$C2, labels = m7$IID, pos = 4, col = m7$color )

m3M2 <- read.table("mds_m3M2.mds", header=TRUE)
m3M2$color <- colors
plot(m3M2$C1, m3M2$C2)
text(m3M2$C1, m3M2$C2, labels = m3M2$IID, col=m3M2$color, pos = 4)

m5M2 <- read.table("mds_m5M2.mds", header=TRUE)
m5M2$color <- colors
plot(m5M2$C1, m5M2$C2)
text(m5M2$C1, m5M2$C2, labels = m5M2$IID, pos = 4, col = m5M2$color )
##
##### find.cluster ##### 
read_m3 <-  gl.read.vcf("m3.snps.vcf")
read_m5 <-  gl.read.vcf("m5.snps.vcf")
read_m7 <-  gl.read.vcf("m7.snps.vcf")


# corremos find.cluster, se escoje un numero maximo de clusters por si existe un efecto placa,
grp3 <- find.clusters(read_m3, max.n.clust = 4, n.pca = 10, choose = FALSE, stat = "BIC", method = "kmeans")
plot(grp3$Kstat, type = "o", xlab = "numero de grupos (K)",
     ylab = "BIC",
     main = "find.clusters")
grp3$grp
dapcm3 <- dapc(read_m3, grp3$grp, n.pca = 10, n.da = 2) 
scatter(dapcm3)


grp5 <- find.clusters(read_m5, max.n.clust = 4, n.pca = 10, choose = FALSE, stat = "BIC", method = "kmeans")
plot(grp5$Kstat, type = "o", xlab = "numero de grupos (K)",
     ylab = "BIC",
     main = "find.clusters")
grp5$grp
dapcm5 <- dapc(read_m5, grp5$grp, n.pca = 10, n.da = 2) 
scatter(dapcm5)


grp7 <- find.clusters(read_m7, max.n.clust = 9, n.pca = 10, choose = FALSE, stat = "BIC", method = "kmeans")
plot(grp7$Kstat, type = "o", xlab = "numero de grupos (K)",
     ylab = "BIC",
     main = "find.clusters")
grp7$grp
dapcm7 <- dapc(read_m7, grp7$grp, n.pca = 10, n.da = 2) 
scatter(dapcm7)


m3_m2 <- gl.read.vcf("m3_M2.snps.vcf")
m5_M2 <- gl.read.vcf("m5_M2.snps.vcf")

grpm3_m2  <- find.clusters(m3_m2 , max.n.clust = 9, n.pca = 10, choose = FALSE, stat = "BIC", method = "kmeans")
plot(grpm3_m2$Kstat, type = "o", xlab = "numero de grupos (K)",
     ylab = "BIC",
     main = "find.clusters")
grpm3_m2$grp
dapcm3_m2  <- dapc(m3_m2 , grpm3_m2$grp, n.pca = 10, n.da = 2) 
scatter(dapcm3_m2)

grpm5_M2 <- find.clusters(m5_M2, max.n.clust = 9, n.pca = 10, choose = FALSE, stat = "BIC", method = "kmeans")
plot(grpm5_M2$Kstat, type = "o", xlab = "numero de grupos (K)",
     ylab = "BIC",
     main = "find.clusters")
grpm5_M2$grp
dapcm5_M2 <- dapc(m5_M2, grpm5_M2$grp, n.pca = 10, n.da = 2) 
scatter(dapcm5_M2)
