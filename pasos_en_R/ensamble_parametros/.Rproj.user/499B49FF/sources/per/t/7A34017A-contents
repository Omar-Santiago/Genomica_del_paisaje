install.packages("adegenet", dependencies = T)

library(ade4)
library(adegenet)

setwd("C:/Users/OmarS/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/IND_PNDL.2/ensamble/ensamble_parametros")


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

