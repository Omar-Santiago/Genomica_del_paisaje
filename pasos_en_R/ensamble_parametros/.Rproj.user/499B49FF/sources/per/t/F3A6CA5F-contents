#####################################################################
# These command lines produce Figure S1, Supplementary material 1 ###
#####################################################################

library(plotly)
require(scales)
require(gridExtra)
library(ggplot2)
library(stringr)
library(tidyr)

# Import SNP dataset information 
snperror<-read.table("5SNP_TMVBEqSam/SNPerrorRates.txt", sep="\t", header = T)
allerror<-read.table("5SNP_TMVBEqSam/LocAllError.txt", sep="\t", header = T)
eucldist<-read.table("5SNP_TMVBEqSam/EuclDist.txt", sep="\t", header = T)

SNPs_5<-cbind(snperror, eucldist[,3], allerror[,c(3:12)])
rm(snperror, eucldist, allerror)

snperror<-read.table("1SNP_TMVBEqSam/SNPerrorRates.txt", sep="\t", header = T)
allerror<-read.table("1SNP_TMVBEqSam/LocAllError.txt", sep="\t", header = T)
eucldist<-read.table("1SNP_TMVBEqSam/EuclDist.txt", sep="\t", header = T)

SNPs_1<-cbind(snperror, eucldist[,3], allerror[,c(3:12)])
rm(snperror, eucldist, allerror)

data<-rbind(SNPs_5, SNPs_1)
rm(SNPs_1, SNPs_5)
# please, notice that rows match each other between 5 and 1 SNP datasets.


split1<-data.frame(str_split_fixed(data$param, "_", 2))
split2<-data.frame(separate(split1, X1, into = c("Treatm", "Parameter"), sep = 2, remove = T))
data<-cbind(split2, data[,c(2:14)])
rm(split1, split2)
names(data)[c(1:6)]<-c("Treat", "Param", "Value", "Pair", "SNP.error.rt", "EuclDist")

# Remove observations with high missing data (>50%)
data <- data[data$MissTotProp < 0.5, ]
# Transform SNP error rates to binary data
data$nonZero <- ifelse(data$SNP.error.rt > 0, 1, 0) 

##################################################################
### Boxplot graphics for SNP error rate and Euclidean distance ###
##################################################################

# Re-order the levels in data
data$Value <- factor(data$Value, levels=c("070","075","080","085","090","095","2","3","4","5","6","7","8","9","10"))

# SNP error rate plot 
p <-ggplot(data, aes(x=as.factor(Value), y=SNP.error.rt, fill=Treat))+
geom_boxplot(varwidth = TRUE) +
#geom_jitter(width = 0.2) +
theme_bw() +
#stat_summary(fun.y=mean, geom="point", shape=19, size=2, color="darkred") +
facet_wrap(. ~ Param, scales = "free_x",)

# Create a color palette
Guspalette =c("#1dabe6", "#1c366a","#e43034","#fc4e51","#af060f","#bdbdbd")

# and just plot
p +labs(x ="Parameter value", y = "SNP error") + 
  scale_fill_manual(values = Guspalette, name="No SNPs/locus") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=10,face="bold"), strip.text = element_text(size = 14))

# and now do the same for Euclidean distances 
p <-ggplot(data, aes(x=as.factor(Value), y=EuclDist, fill=Treat))+
geom_boxplot(varwidth = TRUE) +
#geom_jitter(width = 0.2) +
theme_bw() +
#stat_summary(fun.y=mean, geom="point", shape=19, size=2, color="darkred") +
facet_wrap(. ~ Param, scales = "free_x",)

# and just plot
p +labs(x ="Parameter value", y = "Euclidean distance") + 
  scale_fill_manual(values = Guspalette, name="No SNPs/locus") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=10,face="bold"), strip.text = element_text(size = 14))

###########################
### Number of SNPs plot ###
###########################

# Read and prepare data
count_5<-read.table("5SNP_TMVBEqSam/5SNP_conteos.txt", header = F, sep="\t")
count_5<-data.frame(separate(count_5, V1, into = c("Parameter", "Value"), sep = "_", remove = T))
count_5$Treat<-"5"
count_1<-read.table("1SNP_TMVBEqSam/1SNP_conteos.txt", header = F, sep="\t")
count_1<-data.frame(separate(count_1, V1, into = c("Parameter", "Value"), sep = "_", remove = T))
count_1$Treat<-"1"
data<-rbind(count_5, count_1)
rm(count_5, count_1)
data$Value <- factor(data$Value, levels=c("070","075","080","085","090","095","2","3","4","5","6","7","8","9","10"))

# and just plot 
p <-ggplot(data, aes(x=Value, y=V2, group=as.factor(Treat), color=as.factor(Treat))) + 
      geom_line() +
      geom_point() +
      #geom_jitter(width = 0.2) +
      theme_bw() +
      #stat_summary(fun.y=mean, geom="point", shape=19, size=2, color="darkred") +
      facet_wrap(. ~ Parameter, scales = "free_x")
    
    p +labs(x ="Parameter value", y = "No SNPs") + 
      scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"), name="SNPs/locus") +
      theme(axis.text=element_text(size=12), axis.title=element_text(size=10,face="bold"), strip.text = element_text(size = 14))
      