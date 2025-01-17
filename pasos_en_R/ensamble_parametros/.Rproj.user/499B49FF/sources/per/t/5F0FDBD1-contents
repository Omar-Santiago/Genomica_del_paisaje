################################
### Supplementary material 1 ###
################################

# These command lines were used to create GLMMs with binomial and gamma distributions
# and then, to to determine if there is a significant difference between assemblies.

# The following libraries were used
library(MASS)
library(lme4)
library(lmerTest)
library(nlme)
library(outliers)
library(stringr)
library(tidyr)

# Read the data. Tables were computed using Error_EuclDist_reps.r
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

# Split the column "param", change the column names and remove observations with high missing data
split1<-data.frame(str_split_fixed(data$param, "_", 2))
split2<-data.frame(separate(split1, X1, into = c("Treatm", "Parameter"), sep = 2, remove = T))
data<-cbind(split2, data[,c(2:14)])
rm(split1, split2)
names(data)[c(1:6)]<-c("Treat", "Param", "Value", "Pair", "SNP.error.rt", "EuclDist")
data <- data[data$MissTotProp < 0.5, ]

# Convert SNP error rates to binary data
data$nonZero <- ifelse(data$SNP.error.rt > 0, 1, 0) 
# The data.frame is ready!

##################
### CLUSTERING ###
##################

# Subset the data
data_clust <-data[data$Param == "clust",] # Keep only clustering data 
data_clust$Value<-as.numeric(as.character(data_clust$Value)) # Param value is numeric in this case
data_clust$Treat<-factor(data_clust$Treat) # while treatment should be a factor

# Check and remove outliers: Higher error rates may be linked to specific sample pairs.
# Visual revision of the data
# boxplot(data_clust$SNP.error.rt ~ data_clust$Treat)

datosx1<-data_clust[data_clust$Treat == "X1",] # We detected outliers in the 1SNP dataset
minOutlier<-min(boxplot.stats(datosx1$SNP.error.rt)$out) # Identify outliers
outliers<-factor(as.character((datosx1[datosx1$SNP.error.rt >= minOutlier,])$Pair)) # Identify which pairs are far away from minOutlier
# Remove outliers
data_clustrm <-data_clust[!(data_clust$Pair %in% levels(outliers)),]

### GLMMs for SNP error rate

# Previous inspection lineal models suggested heteroscedasticity (La varianza no es homogenea)
# One solutions would be to use GLMs or generalized least square

# Although the residual distribution is normal, the dependent variable has too much zeros
# Thus, we decided to transform data to binomial values.
# Furthermore, we included the treatment (5 or 1 SNP/locus) and the pair IDs as random variables
# We implemented two models. However, Table S1 was filled with results from the second function.  

m3bin <- glmer(nonZero ~ Value + Treat + (1|Pair),
               family = binomial(link = "logit"), data = data_clustrm)

mbin_clust<-glmmPQL(nonZero ~ Value + Treat,
                    random= list(Pair= ~ 1), 
                    family = binomial(link = logit),
                    data= data_clustrm)

###  GLMMs for Euclidean distances

# Again, check for outliers associated to specific sample pairs
boxplot(data_clust$EuclDist ~ data_clust$Treat)
# We did not observe outliers so, we can start fitting a model
# In this case, we chose a gamma distribution as the dependent variable 
# is left-limited (at zero) and partially skewed.

m3gam <- glmer((EuclDist+1) ~ Value + Treat +(1|Pair),
               family = Gamma(link = log),
               data = data_clust)
  
mgam_clust<-glmmPQL((EuclDist+1) ~ Value + Treat,
                    random= list(Pair= ~ 1), 
                    family = Gamma(link = log),
                    data= data_clust)

#############
### MAXHS ###
#############

### Subset the data
data_maxhs <-data[data$Param == "maxHs",] # Keep just maxHs data 
data_maxhs$Value<-as.numeric(data_maxhs$Value) # Values are categories
data_maxhs$Treat<-factor(data_maxhs$Treat) # Treatment have to be a factor

# Visual revision of the data
# boxplot(data_maxhs$SNP.error.rt ~ data_maxhs$Treat)

datosx1<-data_maxhs[data_maxhs$Treat == "X1",] # We detected outliers in the 1SNP dataset
minOutlier<-min(boxplot.stats(datosx1$SNP.error.rt)$out) # Identify outliers
outliers<-factor(as.character((datosx1[datosx1$SNP.error.rt >= minOutlier,])$Pair)) #Identify pairs

# Remove outliers
data_maxhsrm <-data_maxhs[!(data_maxhs$Pair %in% levels(outliers)),]


### GLMMs for SNP error rate 

m3bin <- glmer(nonZero ~ Value + Treat + (1|Pair),
               family = binomial(link = "logit"),
               data = data_maxhsrm)

mbin_maxhs<-glmmPQL(nonZero ~ Value + Treat,
                    random= ~ 1 | Pair,
                    family = binomial(link = logit),
                    data= data_maxhsrm)

### GLMMs for Euclidean distance
# Visual revision of the data
# boxplot(data_maxhs$EuclDist ~ data_maxhs$Treat)

m3gam <- glmer((EuclDist+1) ~ Value + Treat + (1|Pair),
               family = Gamma(link = log), data = data_maxhs)

mgam_maxHs<-glmmPQL((EuclDist+1) ~ Value + Treat,
                    random= list(Pair= ~ 1), 
                    family = Gamma(link = log),
                    data= data_maxhs)

#################
### MEANDEPTH ###
#################

### Subset the data
data_min <-data[data$Param == "mindepth",] # Keep just maxHs data 
data_min$Value<-as.numeric(data_min$Value) # Values are categories
data_min$Treat<-factor(data_min$Treat) # Treatment have to be a factor

# Visual revision of the data
# boxplot(data_min$SNP.error.rt ~ data_min$Treat)
# Again, only the 1SNP/locus dataset shows outliers  

datosx1<-data_min[data_min$Treat == "X1",] # We detected outliers in the 1SNP dataset
minOutlier<-min(boxplot.stats(datosx1$SNP.error.rt)$out) # Identify outliers
outliers<-factor(as.character((datosx1[datosx1$SNP.error.rt >= minOutlier,])$Pair)) #Identify pairs
# Remove outliers
data_minrm <-data_min[!(data_min$Pair %in% levels(outliers)),]

### GLMM for SNP error rate
m3bin <- glmer(nonZero ~ Value + Treat + (1|Pair),
               family = binomial(link = "logit"),
               data = data_minrm)

mbin_depth<-glmmPQL(nonZero ~ Value + Treat,
                    random= list(Pair= ~ 1), 
                    family = binomial(link = logit),
                    data= data_minrm)


### GLMMs for Euclidean distances
# Visual revision of the data
#boxplot(data_min$EuclDist ~ data_min$Treat)

### GLMMs para Euclidean distances
m3gam <- glmer((EuclDist+1) ~ Value + Treat + (1|Pair),
               family = Gamma(link = log),
               data = data_min)

mgam_depth<-glmmPQL((EuclDist+1) ~ Value + Treat,
                    random= list(Pair= ~ 1), 
                    family = Gamma(link = log),
                    data= data_min)

#################
### Figure S1 ###
#################

# a) SNP error rate
data2<-rbind(data_clust,data_maxhs,data_min)

p<-ggplot(data2, aes(x=as.factor(Value),
                     y=SNP.error.rt, 
                     fill=as.factor(Treat))) +
                     geom_boxplot() + theme_classic() +
                     labs(x ="Parameter value", y = "Freq SNP mismatches") +
                     scale_fill_manual(values = c("white", "grey"))+
                     theme(legend.position = "none",
                        axis.title = element_text(size = 16),
                        axis.text = element_text(size = 12),
                        panel.grid = element_blank(),
                        strip.text = element_text(size=18))+
                     facet_wrap(. ~ Param, scales = "free_x")
                     

# b) Euclidean distances

q<-ggplot(data2, aes(x=as.factor(Value),
                     y=EuclDist, 
                     fill=as.factor(Treat))) +
                     geom_boxplot() + theme_classic() +
                     labs(x ="Parameter value", y = "Euclidean distance") +
                     scale_fill_manual(values = c("white", "grey"))+
                     theme(legend.position = "none",
                        axis.title = element_text(size = 16),
                        axis.text = element_text(size = 12),
                        panel.grid = element_blank(),
                        strip.text = element_text(size=18))+
                     facet_wrap(. ~ Param, scales = "free_x")

### Number of SNPs
        
SNPs<-read.csv("SNP_count_tests.csv")
SNPs<-SNPs[SNPs$Assembly=="TMVB", ] # Keep just TMVB data
Guspalette =c("black", "grey") # Set two colors

tt<-ggplot(SNPs, aes(x = Param_value, y =true_snp,
                 color = factor(Treatment), 
                 group = Treatment)) +
                 geom_point() + geom_line() +
                 facet_grid(. ~ Param, scales = "free_x",) +
                 theme_classic()+
                 labs(x ="Parameter value", y = "SNPs") +
                 scale_color_manual(values = Guspalette) +
                 theme(legend.position = "none",
                    axis.title = element_text(size = 16),
                    axis.text = element_text(size = 12),
                    panel.grid = element_blank(),
                    strip.text = element_text(size=18))
  

pdf("Figure_S1.pdf", width = 8, height = 5)
p
q
tt
dev.off()

# Alternatively, you can use Boxplot_reps.r to create Figure S1