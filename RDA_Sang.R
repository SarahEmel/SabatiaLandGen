#RDA models for genotype-environment association analysis of soil data in Emel et al. submitted
#code adapted from Forester vignette available at https://popgen.nescent.org/2018-03-27_RDA_GEA.html
#Sarah Emel

library(vegan)
library(psych)

#2 SD = alpha < 0.05

#OUTLIER FUNCTION
#identifying SNPs with 3 stdev cutoff (2 tailed p value = 0.0027)
outliers <- function(x,z){
  #find loadings +/- z SD from mean loading:     
  lims <- mean(x) + c(-1, 1) * z * sd(x) 
  #locus names in these tails:
  x[x < lims[1] | x > lims[2]]           
}


#READ IN AND PREPARE THE GENETIC DATA FOR ANALYSIS
gen <- read.csv(file = "Sang_allele_freqs.csv", row.names = 1) #read in genetic data (population allele counts)
dim(gen) #20 pops x 4842 SNPs
gen<-scale(gen, scale = T, center = T)


#READ IN AND SCREEN THE ENVIRONMENTAL PREDICTORS
data<-read.csv(file="envfinalall.csv") #site means
#Make individual names characters (not factors)
data$Pop <- as.character(data$Pop)
row.names(data)<-data$Pop #changing population names column to row names
data<-data[,-1] #remove population column
names(data) <- c("Cluster", "UTM_N", "UTM_E", "K", "P", "pH", "Ca_Mg", "Mg", "Cr", "Cu", "Ni","Pb", "Zn")
str(data) #look at the structure of the dataframe

#Confirm that genotypes and environmental data are in the same order:
identical(rownames(gen), rownames(data))

pairs.panels(data[,4:13], scale = F, cex.cor = 1)
#K is highly correlated with Cu, so remove it. We also want to remove Pb

pred <- data[,c(5:11,13)] #These are the variables to test in 'all'



#RUN THE RDA
#############################################

#All variables
rda.all <- rda(gen ~ ., data = pred, scale = T) 
rda.all #we have as many constrained axes as there are variables in model

#Get adjusted R2
RsquareAdj(rda.all) #expect low explanatory power

summary(eigenvals(rda.all, model = "constrained"))#Eigenvalues for constrained axes = variance explained
screeplot(rda.all)

#check RDA model for significance
signif.full <- anova.cca(rda.all, permutations = 9999, parallel=getOption("mc.cores")) #default is permutation=999
signif.full

#other checks
signif.axis <- anova.cca(rda.all, by="axis", permutations = 9999, parallel=getOption("mc.cores"))
signif.axis #the first axis is significant, the others not even close
vif.cca(rda.all) #check VIF for predictors all values below 10, so no problem


#PLOT THE RDA
levels(data$Cluster) <- c("State Line", "West Chester") #assigning pop to clusters
eco <- data$Cluster
bg <- c("#33a02c","#ffff33") #nice colors
plot(rda.all, type="n", scaling=3)
points(rda.all, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           #the SNPs
points(rda.all, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[eco]) #the populations
#text(rda.all, scaling=3, display="sites", cex=1) #would make this pretter for a fig
text(rda.all, scaling=3, display="bp", col="#0868ac", cex=1)                           #the predictors
legend("bottomleft", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)


#IDENTIFY CANDIDATE SNPS INVOLVED IN LOCAL ADAPTATION
load.rda <- scores(rda.all, choices=1, display="species")  #species scores for the first constrained axis

#look at histograms of loadings
hist(load.rda[,1], main="Loadings on RDA1") #SNPs loading at the tails show a relationship with env
cand1 <- outliers(load.rda[,1],2)
length(cand1) #this is the number of outliers

#Let’s make a single data frame with the axis, SNP name, & loading
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1)) #only first axis significant so no need to combine
colnames(cand1) <- c("axis","snp","loading")
cand.rda <- cand1
cand.rda$snp <- as.character(cand.rda$snp)

write.csv(cand.rda, "rdaoutliers.csv")
