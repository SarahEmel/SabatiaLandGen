#MLPE models for landscape genetic analysis of least cost path and Circuitscape resistance in Emel et al. submitted
#Sarah Emel


library(MuMIn)
library(lme4)

###read in data###
#check filenames/formats after any changes

#genetic distance data
fst<-read.csv(file="FSToutliersremoved.csv",header=T) #Weir and Cockerham's Fst from diveRsity
fst<-fst[lower.tri(fst)]
gdist<-fst

#populations, euclidean distance, least-cost paths
#make sure all pops listed in both pop1 and pop2
pop<-read.csv(file="pops.csv",header=T)
lcpcd<-read.csv(file="Sang_costdist_toimport.csv",header=T)

###create dataframe, subset as necessary###
alldata<-cbind(gdist,pop,lcpcd)
dataseparate<-alldata[which(alldata$Cluster=='S'|alldata$Cluster=='N'),]

############################################
###assign the dataframe for this analysis###
data<-alldata
# data<-dataseparate
############################################

#populations
pop1<-data$pop1
pop2<-data$pop2

#scale and center distances
scaled<-scale(data[,8:ncol(data)],center=T,scale=T)
scaled.data<-cbind(data[,1:7],scaled)
data<-scaled.data

##formula part 1
Zl <- lapply(c("pop1","pop2"), function(nm) Matrix:::fac2sparse(data[[nm]],"d", drop=FALSE))
ZZ <- Reduce("+", Zl[-1], Zl[[1]])

##formula part 2

#MLPE functions
#use MLPE for model fitting, examining coefficients, CIs
#use MLPEnoREML for model comparison
MLPE <- function(variables, data) {
  mod2 <- lme4::lFormula(variables, data = data, REML = TRUE)
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  mod_2 <- lme4::mkMerMod(environment(dfun), opt, mod2$reTrms,fr = mod2$fr)
  mod2$reTrms$Zt <- ZZ
  # Refit the model
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  modelout <- lme4::mkMerMod(environment(dfun), opt, mod2$reTrms,fr = mod2$fr)
  return(modelout)
}

MLPEnoREML <- function(variables, data) {
  mod2 <- lme4::lFormula(variables, data = data, REML = FALSE)
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  mod_2 <- lme4::mkMerMod(environment(dfun), opt, mod2$reTrms,fr = mod2$fr)
  mod2$reTrms$Zt <- ZZ
  # Refit the model
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  modelout <- lme4::mkMerMod(environment(dfun), opt, mod2$reTrms,fr = mod2$fr)
  return(modelout)
}

#heading for the output file
header <- noquote(paste("AICc", "BIC", "R2m", "R2c", "Var",sep=" ")) 
write.table(header, file = "lcpall.txt", append=T, row.names=F, col.names=F, quote = F) 

#do model comparison with MLPEnoREML function
for (i in (8:ncol(data))) {
  fit <- MLPEnoREML(gdist ~ data[, i] + (1|pop1), data=data)
  anova(fit)
  b <- BIC(fit)
  a <- AICc(fit, k = length(unique(pop1))) #this make k = # pops, not # pairs
  r <- r.squaredGLMM(fit)
  rm <- r[1]
  rc <- r[2]
  out<- noquote(paste(a, b, rm, rc, colnames(data)[i], sep =" "))
  write.table(out ,file = "lcpall.txt", append=T,row.names=F,col.names=F, quote = F) 
}


#get confidence intervals for predictors
header2 <- noquote(paste("Var","Est","2.5%", "97.5%", sep=" ")) 
write.table(header2, file = "lcpall_ci.txt", append=T, row.names=F, col.names=F, quote = F)

for (i in (8:ncol(data))) { #may need to adjust code for multivariate models
  fit <- MLPE(gdist ~ data[, i] + (1|pop1), data=data)
  est <- coef(summary(fit))
  ci <- confint(fit, level = 0.95, method = "Wald")
  out<- noquote(paste(colnames(data)[i], est[2,1], ci[4,1], ci[4,2], sep =" "))
  write.table(out ,file = "lcpall_ci.txt", append=T,row.names=F,col.names=F, quote = F) 
}