##################################
# Empirical analysis evaluating 
# Granger Causality between corruption
# perceptions and electoral 
# volatility. It uses a Bayesian
# bivariate normal model (VAR1).
#
# Author: Santiago Olivella
#         (olivella@wustl.edu)
# Date: July 2012.
###################################


#Load libraries and set working directory
library(rjags)
library(mice)
library(mvnormtest)
library(MCMCpack)
library(emdbook)
library(foreign)
library(tseries)
library(plyr)
library(lattice)
load.module("dic")

setwd("~/Dropbox/Corruption/DataAnalysis/")
rm(list=ls())



#Set random seed (for future replicability)
set.seed(831213)


#Load data 
# Update (March 2012) : Use Efficiency Measure instead of corruption.

corvol <- read.dta("CorVolSerieswControls2.dta")
unemployment <- read.csv("unemployment.csv")
corvol <- merge(corvol,unemployment,all=TRUE)

#corvol <- read.dta("GovEff.dta")
#names(corvol)[6] <- "pedersen" 
corvol <- corvol[order(corvol$country,corvol$year),]

corvol$unemploy_lag1 <- ave(corvol$unemploy,corvol$country,FUN=function(x)c(NA,x[-length(x)]))
corvol$inflation_lag1 <- ave(corvol$inflation,corvol$country,FUN=function(x)c(NA,x[-length(x)]))
corvol$chng_gdppc_lag1 <- ave(corvol$chng_gdppc,corvol$country,FUN=function(x)c(NA,x[-length(x)]))

#corvol.sub <- corvol


corvol.sub <- subset(corvol, select=c(country
                                      ,year
                                      ,TermCounter
                                      ,pedersen
                                      ,pedersen_lag1
                                      ,pedersen_lag2
                                      ,type_a
                                      ,type_a_lag1
                                      ,type_a_lag2
                                      ,type_b
                                      ,type_b_lag1
                                      ,type_b_lag2
                                      ,leg_avg
                                      ,leg_avg_lag1
                                      ,leg_avg_lag2
                                      ,polity
                                      ,age
                                      ,enp
                                      ,eng_legal
                                      ,ethling85
                                      ,parl
                                      ,gdppc
                                      ,chng_gdppc
                                      ,chng_gdppc_lag1
                                      ,inflation
                                      ,inflation_lag1
                                      ,unemploy
                                      ,unemploy_lag1
                                      ,pr
                                      ,cl
                                      ,pvsi_main
                                      ,open_econ
                                      ))

corvol.sub$dem <- corvol.sub$polity>=7
# Impute data. M=5
corvol.sub.i <- corvol.sub
corvol.sub.mi <- mice(corvol.sub.i)


# We will get Bayesian estimates
# using each imputed dataset.
# This is equivalent to sampling 
# from the joint posterior of parameters
# and missing values. Effectively,
# it incorporates uncertainity 
# resulting from imputed data.

n.iters <- 1e4 # Gibbs sampler iterations
thin.int <- 15 # Gibbs sampler thinning interval

#Loop over multiply imputed datasets
estimation.fun <- function(target.variableCor,target.variableVol,model,dem.test=FALSE){
  results <- list()
  for (i in 1:5){
    
    corvol.sub.c <- complete(corvol.sub.mi, action = i)
    corvol.sub <- corvol.sub.c
    
    #Calculate logged values
    corvol.sub$cor <- corvol.sub[,match(target.variableCor,names(corvol.sub))]
    corvol.sub$cor_lag1 <- corvol.sub[,match(paste(target.variableCor,"_lag1",sep=""),names(corvol.sub))]
    corvol.sub$vol <- corvol.sub[,match(target.variableVol,names(corvol.sub))]
    corvol.sub$vol_lag1 <- corvol.sub[,match(paste(target.variableVol,"_lag1",sep=""),names(corvol.sub))]
    
    
    ## Predictor matrices with control variables
    # for replication purposes only. 
    country.num <- as.numeric(as.factor(corvol.sub$country))
    
    XCor <- as.matrix(cbind(
      subset(corvol.sub
             ,ifelse(dem.test,logical(dem==dem.test),TRUE)
             ,select=c(
               #Cross-lags
               vol_lag1 
               #Self-lags
               ,cor_lag1
             ))
    )
    )
    XCorEx <- as.matrix(cbind(
      subset(corvol.sub
             ,ifelse(dem.test,logical(dem==dem.test),TRUE)
             ,select=c(
               eng_legal
               ,ethling85
               ,parl
               ,gdppc
               ,pr
               ,cl
               ,open_econ
               ,pvsi_main
             ))
    )
    )
    XCorEx <- scale(as.data.frame(XCorEx))
    
    XVol <- as.matrix(cbind(
      subset(corvol.sub
             ,ifelse(dem.test,logical(dem==dem.test),TRUE)                        
             ,select=c(
               #Cross-lags
               cor_lag1
               #Self-lags
               ,vol_lag1
             ))
    )
    )
    
    XVolEx <- as.matrix(cbind(
      subset(corvol.sub
             ,ifelse(dem.test,logical(dem==dem.test),TRUE)                      
             ,select=c(
               ethling85
               ,age
               ,enp
               ,chng_gdppc
               ,inflation
               ,open_econ
               ,pvsi_main
             ))
    )
    )
    XVolEx <- scale(as.data.frame(XVolEx))
    
    XBeta <- by(subset(corvol.sub
                       ,ifelse(dem.test,logical(dem==dem.test),TRUE)
                       ,select=c(polity,enp,leg_avg))
                ,country.num,colMeans)
    XBeta <- do.call(rbind,XBeta)
                
    # Response matrix. First column: corruption; Second Column: volatility
    Y <- as.matrix(subset(corvol.sub
                          ,ifelse(dem.test,logical(dem==dem.test),TRUE)
                          ,select=c(cor,vol)
                          #,select=c(post_gov_eff,pedersen)
    ))
    
    #Standardize series
    XCor <- scale(as.data.frame(XCor),center=FALSE,scale=apply(as.data.frame(XCor),2,sd,na.rm=TRUE))
    XVol <- scale(as.data.frame(XVol),center=FALSE,scale=apply(as.data.frame(XVol),2,sd,na.rm=TRUE))
    Y <- scale(as.data.frame(Y),center=FALSE,scale=apply(as.data.frame(Y),2,sd,na.rm=TRUE))
    
    # Set parameters fed to JAGS (see CoVolModel.jag)  
    N <- dim(XCor)[1]
    PC <- dim(XCor)[2]
    PV <- dim(XVol)[2]
    PCx <- dim(XCorEx)[2]
    PVx <- dim(XVolEx)[2]
    C <- max(country.num)
    
    
    # Data list (for JAGS)
    corvolData <- switch(model,
                         list(N = N
                       ,beta0=rep(0,PC+PV)
                       ,PC = PC
                       ,PV = PV
                       ,C = C
                       ,country=country.num
                       ,XCor = XCor
                       ,XVol = XVol
                       ,Y = Y
                       ,Ri = solve(cov(Y))
                       ,Rbeta = diag(1,PC+PV)
                       ,TauMuBeta = diag(0.01,PC+PV)
                       ),
                         controls = list(N = N
                                         ,beta0=rep(0,PC+PV)
                                         ,PC = PC
                                         ,PV = PV
                                         ,PCx = PCx
                                         ,PVx = PVx
                                         ,C = C
                                         ,country=country.num
                                         ,XCor = XCor
                                         ,XVol = XVol
                                         ,Y = Y
                                         ,Ri = solve(cov(Y))
                                         ,Rbeta = diag(1,PC+PV)
                                         ,TauMuBeta = diag(0.001,PC+PV)
                                         ,XCorEx = XCorEx
                                         ,XVolEx = XVolEx
                                         ,psi0 = rep(0,PCx+PVx)
                                         ,TauPsi = diag(0.001,PCx+PVx)
                         ),
                         interaction = list(N = N
                                            ,beta0=rep(0,PC+PV-1)
                                            ,PC = PC
                                            ,PV = PV
                                            ,C = C
                                            ,country=country.num
                                            ,XCor = XCor
                                            ,XVol = XVol
                                            ,XBeta = XBeta
                                            ,Y = Y
                                            ,Ri = solve(cov(Y))
                                            ,Rbeta = diag(1,PC+PV-1)
                                            ,TauMuBeta = diag(0.001,PC+PV-1)
                                            ,gamma0 = rep(0,dim(XBeta)[2])
                                            ,TauGamma = diag(0.001,dim(XBeta)[2])
                         )
                         )
    JagsModel <- switch(model
                        ,"CorVolModel.jag"
                        ,controls="CorVolModelEx.jag"
                        ,interaction="CorVolModelInter.jag")
                         
    # Data list for null model
    corvolDataNull <- list(N = N
                           ,beta0=rep(0,PC+PV-2)
                           ,PC = PC-1
                           ,PV = PV-1
                           ,C = C
                           ,country=country.num
                           ,XCor = as.matrix(XCor[,-1])
                           ,XVol = as.matrix(XVol[,-2])
                           ,Y = Y
                           ,Ri = solve(cov(Y))
                           ,Rbeta = diag(1,PC+PV-2)
                           ,TauMuBeta = diag(0.001,PC+PV-2)
                           )
    
    # JAGS call. 2 chains per multiply imputed dataset.
    corvolModel <- jags.model(file=JagsModel,
                              data=corvolData,
                              n.chains=2
                              #,inits=inits.corvol
                              )
    corvolModelNull <- jags.model(file="CorVolModel.jag",
                                  data=corvolDataNull,
                                  n.chains=2
                                  #,inits=inits.corvolNull
                                  )
    update(corvolModel,n.iters/2)
    
    var.names <- switch(model,
                        c(#"alphaCor"
                          #,"alphaVol"
                          "rho"
                          ,"residuals"
                          ,"sigma.acor" 
                          ,"sigma.avol"
                          ,"mu.avol"
                          ,"mu.acor"
                          ,"mu.beta"
                          ,"tau.beta"
                          ,"sigma2"
                        )
                        ,controls=c(#"alphaCor"
                           #,"alphaVol"
                           "rho"
                           ,"residuals"
                           ,"sigma.acor" 
                           ,"sigma.avol"
                           ,"mu.avol"
                           ,"mu.acor"
                           ,"mu.beta"
                           ,"tau.beta"
                           ,"sigma2"
                           ,"psivec")
                        ,interaction=c(#"alphaCor"
                                       #,"alphaVol"
                                       "rho"
                                       ,"residuals"
                                       ,"sigma.acor" 
                                       ,"sigma.avol"
                                       ,"mu.avol"
                                       ,"mu.acor"
                                       ,"mu.beta"
                                       ,"tau.beta"
                                       ,"sigma2"
                                       ,"gammavec"))
    
    #Get posterior samples and store them in results list.
    results[[i]] <- coda.samples(corvolModel,
                                 n.iter=n.iters,
                                 thin=thin.int,
                                 variable.names=var.names)[[sample(2,1)]]
  } ## End of for loop.
  return(list(corvolModel=corvolModel,corvolModelNull=corvolModelNull,XVol=XVol,XCor=XCor,PC=PC,PV=PV,results=results,Y=Y))
}

### All models
#Standard model using pedersen volatility (REPORTED)
std.ped.model <- estimation.fun("leg_avg","pedersen","standard")
res.std.ped <- as.mcmc.list(std.ped.model$results)
summary.res.std.ped <- summary(res.std.ped,c(0.05,0.5,0.95))

## Robustness checks (ONLINE APPENDIX)

#Standard model using pedersen and high polity subset 
#(must uncomment dem==TRUE line in above function)
polity.ped.model <- estimation.fun("leg_avg","pedersen","standard")
res.polity.ped <- as.mcmc.list(polity.ped.model$results)
summary.res.polity.ped <- summary(res.polity.ped,c(0.05,0.5,0.95))

#Standard model using type a volatility
std.typea.model <- estimation.fun("leg_avg","type_a","standard")
res.std.typea <- as.mcmc.list(std.typea.model$results)
summary.res.std.typea <- summary(res.std.typea,c(0.05,0.5,0.95))

#Standard model using type b volatility
std.typeb.model <- estimation.fun("leg_avg","type_b","standard")
res.std.typeb <- as.mcmc.list(std.typeb.model$results)
summary.res.std.typeb <- summary(res.std.typeb,c(0.05,0.5,0.95))

#Interaction model 
int.model <- estimation.fun("leg_avg","pedersen","interaction")
res.int.model <- as.mcmc.list(int.model$results)
summary.res.int.model <- summary(res.int.model,c(0.05,0.5,0.95))

#Controls model 
cot.model <- estimation.fun("leg_avg","pedersen","controls")
res.cot.model <- as.mcmc.list(cot.model$results)
summary.res.cot.model <- summary(res.cot.model,c(0.05,0.5,0.95))


#Summaries for tables
round(summary.res.std.ped$quantiles[c(1:6,345:347,352,357,362,367),],3)
round(summary.res.polity.ped$quantiles[c(1:6,345:347,352,357,362,367),],3)
round(summary.res.std.typea$quantiles[c(1:6,345:347,352,357,362,367),],3)
round(summary.res.std.typeb$quantiles[c(1:6,345:347,352,357,362,367),],3)
round(summary.res.int.model$quantiles[c(1:8,347:349,354,358,362),],3)
round(summary.res.cot.model$quantiles[c(1:21,360:362,367,372,377,382),],3)



### Something like a LR test w.r.t null model
test1.std <- dic.samples(std.ped.model$corvolModel,n.iter=1e3)
test2.std <- dic.samples(std.ped.model$corvolModelNull,n.iter=1e3)
test1.std-test2.std

test1.pol <- dic.samples(polity.ped.model$corvolModel,n.iter=1e3)
test2.pol <- dic.samples(polity.ped.model$corvolModelNull,n.iter=1e3)
test1.pol-test2.pol

test1.typea <- dic.samples(std.typea.model$corvolModel,n.iter=1e3)
test2.typea <- dic.samples(std.typea.model$corvolModelNull,n.iter=1e3)
test1.typea-test2.typea

test1.typeb <- dic.samples(std.typeb.model$corvolModel,n.iter=1e3)
test2.typeb <- dic.samples(std.typeb.model$corvolModelNull,n.iter=1e3)
test1.typeb-test2.typeb

test1.int <- dic.samples(int.model$corvolModel,n.iter=1e3)
test2.int <- dic.samples(int.model$corvolModelNull,n.iter=1e3)
test1.int-test2.int

test1.cont <- dic.samples(cot.model$corvolModel,n.iter=1e3)
test2.cont <- dic.samples(cot.model$corvolModelNull,n.iter=1e3)
test1.cont-test2.cont


# # # #Obtain residuals
residuals.std <- cbind(summary.res.std.ped$quantiles[grep("residuals\\[\\d*,1\\]",
                                                  dimnames(summary.res.std.ped$quantiles)[[1]])
                                             ,2]
                  ,summary.res.std.ped$quantiles[grep("residuals\\[\\d*,2\\]"
                                                   ,dimnames(summary.res.std.ped$quantiles)[[1]])
                                            ,2])
residuals.pol <- cbind(summary.res.polity.ped$quantiles[grep("residuals\\[\\d*,1\\]",
                                                          dimnames(summary.res.polity.ped$quantiles)[[1]])
                                                     ,2]
                       ,summary.res.polity.ped$quantiles[grep("residuals\\[\\d*,2\\]"
                                                           ,dimnames(summary.res.polity.ped$quantiles)[[1]])
                                                      ,2])
residuals.typea <- cbind(summary.res.std.typea$quantiles[grep("residuals\\[\\d*,1\\]",
                                                          dimnames(summary.res.std.typea$quantiles)[[1]])
                                                     ,2]
                       ,summary.res.std.typea$quantiles[grep("residuals\\[\\d*,2\\]"
                                                           ,dimnames(summary.res.std.typea$quantiles)[[1]])
                                                      ,2])
residuals.typeb <- cbind(summary.res.std.typeb$quantiles[grep("residuals\\[\\d*,1\\]",
                                                              dimnames(summary.res.std.typeb$quantiles)[[1]])
                                                         ,2]
                         ,summary.res.std.typeb$quantiles[grep("residuals\\[\\d*,2\\]"
                                                               ,dimnames(summary.res.std.typeb$quantiles)[[1]])
                                                          ,2])
residuals.int <- cbind(summary.res.int.model$quantiles[grep("residuals\\[\\d*,1\\]",
                                                              dimnames(summary.res.int.model$quantiles)[[1]])
                                                         ,2]
                         ,summary.res.int.model$quantiles[grep("residuals\\[\\d*,2\\]"
                                                               ,dimnames(summary.res.int.model$quantiles)[[1]])
                                                          ,2])
residuals.cot <- cbind(summary.res.cot.model$quantiles[grep("residuals\\[\\d*,1\\]",
                                                              dimnames(summary.res.cot.model$quantiles)[[1]])
                                                         ,2]
                         ,summary.res.cot.model$quantiles[grep("residuals\\[\\d*,2\\]"
                                                               ,dimnames(summary.res.cot.model$quantiles)[[1]])
                                                          ,2])
#ACF of residuals
acf(residuals.std[,1],4,ci=0.99)
acf(residuals.std[,2],4,ci=0.99)

# # # # Dimension by dimension R squared
temp.sq <- outer(as.matrix(std.ped.model$Y),colMeans(std.ped.model$Y),"-")
rsqr.std <- 1-(colSums(residuals.std^2)/colSums(cbind(temp.sq[,1,1],temp.sq[,2,2])^2))

temp.sq <- outer(as.matrix(polity.ped.model$Y),colMeans(polity.ped.model$Y),"-")
rsqr.pol <- 1-(colSums(residuals.pol^2)/colSums(cbind(temp.sq[,1,1],temp.sq[,2,2])^2))

temp.sq <- outer(as.matrix(std.typea.model$Y),colMeans(std.typea.model$Y),"-")
rsqr.typea <- 1-(colSums(residuals.typea^2)/colSums(cbind(temp.sq[,1,1],temp.sq[,2,2])^2))

temp.sq <- outer(as.matrix(std.typeb.model$Y),colMeans(std.typeb.model$Y),"-")
rsqr.typeb <- 1-(colSums(residuals.typeb^2)/colSums(cbind(temp.sq[,1,1],temp.sq[,2,2])^2))

temp.sq <- outer(as.matrix(int.model$Y),colMeans(int.model$Y),"-")
rsqr.int <- 1-(colSums(residuals.int^2)/colSums(cbind(temp.sq[,1,1],temp.sq[,2,2])^2))

temp.sq <- outer(as.matrix(cot.model$Y),colMeans(cot.model$Y),"-")
rsqr.cot <- 1-(colSums(residuals.cot^2)/colSums(cbind(temp.sq[,1,1],temp.sq[,2,2])^2))



# # # # Test of normality of residuals.
mshapiro.test(t(residuals.std))

###### Mean (Orthogonalized,opt.) IRFs
ortho.irf <- function(mcmc.li,sigma.ind=NULL,beta.ind){
    chain.ind <- sample(length(mcmc.li),1)
    sample.ind <- sample(dim(mcmc.li[[chain.ind]])[1],1e2)
    mcmc.li <- mcmc.li[[chain.ind]][sample.ind,]
    Omega <- array(NA,c(2,2,7,1e2))
    for(i in 1:1e2){
      sigma.u <- matrix(mcmc.li[i,sigma.ind],ncol=2) #Only for ortho
      P <- chol(sigma.u) #Only for ortho
      phi <- matrix(mcmc.li[i,beta.ind],ncol=2,byrow=TRUE)
      Omega[,,1,i] <- diag(2)  #change to P to get ortho
      for(j in 2:7){
        Omega[,,j,i] <- phi%*%Omega[,,j-1,i]
      }
    }
    return(Omega)
    #return(aperm(aaply(Omega,c(1,2,4),cumsum),c(1,2,4,3))) #For cumulative responses 
}

irfs <- ortho.irf(res.std.ped,sigma.ind=c(348:351),beta.ind=c(4,3,5,6))
irfs.sum <- aaply(irfs,c(1,2,3),quantile,probs=c(0.05,0.5,0.95))

cor.cor <- irfs.sum[1,1,,]
cor.vol <- irfs.sum[2,1,,]
vol.cor <- irfs.sum[1,2,,]
vol.vol <- irfs.sum[2,2,,]

irf.full <- as.data.frame(rbind(cor.cor
                  ,cor.vol
                  ,vol.cor
                  ,vol.vol))
irf.full$Mode <- rep((c("Corruption on Corruption"
                       ,"Corruption on Volatility"
                       ,"Volatility on Corruption"
                       ,"Volatility on Volatility"))
                     ,each=7)
irf.full$Period <- rep(1:7,4)
names(irf.full)[1:3] <- c("Ten","Fifty","Ninety")

xyplot(Fifty~Period | Mode
       , data = irf.full
       , type="l"
       , ylab= "Impulse Response"
       , xlab = "Period"
       ,layout=c(2,2)
       ,xlim=c(1,7)
       #,ylim=c(-0.5,2)
       , index.cond = list(c(2,1,4,3))
       , subscripts=TRUE
       , panel=function(x,y,type,...,subscripts){
         panel.polygon(c(irf.full$Period[subscripts],rev(irf.full$Period[subscripts]))
                       ,c(irf.full$Ninety[subscripts],rev(irf.full$Ten[subscripts]))
                       ,col="gray65"
                       ,border="white")
         panel.xyplot(x,y,type,col="white",lwd=2)
         panel.abline(h=0,lty=2,col="gray45")
       }
       ,strip = function(..., bg,par.strip.text){
         strip.default(..., bg="black",par.strip.text=list(alpha=1
                                                           ,cex=1.05
                                                           ,col="white"
                                                           ,font=2
         ))
       }
)


##Inmediate effects
get.samples <- function(mcmc.li,XCor,XVol,coefs,sim.cor,sim.vol,iter=1e2,scales){
  len <- length(sim.vol)
  pred.vals <- array(NA,c(len,iter,2))
  for(i in 1:iter){
    chain.ind <- sample(length(mcmc.li),1)
    sample.ind <- sample(dim(mcmc.li[[chain.ind]])[1],1)
    betas <- mcmc.li[[chain.ind]][sample.ind,coefs]
    alphas <- mcmc.li[[chain.ind]][sample.ind,1:2]
    XCor.n <- matrix(rep(colMeans(XCor),len),ncol=2,byrow=TRUE)
    XCor.n[,1] <- sim.vol
    XVol.n <- matrix(rep(colMeans(XVol),len),ncol=2,byrow=TRUE)
    XVol.n[,1] <- sim.cor
    pred.vals[,i,1] <- (alphas[1]+XCor.n%*%betas[1:2])*scales[1] 
    pred.vals[,i,2] <- (alphas[2]+XVol.n%*%betas[3:4])*scales[2] 
  }
  return(pred.vals)
}

coef.interest <- grep("mu.beta",rownames(summary.res.std.ped$quantiles))

sims.vals.cor <- seq(min(std.ped.model$Y[,1]),max(std.ped.model$Y[,1]),length.out=100)
sims.vals.vol <- seq(min(std.ped.model$Y[,2]),max(std.ped.model$Y[,2]),length.out=100)
scale.vals <- attributes(std.ped.model$Y)$'scaled:scale'
pred.values <- get.samples(res.std.ped
                           ,std.ped.model$XCor
                           ,std.ped.model$XVol
                           ,coef.interest
                           ,sims.vals.cor
                           ,sims.vals.vol
                           ,scales=scale.vals)


ci.cor <- apply(pred.values[,,1],1,quantile,probs=c(0.05,0.5,0.95))
ci.cor <- as.data.frame(t(ci.cor))
ci.cor$Type <- "Predicted Corruption"
ci.cor$Sim <- sims.vals.vol*scale.vals[2]
ci.vol <- apply(pred.values[,,2],1,quantile,probs=c(0.05,0.5,0.95))
ci.vol <- as.data.frame(t(ci.vol))
ci.vol$Type <- "Predicted Volatility"
ci.vol$Sim <- sims.vals.cor*scale.vals[1]


full.effects <- rbind(ci.cor,ci.vol)
colnames(full.effects)[1:3] <- c("Five","Fifty","NinetyFive")

xyplot(Fifty~Sim | Type
       , data = full.effects
       , type="l"
       , ylab= "Predicted Value"
       , xlab = "Values of Other Variable"
       #, index.cond = list(c(2,1,4,3))
       ,scales = list(relation="free")
       ,ylim=list(c(1,5),c(0,100))
       , subscripts=TRUE
       , panel=function(x,y,type,...,subscripts){
         panel.polygon(c(full.effects$Sim[subscripts],rev(full.effects$Sim[subscripts]))
                       ,c(full.effects$NinetyFive[subscripts],rev(full.effects$Five[subscripts]))
                       ,col="gray65"
                       ,border="white")
         panel.xyplot(x,y,type,col="white",lwd=2)
       }
       ,strip = function(..., bg,par.strip.text){
         strip.default(..., bg="black",par.strip.text=list(alpha=1
                                                           ,cex=1.05
                                                           ,col="white"
                                                           ,font=2
         ))
       }
)