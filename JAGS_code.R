#Kayla Blincow
#3/3/2021

#The purpose of the script is to manually construct the JAGS mixing model to test
#for the relative contributions of different primary production sources to GSB
#diets by size.

#Much of this code draws from the functions and model output from the 
#MixSIAR package created by Brice Semmens, Brian Stock, Andrew Jackson, 
#Eric Ward, Andrew Parnell, and Donald Phillips

#clear my workspace
rm(list = ls())

#load packages
library(tidyverse)
library(R2jags)
library(MCMCpack)


#load necessary data
source <- read.csv("PPsources.csv", header = T)
mix <- read.csv("FinalGSBBulk.csv", header = T)
discr <- read.csv("gsb_TEF.csv", header = T)

#####create objects needed to run the model####

#***INFO NEEDED FROM MIXTURE (GSB) DATA***
#actual isotope values for GSBs
X_iso <- mix[,c("d13C", "d15N")] 

#number of GSB samples
N <- nrow(mix)

#number of factor levels (for random effect of individual)
factor1_levels <- nrow(mix)

#unique factor values for each individual
Factor.1 <- as.numeric(factor(mix$tag_ID))



#***INFO NEEDED FROM SOURCE DATA***
#number of PP sources (2: phytoplankton, kelp)
n.sources <- nrow(source) 

#array of mean source isotope values
MU_array <- array(data = c(source$Meand13C, source$Meand15N),
                  dim = c(2,2))

#array of total number of samples used to calculate the mean source values
n_array <- source$n

#array of standard deviation of source isotope values
SIG_array <- array(data = c(source$SDd13C, source$SDd15N),
                   dim = c(2,2))

#convert that to variance
SIG2_array <- SIG_array^2 

#***INFO NEEDED FROM TROPHIC DISCRIMINATION FACTOR DATA***
#mean fractionation values
frac_mu <- discr[, c("Meand13C", "Meand15N")]

#sd fractionation values
frac_sig2 <- discr[, c("SDd13C", "SDd15N")] #sd fractionation values


#***EXTRA INFO NEEDED***
#number of isotopes informing model (2: d15N, d13C)
n.iso <- 2 

#calculate pooled mean and standard deviation of isotope values and standardize
for(j in 1:n.iso){
  mean.pool <- (N * mean(X_iso[, j], na.rm = T) + 
                  as.vector(as.vector(n_array) %*% 
                              as.vector(MU_array[,j])))/sum(c(as.vector(n_array), N))
  sd.pool <- sqrt((sum((as.vector(n_array) - 1) * as.vector(SIG2_array[, j])) + 
                     as.vector(as.vector(n_array) %*% 
                                 as.vector(MU_array[, j])^2) + (N - 1) * 
                     stats::var(X_iso[,j], na.rm = T) + N * 
                     mean(X_iso[, j], na.rm = T)^2 - 
                     sum(c(as.vector(n_array), N)) * 
                     mean.pool^2)/(sum(c(as.vector(n_array), N)) - 1))
  
  MU_array[, j] <- (MU_array[, j] - mean.pool)/sd.pool
  SIG2_array[, j] <- SIG2_array[, j]/sd.pool^2
  
  X_iso[, j] <- (X_iso[, j] - mean.pool)/sd.pool
  frac_mu[, j] <- frac_mu[, j]/sd.pool
  frac_sig2[, j] <- frac_sig2[, j]/sd.pool^2
}


#set alpha for specifying initial dirichlet prior conditions for p.global
alpha <- rep(1, n.sources) 

#e is a value needed to calculate the isometric log ratio
#(Aitchison-orthonormal basis)
e <- matrix(rep(0, n.sources * (n.sources - 1)), nrow = n.sources, 
            ncol = (n.sources - 1))
for (i in 1:(n.sources - 1)) {
  e[, i] <- exp(c(rep(sqrt(1/(i * (i + 1))), i), -sqrt(i/(i + 1)), 
                  rep(0, n.sources - i - 1)))
  e[, i] <- e[, i]/sum(e[, i])
}

#create necessary empty arrays which will be populated when calculating the ILR
cross <- array(data = NA, dim = c(N, n.sources, n.sources - 1))
tmp.p <- array(data = NA, dim = c(N, n.sources))
cross.fac1 <- array(data = NA, dim = c(factor1_levels, n.sources, n.sources-1))
tmp.p.fac1 <- array(data = NA, dim = c(factor1_levels, n.sources))

#create objects of compiled data
all.data <- c("X_iso", "N", "n.sources", 
              "n.iso", "alpha", "frac_mu", "e", 
              "cross", "tmp.p")
f.data <- c("factor1_levels", "Factor.1", 
            "cross.fac1", "tmp.p.fac1")
s.data <- c("MU_array", "SIG2_array", "n_array")

#specify continuous effect (normalized length data for GSBs)
assign("Cont.1", as.vector(scale(mix$TotalLength)))
c.data <- "Cont.1"


####Specify and Run the Model####
#Specify the model
sink("GSB_MixingModel.txt")
cat(
  "model{
    for(src in 1:n.sources){
      for(iso in 1:n.iso){
        src_mu[src,iso] ~ dnorm(MU_array[src,iso], n_array[src]/SIG2_array[src,iso]);  
        tmp.X[src,iso] ~ dchisqr(n_array[src]);
        src_tau[src,iso] <- tmp.X[src,iso]/(SIG2_array[src,iso]*(n_array[src] - 1));   
      }
    }
    
    # draw p.global (global proportion means) from an uninformative Dirichlet,
    # then ilr.global is the Isometric Log Ratio-transform of p.global
    # ILR Equation from page 206, Egozcue 2003
    p.global[1:n.sources] ~ ddirch(alpha[1:n.sources]);
    for(src in 1:(n.sources-1)){
      gmean[src] <- prod(p.global[1:src])^(1/src);
      ilr.global[src] <- sqrt(src/(src+1))*log(gmean[src]/p.global[src+1]);
    }
    
    fac1.sig ~ dunif(0,20);
    fac1.invSig2 <- 1/(fac1.sig*fac1.sig);
    # draw the fac1 (region) specific ILR terms (random effect)
    for(f1 in 1:factor1_levels) {
      for(src in 1:(n.sources-1)) {
        ilr.fac1[f1,src] ~ dnorm(0,fac1.invSig2);
      }
    }
    
    # Prior on continuous effects (slopes for a linear regression in ilr-space)
      ilr.cont1 ~ dnorm(0,.001)
    
    # DON'T generate individual deviates from the global/region/pack mean (but keep same model structure)
    for(i in 1:N) {
      for(src in 1:(n.sources-1)) {
        ilr.ind[i,src] <- 0;
        ilr.tot[i,src] <- ilr.global[src] + ilr.fac1[Factor.1[i],src] + ilr.cont1[src]*Cont.1[i] + ilr.ind[i,src]; # add all effects together for each individual (in ilr-space)
      }
    }
    
    # Inverse ILR math (equation 24, page 294, Egozcue 2003)
    for(i in 1:N){
      for(j in 1:(n.sources-1)){
        cross[i,,j] <- (e[,j]^ilr.tot[i,j])/sum(e[,j]^ilr.tot[i,j]);
      }
      for(src in 1:n.sources){
        tmp.p[i,src] <- prod(cross[i,src,]);
      }
      for(src in 1:n.sources){
        p.ind[i,src] <- tmp.p[i,src]/sum(tmp.p[i,]);
      }
    }
    
    for(src in 1:n.sources) {
      for(i in 1:N){
        # these are weights for variances
        p2[i,src] <- p.ind[i,src]*p.ind[i,src];
      }
    }
    
    # Transform ilr.fac1 into p.fac1 
    for(f1 in 1:factor1_levels) {
      for(src in 1:(n.sources-1)) {
        ilr.fac1.tot[f1,src] <- ilr.global[src] + ilr.fac1[f1,src];
        cross.fac1[f1,,src] <- (e[,src]^ilr.fac1.tot[f1,src])/sum(e[,src]^ilr.fac1.tot[f1,src]);
      }
      for(src in 1:n.sources) {
        tmp.p.fac1[f1,src] <- prod(cross.fac1[f1,src,]);
      }
      for(src in 1:n.sources){
        p.fac1[f1,src] <- tmp.p.fac1[f1,src]/sum(tmp.p.fac1[f1,]);
      }
    }
    
    
    # for each isotope and population, calculate the predicted mixtures
    for(iso in 1:n.iso) {
      for(i in 1:N) {
        
        mix.mu[iso,i] <- inprod(src_mu[,iso],p.ind[i,]) + inprod(frac_mu[,iso],p.ind[i,]);
      }
    }
    
    # Multiplicative residual error
    for(iso in 1:n.iso){
      resid.prop[iso] ~ dunif(0,20);
    }
    
    
    # Calculate process variance for each isotope and population
    for(iso in 1:n.iso) {
      for(i in 1:N) {
        
        process.var[iso,i] <- inprod(1/src_tau[,iso],p2[i,]) + inprod(frac_sig2[,iso],p2[i,]);
      }
    }
    
    # Construct Sigma, the mixture precision matrix
    for(ind in 1:N){
      for(i in 1:n.iso){
        for(j in 1:n.iso){
          Sigma.ind[ind,i,j] <- equals(i,j)/(process.var[i,ind]*resid.prop[i]);
        }
      }
    }
    
    # Likelihood
    for(i in 1:N) {
      X_iso[i,] ~ dmnorm(mix.mu[,i], Sigma.ind[i,,]);
      loglik[i] <- logdensity.mnorm(X_iso[i,], mix.mu[,i], Sigma.ind[i,,]);
    }
  } # end model
  "
  , fill = TRUE)
sink()

#Bundle data
jags.data <- c(all.data, f.data, s.data, c.data, "frac_sig2")

#Initial values
inits <- function() list(p.global = as.vector(rdirichlet(1, alpha))) 

#Parameters I want to track
jags.params <- c("p.global", "loglik", "p.fac1", "ilr.fac1", "fac1.sig", 
                 "ilr.global", "ilr.cont1", "p.ind", "resid.prop")

#MCMC Settings
ni <- 50000 #iterations
nt <- 5 #thinning by
nb <- 25000 #burn in
nc <- 3 #chains

mod <- jags(jags.data, inits = inits, parameters.to.save = jags.params,
     model.file = "GSB_MixingModel.txt", n.chains = nc, n.burnin = nb, n.thin = nt,
     n.iter = ni)

####posterior checks####

#check that Rhats and n.effs are acceptable
print(mod)

#check the traceplots of parameters for convergence
traceplot(mod)

#attach the model to look at specific parameter values
attach.jags(mod)


####plot the results in terms of proportion of PP source by length####

#create an object with the source names
source_names <- source$Source

#create functions to return upper and lower quantiles for 95% credible intervals
get_high <- function(x){return(quantile(x,.975))}
get_low <- function(x){return(quantile(x,.025))}

#set the number of x values to plot (i.e. length values)
n.plot = 200

#specify the chain length of estimates for transforming from ILR space
chain.len = dim(p.global)[1]

#create a sequence of standardized length values which will form the x axis
Cont1.plot <- seq(from = round(min(scale(mix$TotalLength)), 1), 
                  to = round(max(scale(mix$TotalLength)), 1), 
                  length.out = n.plot)

#create empty arrays that will be populated with ilr estimates
ilr.plot <- array(NA, dim = c(n.plot, n.sources-1, chain.len))
ilr.median <- array(NA, dim = c(n.plot, n.sources-1))
ilr.low <- array(NA, dim = c(n.plot, n.sources-1))
ilr.high <- array(NA, dim = c(n.plot, n.sources-1))

#fill those arrays with appropriate model estimates
for(src in 1:n.sources-1){
  for(i in 1:n.plot){
    ilr.plot[i,src,] <- ilr.global[,src] + ilr.cont1[,src]*Cont1.plot[i]
    ilr.low[i,src] <- get_low(ilr.plot[i,src,])
    ilr.median[i,src] <- median(ilr.plot[i,src,])
    ilr.high[i,src] <- get_high(ilr.plot[i,src,])
  }
}

# Transform regression lines from ILR-space to p-space
e <- matrix(rep(0,n.sources*(n.sources-1)),
            nrow = n.sources, ncol = (n.sources-1))

for(i in 1:(n.sources - 1)){
  e[,i] <- exp(c(rep(sqrt(1/(i*(i+1))),i),-sqrt(i/(i+1)),rep(0,n.sources-i-1)))
  e[,i] <- e[,i]/sum(e[,i])
}

#create a bunch of empty arrays of dummy variables for the inverse ILR calculation
cross.med <- array(data = NA, dim = c(n.plot, n.sources, n.sources-1))

tmp.p.med <- array(data = NA, dim = c(n.plot, n.sources))

p.median <- array(data = NA, dim = c(n.plot, n.sources))

cross.low <- array(data = NA, dim = c(n.plot, n.sources, n.sources-1))

tmp.p.low <- array(data = NA, dim = c(n.plot, n.sources))

p.low <- array(data = NA, dim = c(n.plot, n.sources))

cross.high <- array(data = NA, dim = c(n.plot, n.sources, n.sources-1))

tmp.p.high <- array(data = NA, dim = c(n.plot, n.sources))

p.high <- array(data = NA, dim = c(n.plot, n.sources))

  
#do the inverser ILR calculation
for(i in 1:n.plot){
  for(j in 1:(n.sources-1)){
    cross.med[i,,j] <- (e[,j]^ilr.median[i,j])/sum(e[,j]^ilr.median[i,j]);
    cross.low[i,,j] <- (e[,j]^ilr.low[i,j])/sum(e[,j]^ilr.low[i,j]);
    cross.high[i,,j] <- (e[,j]^ilr.high[i,j])/sum(e[,j]^ilr.high[i,j]);
  }
  for(src in 1:n.sources){
    tmp.p.med[i,src] <- prod(cross.med[i,src,]);
    tmp.p.low[i,src] <- prod(cross.low[i,src,]);
    tmp.p.high[i,src] <- prod(cross.high[i,src,]);
  }
  for(src in 1:n.sources){
    p.median[i,src] <- tmp.p.med[i,src]/sum(tmp.p.med[i,]);
    p.low[i,src] <- tmp.p.low[i,src]/sum(tmp.p.low[i,]);
    p.high[i,src] <- tmp.p.high[i,src]/sum(tmp.p.high[i,]);
  }
}

#assign source names to appropriate median proportion columns
colnames(p.median) <- source_names

#convert x axis back to the original scale 
Cont1.plot <- Cont1.plot*sd(mix$TotalLength) + mean(mix$TotalLength)

#create a dataframe with all our new data for plotting
df <- data.frame(reshape2::melt(p.median)[,2:3], rep(Cont1.plot,n.sources), 
                 reshape2::melt(p.low)[,3], reshape2::melt(p.high)[,3])
colnames(df) <- c("source","median","x","low","high")

# Plot of Diet Proportion vs. Length
png("diet_length.png", height=7, width=7, units='in', res=500)

print(ggplot(data = df, aes(x = x, y = median)) +
        geom_line(aes(x = x, y = median,
                                        group = source, 
                                        colour = source), size=1.5) +
        geom_ribbon(ggplot2::aes(ymin = low, ymax = high, 
                                          group = source, fill = source), 
                             alpha=0.35) +
        ylab("Diet Proportion") +
        xlab("Total Length (cm)") +
        scale_y_continuous(expand = c(0, 0), limits=c(0,1)) +
        theme_bw() +
        theme(panel.border = ggplot2::element_blank(), 
                       panel.grid.major = ggplot2::element_blank(), 
                       panel.grid.minor = ggplot2::element_blank(), 
                       panel.background = ggplot2::element_blank(), 
                       axis.line = ggplot2::element_line(colour = "black"), 
                       axis.title = ggplot2::element_text(size = 16), 
                       axis.text = ggplot2::element_text(size = 14), 
                       legend.text = ggplot2::element_text(size = 14), 
                       legend.position = c(.7,.5), 
                       legend.justification = c(0,1), 
                       legend.title = ggplot2::element_blank()))
dev.off()


