setwd("~/ISPED/M2/STA305/Projet")

library(LaplacesDemon)
library(ggplot2)
library(bayesmeta)
library("metafor")
library(rjags)

# CREATION DES DONNEES SELON CELLES PRESENTES DANS L'ARTICLE
#======================================================================================
#======================================================================================

#Données des études sérologiques
data <- data.frame(
  region = c("Chelsea MA", "Gangelt, Germany", "Geneva, Switzerland (week 1)",
             "Geneva, Switzerland (week 2)", "Los Angeles County", "San Miguel County",
             "Santa Clara County"),
  test.supplier = c("BioMedics", "EUROIMMUN", "EUROIMMUN", "EUROIMMUN", "Premier Biotech",
                    "UBI Group", "Premier Biotech"),
  sample.size = c(200,500,343,417,863,4757,2583),
  number.positive = c(63,70,12,23,35,96,43),
  number.negative = c(137,430,331,440,828,4661,2540)
)

#Données des tests d'anticorps avec résultats faux positifs
data_test_FP <- data.frame(
  test.supplier = c("BioMedics", "EUROIMMUN", "Premier Biotech", "UBI Group"),
  sample.size = c(128,82,401,900),
  number.positive = c(12,3,2,0),
  number.negative = c(116,79,399,900)
)

#Données des tests d'anticorps avec résultats faux négatifs
data_test_FN <- data.frame(
  test.supplier = c("BioMedics", "EUROIMMUN", "Premier Biotech", "UBI Group"),
  sample.size = c(397,30,197,30),
  number.positive = c(352,20,178,30),
  number.negative = c(45,10,19,0)
)

#creation de variables
data$n_obs<-data$number.positive
data$N_obs<-data$sample.size

data_2<-data[, c(1,2,6,7)]

data_test_FP$n_fpos<-data_test_FP$number.positive
data_test_FP$N_neg<-data_test_FP$number.negative

data_test_FN$n_fneg<-data_test_FN$number.negative
data_test_FN$N_pos<-data_test_FN$number.positive

data_tot1<-merge(data_test_FP, data_test_FN, by="test.supplier")
data2<-data_tot1[,c(1,5,6,10,11)]

data_tot2<-merge(data2, data_2, by="test.supplier")

# save.image()
#======================================================================================
#======================================================================================


library(rjags)
jags_meta <- jags.model(file = "model_meta.txt", 
                             data = list(n_obs = data_tot2$n_obs,
                                         N_obs = data_tot2$N_obs,
                                         N_neg = data2$N_neg,
                                         N_pos = data2$N_pos,
                                         N1 = nrow(data_2), N2 = nrow(data_2)),
                             n.chains = 4)


library(rstan)
liste_stan = list(n_obs = data_tot2$n_obs,
                  N_obs = data_tot2$N_obs,
                  N_neg = data2$N_neg,
                  N_pos = data2$N_pos,
                  N1 = nrow(data_2), N2 = nrow(data_2))
fit <- stan("model_meta.stan",data = liste_stan)




model{
  
  
  n_obs ~ dbin(N_obs, p_obs)
    
  p_obs = p_prev * (1 - p_fneg) + (1 - p_prev) * p_fpos
  n_fpos ~ dbin(N_neg, p_fpos)
  n_fneg ~ dbin(N_pos, p_fneg) 
    
  
  #Priors
  p_prev <- 1/(1+exp(-alpha_prev))
    
  alpha_prev ~ dnorm(mu_prev, sigma_prev)
  
  
  p_fpos <- 1/(1+exp(-alpha_fpos))
  
  p_fneg <- 1/(1+exp(-alpha_fneg))
    
  alpha_fpos ~ dnorm(mu_fpos, sigma_fpos)
  alpha_fneg ~ dnorm(mu_fpos, sigma_fneg)

  
  #Hyper-prior
  mu_prev ~ dnorm(0,10)
  mu_fpos ~ dnorm(0,10)
  mu_fneg ~ dnorm(0,10)
  
  sigma_prev ~ dexp(1)
  sigma_fpos ~ dexp(1)
  sigma_fneg ~ dexp(1)
  
}

#==========================================================================

data <- data.frame(
  region = c("Chelsea MA", "Gangelt, Germany", "Geneva, Switzerland (week 1)",
             "Geneva, Switzerland (week 2)", "Los Angeles County", "San Miguel County",
             "Santa Clara County"),
  test.supplier = c("BioMedics", "EUROIMMUN", "EUROIMMUN", "EUROIMMUN", "Premier Biotech",
                    "UBI Group", "Premier Biotech"),
  sample.size = c(200,500,343,417,863,4757,2583),
  number.positive = c(63,70,12,23,35,96,43),
  number.negative = c(137,430,331,440,828,4661,2540)
)

data.es <- escalc(measure="PR", xi=number.positive, ni=sample.size, slab = region, data = data)
data.es
?escalc()

data_bayesmeta <- bayesmeta(y=data.es$yi, sigma=sqrt(data.es$vi), label=data.es$region,
                            tau.prior = function(t) {
                              dunif(t, max = 4)
                            }, mu.prior = c(0, 1), interval.type = "central")
summary(data_bayesmeta)
plot(data_bayesmeta, which=2)

#pfpos : FPR
#pprev

data.es.jags <- jags.model(file="model_seul.txt", data=list(logPR=data.es$yi,
                                                            sigma=sqrt(data.es$vi),
                                                            N=length(data.es$yi)),
                           n.chains=4)
# nobs,N1,N2,nfpos,nfneg pr ttes les locations et tous les tests

jags_res <- coda.samples(model = data.es.jags, variable.names = c("mu","tau"),
                         n.iter = 10000)

plot(jags_res)
summary(jags_res)


#================================================================================

# model{
#   
#   # Modèle d'échantillonnage/ Vraisemblance
#   for (i in 1:N1){
#     j=ifelse(i==1,1,ifelse(i<=4,2,ifelse(i==6,4,3))
#              n_nobs[i]~dbin(N_obs[i], p_obs[i])
#              p_obs[i] = p_prev[i]*(1-p_fneg[j]+(1-p_prev[i])*p_fpos[j]
#                                    n_fpos[j]~dbin(N_neg,p_fpos)
#                                    n_fneg[j] ~ dbin(N_pos,p_fneg)
#   } }
# 
# #Priors
# for (i in 1:N1){
#   p_prev[i] = exp(alpha_prev[i])/(1+exp(alpha_prev[i]))
#   alpha_prev[i]~dnorm(mu_prev, sigma_prev)
# } 
# for (j in 1:N2) {
#   alpha_fpos[j]~dnorm(mu_fpos, sigma_fpos)
#   alpha_fneg[j]~dnorm(mu_fneg, sigma_fneg) 
#   p_fpos[j] = exp(alpha_fpos)/(1+exp(alpha_fpos))
#   p_fneg[j] = exp(alpha_fneg)/(1+exp(alpha_fneg))
# }
# 
# #hyper-priors
# mu_prev ~ dnorm(0, 10)
# mu_fpos ~ dnorm(0, 10)
# mu_fneg ~ dnorm(0, 10)
# sigma_prev ~ dexp(1)
# sigma_fpos ~ dexp(1)
# sigma_fneg ~ dexp(1)
# }




data.es.jags <- jags.model(file="model_meta.txt", data=list(N_obs=data$N_obs,
                                                            n_obs=data$n_obs,
                                                            N1=7,
                                                            N2=4,
                                                            n_fpos=data_test_FP$n_fpos,
                                                            n_fneg=data_test_FN$n_fneg,
                                                            N_pos=data_test_FN$N_pos,
                                                            N_neg=data_test_FP$N_neg),
                           n.chains=4)

# nobs,N1,N2,nfpos,nfneg pr ttes les locations et tous les tests

jags_res <- coda.samples(model = data.es.jags, variable.names = c("p_prev"),
                         n.iter = 10000)

b <- window(jags_res, start=5001)
summary(b)

plot(jags_res)
summary(jags_res)
hist(unlist(jags_res[1][1][,1]), freq=F)
a <- jags_res[1][1][,1]
unlist(a)
# posterior = prior * vraisemblance 
#ordonnée pfpos du test j correspondant pour les densités plot multicolors


options(error=NULL)

#===========

model{
  
  # Modèle d'échantillonnage/ Vraisemblance
  
  n_obs[1]~dbin(N_obs[1], p_obs[1])
  n_obs[2]~dbin(N_obs[2], p_obs[2])
  n_obs[3]~dbin(N_obs[3], p_obs[3])
  n_obs[4]~dbin(N_obs[4], p_obs[4])
  n_obs[5]~dbin(N_obs[5], p_obs[5])
  n_obs[6]~dbin(N_obs[6], p_obs[6])
  n_obs[7]~dbin(N_obs[7], p_obs[7])
  
  p_obs[1] = p_prev[1]*(1-p_fneg[1])+(1-p_prev[1])*p_fpos[1]
  p_obs[2] = p_prev[2]*(1-p_fneg[2])+(1-p_prev[2])*p_fpos[2]
  p_obs[3] = p_prev[3]*(1-p_fneg[2])+(1-p_prev[3])*p_fpos[2]
  p_obs[4] = p_prev[4]*(1-p_fneg[2])+(1-p_prev[4])*p_fpos[2]
  p_obs[5] = p_prev[5]*(1-p_fneg[3])+(1-p_prev[5])*p_fpos[3]
  p_obs[6] = p_prev[6]*(1-p_fneg[4])+(1-p_prev[6])*p_fpos[4]
  p_obs[7] = p_prev[7]*(1-p_fneg[3])+(1-p_prev[7])*p_fpos[3]
  
  
  for (j in 1:N2){
    n_fpos[j]~dbin(N_neg,p_fpos)
    n_fneg[j] ~ dbin(N_pos,p_fneg)
  }
  
  #Priors
  for (i in 1:N1){
    p_prev[i] = exp(alpha_prev[i])/(1+exp(alpha_prev[i]))
    alpha_prev[i]~dnorm(mu_prev, sigma_prev)
  }
  for (j in 1:N2) {
    alpha_fpos[j]~dnorm(mu_fpos, sigma_fpos)
    alpha_fneg[j]~dnorm(mu_fneg, sigma_fneg) 
    p_fpos[j] = exp(alpha_fpos[j])/(1+exp(alpha_fpos[j]))
    p_fneg[j] = exp(alpha_fneg[j])/(1+exp(alpha_fneg[j]))
  }
  
  #hyper-priors
  mu_prev ~ dnorm(0, 10)
  mu_fpos ~ dnorm(0, 10)
  mu_fneg ~ dnorm(0, 10)
  sigma_prev ~ dexp(1)
  sigma_fpos ~ dexp(1)
  sigma_fneg ~ dexp(1)
}

mu_prev <- dnorm(0, 10)
mu_fpos <- dnorm(0, 10)
mu_fneg <- dnorm(0, 10)
sigma_prev <- dexp(1)
sigma_fpos <- dexp(1)
sigma_fneg <- dexp(1)


p_obs[1] = p_prev[1]*(1-p_fneg[1])+(1-p_prev[1])*p_fpos[1]

p_prev[i] = exp(alpha_prev[i])/(1+exp(alpha_prev[i]))

hist(rnorm(100, 0,10))


??dbin
