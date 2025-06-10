# PACKAGES
library(dplyr)
library(tidyr)
library(readr)
library(brms)
library(tidyverse)
library(ggplot2)
library(mgcv)
library(rjags)

#DATA
data <- read_csv("Outputs/Combined_data.csv")

# OPTIONAL STEP TO EXCLUDE OUTLIERS HERE
data <- data %>%
  mutate(
    LMA_leaf = ifelse(LMA_leaf_QC_dist > 0, NA_real_, LMA_leaf),
    leaf_thickness = ifelse(thickness_QC_dist > 0, NA_real_, leaf_thickness),
    leaf_area = ifelse(leaf_area_QC_dist > 0, NA_real_, leaf_area),
    leaf_drymass = ifelse(leaf_drymass_QC_dist > 0, NA_real_, leaf_drymass),
    Narea = ifelse(Narea_QC_dist > 0, NA_real_, Narea),
    Nmass = ifelse(Nmass_QC_dist > 0, NA_real_, Nmass),
    Cmass = ifelse(Cmass_QC_dist > 0, NA_real_, Cmass),
    PTLP = ifelse(PTLP_QC_dist > 0, NA_real_, PTLP),
    Vcmax_Estimate_final = ifelse(Vcmax_Estimate_final_QC_dist > 0, NA_real_, Vcmax_Estimate_final),
    Jmax_Estimate_final = ifelse(Jmax_Estimate_final_QC_dist > 0, NA_real_, Jmax_Estimate_final),
    TPU_Estimate_final = ifelse(TPU_Estimate_final_QC_dist > 0, NA_real_, TPU_Estimate_final)
  )

species <- unique(data$sp)

data_1sp <- data %>%
  filter(sp == species[2]) %>%
  subset(!is.na(Age_leaf_mid))

summary(data_1sp[, "Age_leaf_range"])

#do the age is uniformly distributed between age min and age max, or normal around age mid ?

model1 <- "model{


    # Likelihood model - censored in Xobs
    for(k in 1:N){
    
  yobs[k] ~ dnorm(y[k], tau_y)
	y[k] <- a + b*x[k]
  
  x[k] ~ dunif(0, lifespan)
  
	p[k] <- pnorm(cut[k], x[k], tau_x)  
	#?? quelle erreur de mesure, normale avec sigma autour d'1 jour to start
	#cut[k] in data c(age_min[k], age_max[k])
	
  z[k] ~ dbern(p[k]) #in data z = 1 because interval censored
        
    }

    # Priors
    tau_x ~ dgamma(1, 0.001)
    tau_y ~ dgamma(1, 0.001)
         
    a ~ dnorm(0, 1/pow(sd, 2))
    b ~ dnorm(1, 1/pow(10, 2))

}"



inits1 <- list(list(sigma = runif(1, 0, 5),
                   alpha0 = rnorm(1,mean = 0, sd = 1),
                   alpha1 = rnorm(1,mean = 0, sd = 1),
                   alpha2 = rnorm(1,mean = 0, sd = 1)),
               list(sigma =  runif(1, 0, 5),
                    alpha0 = rnorm(1,mean = 0, sd = 1),
                    alpha1 = rnorm(1,mean = 0, sd = 1),
                    alpha2 = rnorm(1,mean = 0, sd = 1)),
               list(sigma =  runif(1, 0, 5),
                    alpha0 = rnorm(1,mean = 0, sd = 1),
                    alpha1 = rnorm(1,mean = 0, sd = 1),
                    alpha2 = rnorm(1,mean = 0, sd = 1))
)

data_jags1 <- list(N = nrow(data_1sp),
                  y_trait = log(data_1sp$Vcmax_Estimate_final) ,
                  age_min = data_1sp$Age_leaf_min,
                  age_max = data_1sp$Age_leaf_max)

m1 <- jags.model(textConnection(model1), inits = inits1, data = data_jags1, n.chains = 3)
update(m1, 2000)
mcmc1 <- coda.samples(m1, variable.names = c("sigma", "alpha0", "alpha1", "alpha2"), 
                             n.iter = 50000, thin = 50)

