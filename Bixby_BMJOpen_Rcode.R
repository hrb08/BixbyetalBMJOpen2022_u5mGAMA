# R Nimble code for estimating under-five mortality rates of neighbourhoods of the Greater Accra Metropolitan Area
# Bixby H, Bennett JE, Bawah AA, et al.
# Quantifying within-city inequalities in child mortality across neighbourhoods in Accra, Ghana: a Bayesian spatial analysisBMJ Open 2022;12:e054030

# The full microdata of the 2010 Ghana Population and Housing Census were accessed via the Ghana Statistical Service and are not publicly available. 
# Data from a random 10% sample of households enumerated in the 2010 Population and Housing census are publicly available and can be downloaded from 
# Ghana Statistical Services online data catalogue (https://www2.statsghana.gov.gh/nada/index.php/catalog/51).

# The shapefile used to generate the neighbourhood adjacency matrix ('shp_nhd.adj') was provided by the GIS team of the Ghana Statistical Service.

# The full microdata of the 2010 Ghana census included summary birth history data of all GAMA women aged 25-49 years.
# We summarised counts of child births and deaths summarised by neighbourhood and five-year age group of women and used the Maternal Age Cohort method (Rajaratnam et al 2012) to convert the summaries into probabilities of death under five (MAC 5q0). 
# Using this method each probability is assigned a reference year. Full details of the data pre-processing are provided in the Methods.
# The input data for the model was stored in the dataframe  'model_dat25plus'.

#-------------------------------
# constants
constants <- 
  list(
    N_obs = nrow(model_dat25plus), # the number of neighbourhood-year units.
    N_area = length(unique(model_dat25plus$nhd_idx)), # the number of neighbourhood areas.
    L = length(shp_nhd.adj$adj), # the number of adjacent neighbourhood pairs.
    adj = shp_nhd.adj$adj, # a vector giving the IDs of adjacent neighbourhoods for each neighbourhood area.
    weights = shp_nhd.adj$weights, # a vector of the same length as 'adj' giving weights associated with each pair of neighbourhood areas (equal to 1 for all pairs).
    num = shp_nhd.adj$num, # a vector of length N_area giving the number of adjacent neighbourhoods for each neighbourhood area.
    ceb = model_dat25plus$ceb, # a vector of length N_obs giving the number of children ever born for each five-year age group of women and neighbourhood.
    timeref =  model_dat25plus$timeref - mean(model_dat25plus$timeref), # a vector of length N_obs giving the reference year for each MAC 5q0 estimate (mean-centred).
    area = model_dat25plus$nhd_idx # a vector giving the numerical ID of each neighbourhood.
  )

# data
data <- list(y = model_dat25plus$fiveq0_probit) # a vector giving the MAC 5q0 estimates for each neighbourhood and reference year, transformed to the probit scale.

# ------BYM model code----------
# global terms: intercept (alpha) and linear slope (beta) 
# neighbourhood-specific terms: spatially structured intercept (U) and linear slope (Q) and spatially unstructured intercept (V) and linear slope (Z))
# sigma2_i is the neighbourhood-year specific variance term weighted by the number of births observed

code <- 
  nimbleCode({
    # Likelihood
    for(i in 1:N_obs) {
      y[i]~ dnorm( mean = mu[i],var = sigma2_i[i]) 
      mu[i] <- alpha + U[area[i]] + V[area[i]] + (beta + Q[area[i]] + Z[area[i]])*timeref[i]  
      sigma2_i[i] <- rho^2/delta[i]
      delta[i] <-  (ceb[i]+1) 
    }
    # Priors
    # global terms
    alpha ~ dflat()
    beta  ~ dnorm(0,tau =0.0001)
    rho ~  dinvgamma(2, rate = 2) 
    # neighbourhood-specific terms
    # spatially unstructured neighbourhood random intercept and slope
    for(j in 1:N_area) {
      V[j] ~ dnorm(0, tau =prec.v) 
      Z[j] ~ dnorm(0, tau =prec.z) 
    }
    prec.v ~ dgamma(0.5, rate = 0.0005)
    prec.z ~ dgamma(0.5, rate = 0.0005)
    #spatially structured neighbourhood random intercept and slope
    U[1:N_area] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N_area], prec.u, zero_mean = 1)
    Q[1:N_area] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N_area], prec.q, zero_mean = 1)
    prec.u ~ dgamma(0.5, rate = 0.0005)
    prec.q ~ dgamma(0.5, rate = 0.0005)
  })

# initial values
inits<- list(
  alpha = runif(1, 0.01,0.1),
  rho = rnorm(1, mean = 0.25, sd = 0.005), 
  prec.u = runif(1, 10,120),
  prec.v = runif(1, 10,120),
  prec.q = runif(1, 10,120),
  prec.z = runif(1, 10,120),
  U = rnorm(N_area, mean= 0.01,sd = 0.005),
  V = rnorm(N_area, mean= 0.01,sd = 0.005),
  Q = rnorm(N_area, mean= 0.1,sd = 0.05),
  Z = rnorm(N_area, mean= 0.1,sd = 0.05),
  beta = rnorm(1, mean=0.02, sd=0.001)
)

# setup and run model
nIts = 105000
nChn = 5
nBurn = 5000
nThin = 100
# construct model
Model <- nimbleModel(code = code, inits = inits, constants = constants, data = data)
# compile model
cModel <- compileNimble(Model)
# configure MCMC - set monitors
config <-  configureMCMC(Model, monitors = c( "mu","sigma2_i", names(inits)))
# build MCMC
MCMC <- buildMCMC(config,enableWAIC = TRUE)
# compile MCMC
cMCMC <- compileNimble(MCMC, project = Model, resetFunctions = TRUE)
# run model
samples <- runMCMC(cMCMC, nchains = nChn, niter = nIts,  nburnin = nBurn, thin = nThin, summary = TRUE, WAIC = TRUE, setSeed = c(1:nChn))

