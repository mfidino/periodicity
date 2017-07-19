model{
####################################
# Fidino 2017 Stochastic time model#
####################################
#
########
#Priors#
########
#
# Initial occupancy
#
psinit ~ dbeta(spa, spb)
#
# Random effects
#
for(i in 1:(nyear - 1)){
  Ut[i] ~ dnorm(m0, Ut_tau) # colonization
  Gt[i] ~ dnorm(d0, Gt_tau) # persistence
}
for(yr in 1:(nyear)){
  omega[yr] ~ dnorm(0, omega_tau) # year effect detection
}
for(site in 1:(nsite)){
  epsilon[site] ~ dnorm(0, epsilon_tau) # site effect detection
}
#
# intercepts
#
m0 ~ dnorm(0, 0.30) # colonization
d0 ~ dnorm(0, 0.30) # persistence
f0 ~ dnorm(0, 0.30) # detection
#
# covariate effects
#
m1 ~ dnorm(0, 0.30) # colonization
d1 ~ dnorm(0, 0.30) # persistence
f1 ~ dnorm(0, 0.30) # detection
#
# half-Cauchy hyperpriors for sd of random effects
# 
Ut_sd ~ dt(0,25,1) T(0,) # colonization
Gt_sd ~ dt(0,25,1) T(0,) # persistence
omega_sd ~ dt(0,25,1) T(0,) # detection time
epsilon_sd ~ dt(0,25,1) T(0,) # detection site
#
# convert sd to precision
#
Ut_tau <- 1/pow(Ut_sd,2) # colonization
Gt_tau <- 1/pow(Gt_sd,2) # persistence
omega_tau <- 1/pow(omega_sd,2) # detection time
epsilon_tau <- 1/pow(epsilon_sd,2) # detection site
#
####################
#Latent state model#
####################
#
for(k in 1:(nsite)){
  z[k,1] ~ dbern(psinit) 
  for(t in 2:(nyear)){
    logit(psi[k,t]) <- (z[k,t-1] * (Gt[t-1] + d1 * x[k])) + 
       ((1-z[k,t-1]) * (Ut[t-1] + m1 * x[k])) 
    z[k,t] ~ dbern(psi[k,t]) # likelihood
  }
}
#
#######################
#Observational process#
#######################
#
for(k in 1:(nsite)){
  for(t in 1:(nyear)){
    logit(p[k,t]) <- f0 + f1 * x[k] + omega[t] + epsilon[k]
    y[k,t] ~ dbin(z[k,t]*p[k,t], jmat[k,t])
  }
}
}