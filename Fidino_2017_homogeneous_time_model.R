model{
#####################################
# Fidino 2017 Homogeneous time model#
#####################################
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
  Gt[i] ~ dnorm(d0, Gt_tau)
}
for(yr in 1:(nyear)){
  omega[yr] ~ dnorm(0, omega_tau)
}
for(site in 1:(nsite)){
  epsilon[site] ~ dnorm(0, epsilon_tau)
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
# half-cauchy hyperpriors for sd of random effects
# 
Gt_sd ~ dt(0,25,1) T(0,) # persistence
omega_sd ~ dt(0,25,1) T(0,) # detection time
epsilon_sd ~ dt(0,25,1) T(0,) # detection site
#
# convert sd to precision
#
Gt_tau <- 1/pow(Gt_sd,2) # persistence
omega_tau <- 1/pow(omega_sd,2) # detection time
epsilon_tau <- 1/pow(epsilon_sd,2) # detection site
#
####################
#latent state model#
####################
#
for(k in 1:(nsite)){
  z[k,1] ~ dbern(psinit) # initial occupancy
  for(t in 2:(nyear)){
    logit(psi[k,t]) <- (z[k,t-1] * (Gt[t-1] + d1 * x[k])) + 
      ((1-z[k,t-1]) * (m0 + m1 * x[k])) 
    z[k,t] ~ dbern(psi[k,t])
  }
}
#
#######################
#Observational process#
#######################
#
for(j in 1:(nsite)){
  for(t in 1:(nyear)){
    logit(p[j,t]) <- f0 + f1 * x[j] + omega[t] + epsilon[j] # detection
    y[j,t] ~ dbin(z[j,t]*p[j,t], jmat[j,t])
  }
}
}