model{
##############################################
# Fidino 2017 Periodic time model: boom bust #
##############################################
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
  Gt[i] ~ dnorm(d0, Gt_tau)
}
for(yr in 1:(nyear)){
  omega[yr] ~ dnorm(0, omega_tau)
}
for(site in 1:(nsite)){
  epsilon[site] ~ dnorm(0, epsilon_tau)
}
#
# Intercepts
#
m0 ~ dnorm(0, 0.30) # colonization
d0 ~ dnorm(0, 0.30) # persistence
f0 ~ dnorm(0, 0.30) # detection
#
# Covariate effects
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
Ut_tau <- 1/pow(Gt_sd,2) # colonization
Gt_tau <- 1/pow(Gt_sd,2) # persistence
omega_tau <- 1/pow(omega_sd,2) # detection time
epsilon_tau <- 1/pow(epsilon_sd,2) # detection site
#
# Fourier priors, delta supplied as data
#
A ~ dgamma(1, 1) # Amplitude
#
# Transform A and delta to B1 and B2
#
B1 <- A * cos((pi*2*delta)/2)
B2 <- A * sin((pi*2*delta)/2)
#
####################
#latent state model#
####################
#
for(k in 1:(nsite)){
  z[k,1] ~ dbern(psinit)
  for(t in 2:nyear){
    logit(psi[k,t]) <- (z[k,t-1] * (Gt[t-1] + d1 * x[k])) + # persistence
      ((1-z[k,t-1]) * m0 + m1 * x[k] + B1 * C[t-1] + B2 * S[t-1]) # colonization
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