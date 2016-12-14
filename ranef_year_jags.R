model{
  
  #priors
  psinit ~ dbeta(spa, spb)
  
  # add a random year effect
  for(i in 1:(nyear - 1)){
    gy[i] ~ dnorm(0, gy_tau)
    py[i] ~ dnorm(0, py_tau)
    
  }
  for(yr in 1:nyear){
    ly[yr] ~ dnorm(0, ly_tau)
  }
  
  # intercepts
  g0 ~ dnorm(0, 0.20)
  p0 ~ dnorm(0, 0.20)
  lp ~ dnorm(0, 0.20)
  
  # covariate effect
  g1 ~ dnorm(0, 0.20)
  p1 ~ dnorm(0, 0.20)
  
  # half-cauchy hyperpriors for sd
  # of random effects
  gy_sd ~ dt(0,25,1) T(0,)
  py_sd ~ dt(0,25,1) T(0,)
  ly_sd ~ dt(0,25,1) T(0,)
  
  gy_tau <- 1/pow(gy_sd,2)
  py_tau <- 1/pow(py_sd,2)
  ly_tau <- 1/pow(ly_sd,2)
  
  # latent state
  for(k in 1:nsite){
    z[k,1] ~ dbern(psinit)
    for(t in 2:nyear){
      logit(psi[k,t]) <- (z[k,t-1] * (p0 + py[t-1]) + p1 * cov[k]) + 
        ((1-z[k,t-1]) * (g0 + gy[t-1]) + g1 * cov[k])
      z[k,t] ~ dbern(psi[k,t])
    }
  }
  
  # observational process
  
  for(j in 1:nsite){
    for(t in 1:nyear){
      logit(p[j,t]) <- lp + ly[t]
      y[j,t] ~ dbin(z[j,t]*p[j,t], jmat[j,t])
      y_pred[j,t] ~ dbin(z[j,t] * p[j,t], jmat[j,t])
    }
  }
  
  
}