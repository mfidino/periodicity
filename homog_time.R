model{
  
  #priors
  psinit ~ dbeta(spa, spb)
  
  # add a random year effect
  for(i in 1:(nyear - 1)){
    py[i] ~ dnorm(p_mu, py_tau)
    
  }
  for(yr in 1:nyear){
    ly[yr] ~ dnorm(0, ly_tau)
  }
  for(site in 1:nsite){
    ls[site] ~ dnorm(0, ls_tau)
  }
  
  # intercepts
  g_mu ~ dnorm(0, 0.30)
  p_mu ~ dnorm(0, 0.30)
  lp ~ dnorm(0, 0.30)
  
  # covariate effect
  g1 ~ dnorm(0, 0.30)
  p1 ~ dnorm(0, 0.30)
  l1 ~ dnorm(0, 0.30)
  
  # half-cauchy hyperpriors for sd
  # of random effects
  py_sd ~ dt(0,25,1) T(0,)
  ly_sd ~ dt(0,25,1) T(0,)
  ls_sd ~ dt(0,25,1) T(0,)
  
  py_tau <- 1/pow(py_sd,2)
  ly_tau <- 1/pow(ly_sd,2)
  ls_tau <- 1/pow(ls_sd,2)
  
  # latent state
  for(k in 1:nsite){
    z[k,1] ~ dbern(psinit)
    z_pred[k,1] ~ dbern(psinit)
    for(t in 2:nyear){
      logit(psi[k,t]) <- (z[k,t-1] * (py[t-1] + p1 * cov[k])) + 
        ((1-z[k,t-1]) * (g_mu + (g1 * cov[k]))) 
      z[k,t] ~ dbern(psi[k,t])
      z_pred[k,t] ~ dbern(psi[k,t])
    }
  }
  
  # observational process
  
  for(j in 1:nsite){
    for(t in 1:nyear){
      logit(p[j,t]) <- lp + ly[t] + l1 * cov[j] + ls[j]
      y[j,t] ~ dbin(z[j,t]*p[j,t], jmat[j,t])
      y_pred[j,t] ~ dbin(z[j,t] * p[j,t], jmat[j,t])
    }
  }
  
  
}