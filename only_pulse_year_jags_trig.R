model{
  
  #priors
  psinit ~ dbeta(spa, spb)
  
  # add a random year effect
  for(i in 1:(nyear - 1)){
  #  gy[i] ~ dnorm(g_mu, gy_tau)
     py[i] ~ dnorm(p_mu, py_tau)
  #  
  #}
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
  p1 ~ dnorm(0, 0.30)
  g1 ~ dnorm(0, 0.30)
  l1 ~ dnorm(0, 0.30)
 

  
  # fourier stuff
  theta ~ dgamma(1, 1)
  dprobs ~ dcat(probs[1:P])
  for(i in 1:P){
    probs[i] <- 1/P
  }
  
  dp <- dprobs - 1
  
  # half-cauchy hyperpriors for sd
  # of random effects
  #gy_sd ~ dt(0,25,1) T(0,)
  py_sd ~ dt(0,25,1) T(0,)
  ly_sd ~ dt(0,25,1) T(0,)
  ls_sd ~ dt(0,25,1) T(0,)
  
 # gy_tau <- 1/pow(gy_sd,2)
  py_tau <- 1/pow(py_sd,2)
  ly_tau <- 1/pow(ly_sd,2)
  ls_tau <- 1/pow(ls_sd,2)
  
  # latent state
  for(k in 1:nsite){
    z[k,1] ~ dbern(psinit)
    #z_pred[k,1] ~ dbern(psinit)
    for(t in 2:nyear){
      logit(psi[k,t]) <- (z[k,t-1] * (py[t-1])) +
      (z[k,t-1] * (p1 * cov[k])) + 
        ((1-z[k,t-1]) * g_mu)+
        ((1-z[k,t-1]) * g1 * cov[k])+
        ((1 - z[k,t-1]) * (theta * cos((pi*1*dp)/2) * C[t-1,1] + theta * sin((pi*1*dp)/2) * S[t-1,1]))+
        ((1 - z[k,t-1]) * (theta * cos((pi*2*dp)/2) * C[t-1,2] + theta * sin((pi*2*dp)/2) * S[t-1,2]))+
        ((1 - z[k,t-1]) * (theta * cos((pi*3*dp)/2) * C[t-1,3] + theta * sin((pi*3*dp)/2) * S[t-1,3]))+
        ((1 - z[k,t-1]) * (theta * cos((pi*4*dp)/2) * C[t-1,4] + theta * sin((pi*4*dp)/2) * S[t-1,4]))
      
      z[k,t] ~ dbern(psi[k,t])
      #z_pred[k,t] ~ dbern(psi[k,t])
    }
  }
  
  # observational process
  
  for(j in 1:nsite){
    for(t in 1:nyear){
      logit(p[j,t]) <- lp + ly[t] + l1 * cov[j] + ls[j]
      y[j,t] ~ dbin(z[j,t]*p[j,t], jmat[j,t])
      y_pred[j,t] ~ dbin(z[j,t]*p[j,t], jmat[j,t])
      #pz[j,t] <- z[j,t]*p[j,t]
    }
  }
  
  
}