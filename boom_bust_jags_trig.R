model{
  
  #priors
  psinit ~ dbeta(spa, spb)
  
  # add a random year effect
  for(i in 1:(nyear - 1)){
    gy[i] ~ dnorm(g_mu, gy_tau)
    py[i] ~ dnorm(0, py_tau)
    
  }
  for(yr in 1:nyear){
    ly[yr] ~ dnorm(0, ly_tau)
  }
  
  # intercepts
  g_mu ~ dnorm(0, 0.20)
  p0 ~ dnorm(0, 0.20)
  lp ~ dnorm(0, 0.20)
  
  # fourier stuff
  a1_phi ~ dunif(0.01, 10)
  a1_gam ~ dunif(0.01, 10)
  dprobs_phi ~ dcat(probs_phi[1:P])
  dprobs_gam ~ dcat(probs_gam[1:P])
  for(i in 1:P){
    probs_phi[i] <- 1/P
    probs_gam[i] <- 1/P
  }
  
  a2_phi <- dprobs_phi - 1
  a2_gam <- dprobs_gam - 1
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
      logit(psi[k,t]) <- (z[k,t-1] * (py[t-1])) +
        ((z[k,t-1]) * (a1_phi * cos((pi*2*a2_phi)/2) * C[t-1] + a1_phi * sin((pi*2*a2_phi)/2) * S[t-1])) +
        ((1-z[k,t-1]) * gy[t-1])+
        ((1 - z[k,t-1]) * (a1_gam * cos((pi*2*a2_gam)/2) * C[t-1] + a1_gam * sin((pi*2*a2_gam)/2) * S[t-1]))
        
        
      
      z[k,t] ~ dbern(psi[k,t])
    }
  }
  
  # observational process
  
  for(j in 1:nsite){
    for(t in 1:nyear){
      logit(p[j,t]) <- lp + ly[t]
      y[j,t] ~ dbin(z[j,t]*p[j,t], jmat[j,t])
    }
  }
  
  
}