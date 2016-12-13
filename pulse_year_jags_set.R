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
  
  # fourier stuff
  theta ~ dunif(0.01, 10)
  dprobs ~ dcat(probs[1:4])
  for(i in 1:4){
    probs[i] <- 0.25
  }
  
  dp <- dprobs - 1
  
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
      logit(psi[k,t]) <- (z[k,t-1] * (p0 + py[t-1])) +
        ((1-z[k,t-1]) * (g0 + gy[t-1]))+
        (1 - z[k,t-1] * (theta/(1*pi))*sin((pi*1)/4)*cos((pi*1 *((t-1)-dp))/2))+
        (1 - z[k,t-1] * (theta/(2*pi))*sin((pi*2)/4)*cos((pi*2 *((t-1)-dp))/2))+
        (1 - z[k,t-1] * (theta/(3*pi))*sin((pi*3)/4)*cos((pi*3 *((t-1)-dp))/2))
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