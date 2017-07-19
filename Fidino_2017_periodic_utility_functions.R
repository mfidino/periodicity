############################
#
# Utility script for Fidino and Magle 2017:
# Using Fourier series to estimate periodic
# trends in dynamic occupancy models.
#
#

########################################################################
# Copyright 2017 Mason Fidino

# Permission is hereby granted, free of charge, 
# to any person obtaining a copy of this software 
# and associated documentation files, (the "Software")
# to deal in the Software without restriction, including 
# without limitation the rights to use, copy, modify, merge, 
# publish, distribute, sublicense, and/or sell copies of the 
# Software, and to permit persons to whom the Software is 
# furnished to do so, subject to the following conditions:
  
# The above copyright notice and this permission notice 
# shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY 
# KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE 
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE 
# AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
# IN THE SOFTWARE.
########################################################################
#
#
### DESCRIPTION ###### DESCRIPTION ###### DESCRIPTION ###### DESCRIPTION ###
# This script contains utility functions to fit homogeneous time,
# stochastic time, and periodic time dynamic occupancy models
### DESCRIPTION ###### DESCRIPTION ###### DESCRIPTION ###### DESCRIPTION ###

# package_load:
# A general function to load packages, and if not on a computer
# to download them first
package_load<-function(packages = NA, quiet=TRUE, verbose=FALSE, warn.conflicts=FALSE){
  
  # download required packages if they're not already
  pkgsToDownload<- packages[!(packages  %in% installed.packages()[,"Package"])]
  if(length(pkgsToDownload)>0)
    install.packages(pkgsToDownload, repos="http://cran.us.r-project.org", quiet=quiet, verbose=verbose)
  
  # then load them
  for(i in 1:length(packages))
    require(packages[i], character.only=T, quietly=quiet, warn.conflicts=warn.conflicts)
}

# df_2_array:
# Converts a dataframe to a three dimensional array.
# needed to convert our camera trap data
# to an array for use in JAGS
df_2_array <- function(my_df = NULL){
  
  require(reshape2)
  my_array <- acast(my_df, species~site~season, value.var = "count")
  dimnames(my_array) <- NULL
  return(my_array)
}

# betaABfromMeanSD:
# Calculates shape parameters of a Beta distribution
# from a mean and standard deviation
betaABfromMeanSD = function( mean , sd ) {
  if ( mean <=0 | mean >= 1) stop("must have 0 < mean < 1")
  if ( sd <= 0 ) stop("sd must be > 0")
  kappa = mean*(1-mean)/sd^2 - 1
  if ( kappa <= 0 ) stop("invalid combination of mean and sd")
  a = mean * kappa
  b = ( 1.0 - mean ) * kappa
  return( list( a=a , b=b ) )
}

# inits_STM:
# Generates initial values of the stochastic time model
# works for up to 8 chains
inits_STM <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = z,
      psinit = rbeta(1,1,1), # initial occupancy
      m0 = runif(1, -3, 3), # colonization
      m1 = runif(1, -3, 3),
      Ut = runif(8, -3, 3),
      Ut_sd = rgamma(1,1,1),
      d0 = runif(1, -3, 3), #persistence
      d1 = runif(1, -3, 3),
      Gt = runif(8, -3, 3),
      Gt_sd = rgamma(1,1,1),
      f0 = runif(1, -3, 3), # detection
      f1 = runif(1, -3, 3),
      omega = runif(9, -3, 3),
      omega_sd = rgamma(1,1,1),
      epsilon = runif(95, -3, 3),
      epsilon_sd = rgamma(1,1,1),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Marsaglia-Multicarry",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(switch(chain,           
                "1" = gen_list(chain),
                "2" = gen_list(chain),
                "3" = gen_list(chain),
                "4" = gen_list(chain),
                "5" = gen_list(chain),
                "6" = gen_list(chain),
                "7" = gen_list(chain),
                "8" = gen_list(chain)
  )
  )
}

# inits_HTM:
# Generates initial values of the homogeneous time model
# works for up to 8 chains
inits_HTM <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = z,
      psinit = rbeta(1,1,1), # initial occupancy
      m0 = runif(1, -3, 3), # colonization
      m1 = runif(1, -3, 3),
      d0 = runif(1, -3, 3), #persistence
      d1 = runif(1, -3, 3),
      Gt = runif(8, -3, 3),
      Gt_sd = rgamma(1,1,1),
      f0 = runif(1, -3, 3), # detection
      f1 = runif(1, -3, 3),
      omega = runif(9, -3, 3),
      omega_sd = rgamma(1,1,1),
      epsilon = runif(95, -3, 3),
      epsilon_sd = rgamma(1,1,1),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Marsaglia-Multicarry",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(switch(chain,           
                "1" = gen_list(chain),
                "2" = gen_list(chain),
                "3" = gen_list(chain),
                "4" = gen_list(chain),
                "5" = gen_list(chain),
                "6" = gen_list(chain),
                "7" = gen_list(chain),
                "8" = gen_list(chain)
  )
  )
}

# inits_PTM_pulse:
# Generates initial values of the periodic time model
# for a single season pulse
# works for up to 8 chains
inits_PTM_pulse <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = z,
      psinit = rbeta(1,1,1), # initial occupancy
      m0 = runif(1, -3, 3), # colonization
      m1 = runif(1, -3, 3),
      d0 = runif(1, -3, 3), #persistence
      d1 = runif(1, -3, 3),
      Gt = runif(8, -3, 3),
      Gt_sd = rgamma(1,1,1),
      f0 = runif(1, -3, 3), # detection
      f1 = runif(1, -3, 3),
      omega = runif(9, -3, 3),
      omega_sd = rgamma(1,1,1),
      epsilon = runif(95, -3, 3),
      epsilon_sd = rgamma(1,1,1),
      A = rgamma(1,1,1),
      delta_plus_one = sample(1:4,1),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Marsaglia-Multicarry",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(switch(chain,           
                "1" = gen_list(chain),
                "2" = gen_list(chain),
                "3" = gen_list(chain),
                "4" = gen_list(chain),
                "5" = gen_list(chain),
                "6" = gen_list(chain),
                "7" = gen_list(chain),
                "8" = gen_list(chain)
  )
  )
}
# inits_PTM_boom:
# Generates initial values of the periodic time model
# for boom bust pattern
# works for up to 8 chains

inits_PTM_boom <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = z,
      psinit = rbeta(1,1,1), # initial occupancy
      m0 = runif(1, -3, 3), # colonization
      m1 = runif(1, -3, 3),
      d0 = runif(1, -3, 3), #persistence
      d1 = runif(1, -3, 3),
      Gt = runif(8, -3, 3),
      Gt_sd = rgamma(1,1,1),
      f0 = runif(1, -3, 3), # detection
      f1 = runif(1, -3, 3),
      omega = runif(9, -3, 3),
      omega_sd = rgamma(1,1,1),
      epsilon = runif(95, -3, 3),
      epsilon_sd = rgamma(1,1,1),
      A = rgamma(1,1,1),
      delta_plus_one = sample(1:2,1),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Marsaglia-Multicarry",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(switch(chain,           
                "1" = gen_list(chain),
                "2" = gen_list(chain),
                "3" = gen_list(chain),
                "4" = gen_list(chain),
                "5" = gen_list(chain),
                "6" = gen_list(chain),
                "7" = gen_list(chain),
                "8" = gen_list(chain)
  )
  )
}

# make_c_s_mat: 
# creates C and S matrices for Periodic pulse
# time_steps is a vector of discrete time steps (e.g., 1:8)
# P is the period of the pulse
make_c_s_mat <- function(time_steps = NULL, P = NULL){
  C <- matrix(0, nrow = length(time_steps), ncol = P)
  S <- matrix(0, nrow = length(time_steps), ncol = P)
  for(n in 1:(P)){
    C[,n] <- (cos((2*pi*n*time_steps)/P) * sin((pi*n)/P))*(2/(pi*n))
    S[,n] <- (sin((2*pi*n*time_steps)/P) * sin((pi*n)/P))*(2/(pi*n))
  }
  return(list(C=C, S=S))
}



# make_c_s: 
# creates C and S for periodic boom bust
# time_steps is a vector of discrete time steps (e.g., 1:8)
# P is the period of the pulse
make_c_s <- function( time_steps = NULL, P = NULL){
  
  C <- cos((2 * pi * time_steps)/P)
  S <- sin((2 * pi * time_steps)/P)
  
  return(list(C=C, S=S))
}

# fit_models:
# Fit all 3 models to a species
# and calculate the CPO
# saves posterior as a matrix
# in working directory
# ARGUMENTS 
# models = a character vector of model names
# data_list = a list with data for models
# init_function = a list of functions that generate initial values
#                 for each model in same order as models argument
# to_monitor = a list of character vectors of parameters to track
# species = a character object of the species name, for saving the
#           posterior (added to the file name)
fit_models <- function(models = NULL, data_list = NULL,
                       init_function = NULL, 
                       to_monitor = NULL, species = NULL){
  # for CPO summary stat
  loss_score  <- rep(0, length(models))
  for(i in 1:length(models)){
    mout <- run.jags( model= models[i] , 
                      monitor= to_monitor[[i]] , 
                      data=data_list ,  
                      inits=init_function[[i]] , 
                      n.chains=detectCores()-2 ,
                      adapt=1000,
                      burnin=4000 , 
                      sample=ceiling(200000/(detectCores()-2)) ,
                      thin=2,
                      summarise=FALSE ,
                      plots=FALSE,
                      method = "parallel")
    # convert mcmc.list to matrix
    mmat <- as.matrix(as.mcmc.list(mout), chains = TRUE)
    # get z values
    tz <- as.matrix(mmat[,grep("z\\[", colnames(mmat))])
    # get detection probability
    tp <- as.matrix(mmat[,grep("p\\[", colnames(mmat))])
    tzr <- nrow(tz) # number of sites * surveys
    tzc<- ncol(tz) # number of mcmc samples
    z2 <- as.numeric(tz) * as.numeric(tp) # for CPO
    rm(tp, tz, mout) # clear space
    # number of active days per site and season,
    # ordered like z2
    jnum <- as.numeric(matrix(as.numeric(data_list$jmat), 
                              nrow = tzr, ncol = tzc, byrow = TRUE))
    # number of detections per site and season,
    # ordered like z2
    ynum <- as.numeric(matrix(as.numeric(data_list$y), 
                              nrow = tzr, ncol = tzc, byrow = TRUE))
    # likelihood of observed data given model
    ppd <- dbinom(ynum, jnum, z2)
    # convert to matrix
    mnm <- strsplit(models[i], "\\.")[[1]][1] # get name of model
    ppd <- matrix(ppd, nrow = tzr, ncol = tzc)
    write.table(ppd, paste0(mnm,"_", species,"_ppd.txt"), row.names = FALSE, 
                sep = "\t")
    # don't compute CPO for sites that sampling did not occur.
    to_na <- which(is.na(ppd[1,]))
    # calculate CPO for each site and season
    cpo_vec <- nrow(ppd)/apply(1/ppd, 2, sum, na.rm = TRUE)
    cpo_vec[to_na] <- NA
    cpo_ot <- matrix(cpo_vec, ncol = data_list$nyear, nrow = data_list$nsite)
    cpo_ot <- apply(cpo_ot, 2, function(x) -sum(log(x), na.rm = TRUE))
    write.table(cpo_ot,paste0(mnm,"_", species,"_cpovec.txt") , 
      row.names = FALSE, sep = "\t")
    rm(ppd) # clear space
    loss_score[i] <- -sum(log(cpo_vec), na.rm = TRUE) # summary CPO
    
    # write posterior in working directory
    write.table(mmat, paste0(mnm,"_", species,"_post", ".txt" ), row.names = FALSE,
                sep = "\t")
    rm(mmat) # clear space
  }
  return(data.frame(CPO = loss_score, species = rep(species, length(models)),
                    model = models))
}
  

  
# colfun:
# determines number of colonization events from
# two columns in the z matrix, used in
# zest_posterior function
colfun <- function(x) {
  cold <- x[names(x)=="01"]
  if(length(cold)==0) cold <- 0
  none <- x[names(x)=="00"]
  if(length(none)==0) none <- 0
  ans <- cold/(cold + none)
  if(is.na(ans)) ans <- 0
  return(ans)
}

# zest_posterior
# calculates number of colonization events
# per time step from mcmc matrix of a model
# does this in parallel to speed up computation
zest_posterior <- function(zmat = NULL){
  cl <- makeCluster(detectCores()-2)
  prop_coln <- matrix(0, ncol = dim(zmat)[3]-1, nrow = 3)
  for(j in 2:dim(zmat)[3]){
    z_events <- paste(zmat[,,j-1], zmat[,,j],sep = "") 
    z_events_matrix <- matrix(z_events, nrow = nrow(zmat), ncol = 95)
    z_events_table <- parApply(cl, z_events_matrix, 1, table)
    if(class(z_events_table)=='matrix'){
    z_events_table <- unlist(apply(z_events_table, 2, list), recursive = FALSE)
    }
    prop_coln[,j-1] <- quantile(parSapply(cl, z_events_table, colfun), 
                                probs = c(0.025,0.5,0.975), na.rm = TRUE)
  }
  stopCluster(cl)
  return(prop_coln)
}


# zest_observed:
# similar to zest_posterior, but does it with
# raw detection / non-detection data
  zest_observed <- function(z = NULL) {
    prop_coln <- rep(0, ncol(z)-1)
    for(j in 2:ncol(z)){
      z_events <- paste(z[,j-1], z[,j],sep = "") 
      z_table <- table(z_events)
      prop_coln[j-1] <- colfun(z_table)
    }
    return(prop_coln)
  }
  
  #ex:
  # inverse logit function
  ex <- function(x){
    1 / (1 + exp(-x))
  }
  
  # periodic_pulse_predict:
  # can be used to construct a periodic pulse wave.
  periodic_pulse_predict <- function(A = NULL, delta = NULL, P = NULL, 
                                     time_steps = NULL, m0 = NULL){
    pulse <- matrix(0, ncol = 75, nrow = length(time_steps))
    
    for(n in 1:75){
      pulse[,n] <- (2*A/(n*pi))*sin((pi*n)/P)*
        cos((2 *pi*n * (time_steps - delta))/(P))
    }
    
    ans <- m0 + rowSums(pulse)
    return(ans)
  }
  
  
  #make_pred_pulse:
  # Uses periodic_pulse_predict in parallel
  # to make predictions from the mcmc matrix
  # of a periodic pulse model
  make_pred_PTM_pulse <- function(mcmc_mat, delta = NULL, P = NULL, time_steps = NULL) {
    cl <- makeCluster(detectCores()-2)
    my.env <- new.env()
    my.env$periodic_pulse_predict <- function(A = NULL, delta = NULL, P = NULL, 
                                                time_steps = NULL, m0 = NULL){
      pulse <- matrix(0, ncol = 75, nrow = length(time_steps))
      for(n in 1:75){
        pulse[,n] <- (2*A/(n*pi))*sin((pi*n)/P)*
          cos((2 *pi*n * (time_steps - delta))/(P))
      }
      ans <- m0 + rowSums(pulse)
      return(ans)
    }
    my.env$P <- P
    my.env$time_steps <- time_steps
    my.env$delta <- delta
    A <- mcmc_mat[,grep("^A", colnames(mcmc_mat))]
    m0 <- mcmc_mat[,grep("m0", colnames(mcmc_mat))]
    params <- cbind(A, m0)
    clusterExport(cl, c("periodic_pulse_predict",
      "P","time_steps","delta"), my.env)
    ppul_to_apply <- function(x){periodic_pulse_predict(x[1],
      delta,P,time_steps,x[2])}
    pulse_pred <- t(parApply(cl, params, 1, ppul_to_apply))
    ex <- function(x) 1 / (1 + exp(-x))
    answer <- ex(apply(pulse_pred, 2, quantile, probs = c(0.025,0.5,0.975)))
    stopCluster(cl)
    return(answer)
  }
  
  # periodic_pulse_predict:
  # can be used to construct a sinusoidal wave
  periodic_boom_predict <- function(m0, A, delta, time_steps, P){
    m0 + A*cos(((2*pi)/P)* (time_steps - delta))
  }
  
  #make_pred_boom:
  # Uses periodic_boom_predict in parallel
  # to make predictions from the mcmc matrix
  # of a periodic boom bust model
  
make_pred_PTM_boom <- function(mcmc_mat, delta = NULL, P = NULL, time_steps = NULL) {
  cl <- makeCluster(detectCores()-2)
  my.env <- new.env()
  my.env$periodic_boom_predict <- function(m0, A, delta, time_steps, P){ 
    m0 + A*cos((2*pi*(time_steps - delta))/P)
  }
  my.env$P <- P
  my.env$time_steps <- time_steps
  my.env$delta <- delta
  clusterExport(cl, c("periodic_boom_predict","P","time_steps","delta"), my.env)
  A <- mcmc_mat[,grep("^A", colnames(mcmc_mat))]
  m0 <- mcmc_mat[,grep("m0", colnames(mcmc_mat))]
  params <- cbind(m0, A)
  pboom_to_apply <- function(x){periodic_boom_predict(x[1],x[2],
    delta,time_steps,P)}
  boom_pred <- t(parApply(cl, params, 1, pboom_to_apply))
  ex <- function(x) 1 / (1 + exp(-x))
  answer <- ex(apply(boom_pred, 2, quantile, probs = c(0.025,0.5,0.975)))
  stopCluster(cl)
  rm(my.env)
  return(answer)
}
  
# make_pred_STM
# makes predictions from stochastic time model
# calculates median, 95% CI, and 95% posterior
# predictive distribution
make_pred_STM <- function(mcmc_mat = NULL) {
  # collect parameters
  Ut_sd <- mcmc_mat[,grep("Ut_sd", colnames(mcmc_mat))]
  m0 <- mcmc_mat[,grep("m0", colnames(mcmc_mat))]
  Ut <- mcmc_mat[,grep("Ut\\[", colnames(mcmc_mat))]
  nyear <- ncol(Ut)
  # sample from posterior predictive distribution
  post_pred <- rnorm(1e6, median(m0), Ut_sd)
  # calculate 95% posterior pred on probability scale
  pred_95 <- ex(quantile(post_pred, probs = c(0.025,0.975)))
  # calculate colonization rate at each time step
  Ut_est <- ex(apply(Ut, 2, quantile, probs = c(0.025,0.5,0.975)))
  answer <- matrix(0, ncol = nyear, nrow = 5)
  # combine posterior pred and Ut_est
  answer[1,]<- pred_95[1]
  answer[2:4,] <- Ut_est
  answer[5,]<- pred_95[2]
  return(answer)
}

# make_pred_HTM
# makes predictions from homogeneous time model
# calculates median and 95% CI
make_pred_HTM <- function(mcmc_mat = NULL, nyear = NULL) {
  # collect mu 
  m0 <- mcmc_mat[,grep("m0", colnames(mcmc_mat))]
  # get 
  m0_est <- ex(quantile(m0, probs = c(0.025,0.5,0.975)))
  answer <- matrix(m0_est, ncol = nyear, nrow = 3)
  return(answer)
}
