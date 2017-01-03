# Utility sript to bring in the y, z array and j matrix to be used for 
### occupancy analysis.
package_load<-function(packages = c("reshape2"), quiet=TRUE, verbose=FALSE, warn.conflicts=FALSE){
  
  # download required packages if they're not already
  pkgsToDownload<- packages[!(packages  %in% installed.packages()[,"Package"])]
  if(length(pkgsToDownload)>0)
    install.packages(pkgsToDownload, repos="http://cran.us.r-project.org", quiet=quiet, verbose=verbose)
  
  # then load them
  for(i in 1:length(packages))
    require(packages[i], character.only=T, quietly=quiet, warn.conflicts=warn.conflicts)
}

# change the dataframe back to an array
df_2_array <- function(my_df = NULL){
  
  require(reshape2)
  my_array <- acast(my_df, species~site~season, value.var = "count")
  dimnames(my_array) <- NULL
  return(my_array)
}


betaABfromMeanSD = function( mean , sd ) {
  if ( mean <=0 | mean >= 1) stop("must have 0 < mean < 1")
  if ( sd <= 0 ) stop("sd must be > 0")
  kappa = mean*(1-mean)/sd^2 - 1
  if ( kappa <= 0 ) stop("invalid combination of mean and sd")
  a = mean * kappa
  b = ( 1.0 - mean ) * kappa
  return( list( a=a , b=b ) )
}


inits_ranef <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = z,
      g_mu = runif(1, -3, 3),
      p_mu = runif(1, -3, 3),
      py = runif(8, -3, 3),
      gy = runif(8, -3, 3),
      ly = runif(9, -3, 3),
      lp = runif(1, -3, 3),
      ls = runif(95, -3, 3),
      p1 = runif(1, -3, 3),
      g1 = runif(1, -3, 3),
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
inits_homog <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = z,
      g_mu = runif(1, -3, 3),
      py = runif(8, -3, 3),
      ly = runif(9, -3, 3),
      ls = runif(95,-3,3),
      lp = runif(1, -3, 3),
      p1 = runif(1, -3, 3),
      g1 = runif(1, -3, 3),
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


dl <- 

mase <- function(fmat = NULL, dl = NULL, type = "naive", fcast = FALSE){
  # our y values as a vector
  my <- as.numeric(dl$y)
  # our y values as a matrix
  yy <- dl$y
  if(!any(fmat>1)){
    my[which(my>1)] <- 1
    yy[which(yy>1)] <- 1
  }
  #if(fcast){
   # dl$nyear <- 5
    #if(type == "naive") yy <- yy[,9:13]
    #if(type == "season")yy <- yy[,6:13]
    #my <- as.numeric(yy)
  #}
  
  if(type == "naive"){
    # naive matrix for denominator
    naive <- matrix(NA, ncol = dl$nyear -1, nrow = dl$nsite)
    for(i in 1:ncol(naive)){
      naive[,i] <- yy[,i+1] - yy[,i] # fill matrix
    }
    # take absolute value
    naive <- abs(naive)
    # calculate naive forecast error
    if(fcast){
      nf <- sum(naive[,9:12], na.rm = TRUE) * (dl$nyear/(dl$nyear-1))
    }else{
      nf <- sum(naive, na.rm = TRUE) * (dl$nyear/(dl$nyear-1))
    }
    
    # determine which values to na in our pred
    to_na <- which(is.na(my)==TRUE)
    # na the y_pred values
    fmat[,to_na] <- NA
    # subtract our y values
    ets <- sweep(fmat, 2, my)
    # take absolute value
    ets <- abs(ets)
    # sum error for every sample and divide
    if(fcast){
      MASE <- rowSums(ets[,901:1300], na.rm = TRUE) / nf
    }else{
    MASE <- rowSums(ets, na.rm = TRUE) / nf
    }
    return(MASE)
  }
  if(type == "season"){

    # seasonal matrix for denominator
    myc <- (dl$nyear - dl$P)
    #if(fcast) myc <- 4
    seas <- matrix(NA, ncol = myc, nrow = dl$nsite)
    for(i in 1:ncol(seas)){
      seas[,i] <- yy[,i+dl$P] - yy[,i] # fill matrix
    }
    # take absolute value
    seas <- abs(seas)
    # calculate naive forecast error
    if(fcast){
      nf <- sum(seas[,6:9], na.rm = TRUE) * (dl$nyear/(dl$nyear-dl$P))
    }else{
      nf <- sum(seas, na.rm = TRUE) * (dl$nyear/(dl$nyear-dl$P))
    }
    
    # determine which values to na in our pred
    to_na <- which(is.na(my)==TRUE)
    # na the y_pred values
    fmat[,to_na] <- NA
    # subtract our y values
    ets <- sweep(fmat, 2, my)
    # take absolute value
    ets <- abs(ets)
    # sum error for every sample and divide
    if(fcast){
      MASE <- rowSums(ets[,901:1300], na.rm = TRUE) / nf
    }else{
      MASE <- rowSums(ets, na.rm = TRUE) / nf
    }
    return(MASE)
    
  }
}


inits_pulse <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = z,
      g1 = runif(1, -3, 3),
      p1 = runif(1, -3, 3),
      py = runif(8, -3, 3),
      gy = runif(8, -3, 3),
      ly = runif(9, -3, 3),
      ls = runif(100,-3,3),
      lp = runif(1, -3, 3),
      g_mu = runif(1, -3, 3),
      p_mu = runif(1, -3, 3),
      theta = runif(1, 0.01, 3),
      dprobs = 2,
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

inits_boom <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = z,
      g1 = runif(1, -3, 3),
      p1 = runif(1, -3, 3),
      py = runif(8, -3, 3),
      gy = runif(8, -3, 3),
      ly = runif(9, -3, 3),
      ls = runif(95,-3,3),
      lp = runif(1, -3, 3),
      g_mu = runif(1, -3, 3),
      p_mu = runif(1, -3, 3),
      a1_gam = runif(1, 0.01, 3),
      dprobs_gam = 2,
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

inits_only_pulse <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = z,
      g1 = runif(1, -3, 3),
      p1 = runif(1, -3, 3),
      py = runif(8, -3, 3),
      ly = runif(9, -3, 3),
      lp = runif(1, -3, 3),
      ls = runif(100,-3,3),
      g_mu = runif(1, -3, 3),
      p_mu = runif(1, -3, 3),
      theta = runif(1, 0.01, 3),
      dprobs = 2,
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
inits_only_boom <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = z,
      g1 = runif(1, -3, 3),
      p1 = runif(1, -3, 3),
      py = runif(8, -3, 3),
      ly = runif(9, -3, 3),
      lp = runif(1, -3, 3),
      ls = runif(95,-3,3),
      g_mu = runif(1, -3, 3),
      p_mu = runif(1, -3, 3),
      a1_gam = runif(1, 0.01, 3),
      dprobs_gam = 2,
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

make_c_s_mat <- function( ti = NULL, n = NULL){
  C <- matrix(0, nrow = length(ti), ncol = n)
  S <- matrix(0, nrow = length(ti), ncol = n)
  for(i in 1:n){
    C[,i] <- (cos((pi*i*ti)/2) * sin((pi*i)/4))/(pi*i)
    S[,i] <- (sin((pi*i*ti)/2) * sin((pi*i)/4))/(pi*i)
  }
  return(list(C=C, S=S))
}



  pploss <- function(post = NULL, yobs = NULL){
    yobs <- as.numeric(yobs)
    a1 <- sweep(post, 2, yobs)^2
    togo <- which(is.na(yobs))
    am <- apply(a1[,-togo], 2, mean)
    
    a2 <- apply(post[,-togo], 2, var)
    
    dsel <- sum(am) + sum(a2)
    return(dsel)
  }
  
  
  make_c_s <- function( ti = NULL, p = NULL){
    
    C <- cos((2 * pi * ti)/p)
    S <- sin((2 * pi * ti)/p)
    
    return(list(C=C, S=S))
  }

  
  fit_models <- function(models = NULL, dl = NULL,
                         inl = NULL, to_monitor = NULL, species = NULL){
    
    loss_score <- loss_score2 <- rep(0, length(models))
    for(i in 1:length(models)){
      
      mout <- run.jags( model= models[i] , 
                        monitor= to_monitor[[i]] , 
                        data=dl ,  
                        inits=inl[[i]] , 
                        n.chains=detectCores()-1 ,
                        adapt=3000,
                        burnin=3000 , 
                        sample=ceiling(10000/(detectCores()-1)) ,
                        thin=5,
                        summarise=FALSE ,
                        plots=FALSE,
                        method = "parallel")
      mmat <- as.matrix(as.mcmc.list(mout), chains = TRUE)
      tz <- as.matrix(mmat[,grep("z\\[", colnames(mmat))])
      tp <- as.matrix(mmat[,grep("p\\[", colnames(mmat))])
      jnum <- as.numeric(matrix(as.numeric(dl$jmat), nrow = nrow(tz), ncol = ncol(tz), byrow = TRUE))
      ynum <- as.numeric(matrix(as.numeric(dl$y), nrow = nrow(tz), ncol = ncol(tz), byrow = TRUE))
      z2 <- as.numeric(tz) * as.numeric(tp)
      ppd <- dbinom(ynum, jnum, z2)
      ppd <- matrix(ppd, nrow = nrow(tz), ncol = ncol(tz))
      cpo_vec <- nrow(ppd)/apply(1/ppd, 2, sum, na.rm = TRUE)
      tg <- which(is.na(as.numeric(dl$y)))
      loss_score[i] <- -sum(log(cpo_vec[-tg]))
      loss_score2[i] <- pploss(mmat[,grep("y_pred", colnames(mmat))], dl$y)
      mnm <- strsplit(models[i], "\\.")[[1]][1]
      write.table(mmat, paste0("./model_outputs/",species,"_", mnm,".txt" ), row.names = FALSE,
                  sep = "\t")
      rm(mmat)
    }
    return(data.frame(CPO = loss_score, pploss = loss_score2, species = rep(species, length(models)),
                      model = models))
  }
  
  
  plsg <- function(theta = NULL, delta = NULL, p = NULL, t = NULL, a0 = NULL){
    
    n <- 1:3  
    a1 <- (theta/(1*pi))*sin((pi*1)/p)*cos((2 *pi*1 * (t - delta))/(p))
    a2 <- (theta/(2*pi))*sin((pi*2)/p)*cos((2 *pi*2 * (t - delta))/(p))
    a3 <- (theta/(3*pi))*sin((pi*3)/p)*cos((2 *pi*3 * (t - delta))/(p))
    a4 <- (theta/(4*pi))*sin((pi*4)/p)*cos((2 *pi*4 * (t - delta))/(p))
    a5 <- (theta/(5*pi))*sin((pi*5)/p)*cos((2 *pi*5 * (t - delta))/(p))
    a6 <- (theta/(6*pi))*sin((pi*6)/p)*cos((2 *pi*6 * (t - delta))/(p))
    a7 <- (theta/(7*pi))*sin((pi*7)/p)*cos((2 *pi*7 * (t - delta))/(p))
    a8 <- (theta/(8*pi))*sin((pi*8)/p)*cos((2 *pi*8 * (t - delta))/(p))
    a9 <- (theta/(9*pi))*sin((pi*9)/p)*cos((2 *pi*9 * (t - delta))/(p))
    a10 <- (theta/(10*pi))*sin((pi*10)/p)*cos((2 *pi*10 * (t - delta))/(p))
    a11 <- (theta/(11*pi))*sin((pi*11)/p)*cos((2 *pi*11 * (t - delta))/(p))
    
    ans <- a0 + a1 +a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 
    return(ans)
  }
  
  p2 <- function(a0, a1, a2, ti, p){
    a0 + a1*cos(((2*pi)/p)* (ti - a2))
  }