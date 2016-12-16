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
      g0 = runif(1, -3, 3),
      p0 = runif(1, -3, 3),
      py = runif(12, -3, 3),
      gy = runif(12, -3, 3),
      ly = runif(13, -3, 3),
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

mase <- function(fmat = NULL, dl = NULL, type = "ns", fcast = FALSE){
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
      sum(naive[,9:12], na.rm = TRUE) * (dl$nyear/(dl$nyear-1))
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
      py = runif(12, -3, 3),
      gy = runif(12, -3, 3),
      ly = runif(13, -3, 3),
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

make_c_s_mat <- function( ti = NULL, n = NULL){
  C <- matrix(0, nrow = length(ti), ncol = n)
  S <- matrix(0, nrow = length(ti), ncol = n)
  for(i in 1:n){
    C[,i] <- (cos((pi*i*ti)/2) * sin((pi*i)/4))/(pi*i)
    S[,i] <- (sin((pi*i*ti)/2) * sin((pi*i)/4))/(pi*i)
  }
  return(list(C=C, S=S))
}

fmat <- rany
fmat <- pulsey
a1 <- apply(rany, 2, median)
a2 <- apply(pulsey, 2, median)
plot(a1 - a2)
dl <- data_list
  mase2 <- function(fmat = NULL, dl = NULL, type = "ns"){
    if(type == "ns"){
      # our y values as a vector
      my <- as.numeric(dl$y)
      # our y values as a matrix
      yy <- dl$y
      # naive matrix for denominator
      naive <- matrix(NA, ncol = dl$nyear -1, nrow = dl$nsite)
      for(i in 1:ncol(naive)){
        naive[,i] <- yy[,i+1] - yy[,i] # fill matrix
      }
      # take absolute value
      naive <- abs(naive)
      # calculate naive forecast error
      nf <- sum(naive, na.rm = TRUE) * (dl$nyear/(dl$nyear-1))
      # determine which values to na in our pred
      to_na <- which(is.na(my)==TRUE)
      # na the y_pred values
      fmat[,to_na] <- NA
      # subtract our y values
      ets <- sweep(fmat, 2, my)
      # take absolute value
      ets <- abs(ets)
      # sum error for every sample and divide
      MASE <- rowSums(ets, na.rm = TRUE) / nf
      return(MASE)
    }else{
      fmat <- fmat[,-c(1:100)]
      # our y values as a vector
      my <- as.numeric(dl$y)[-c(1:100)]
      # our y values as a matrix
      yy <- dl$y
      # naive matrix for denominator
      myc <- (dl$nyear -1 - dl$P)
      seas <- matrix(NA, ncol = myc, nrow = dl$nsite)
      for(i in 1:ncol(seas)){
        seas[,i] <- yy[,i+dl$P+1] - yy[,i + 1] # fill matrix
      }
      # take absolute value
      seas <- abs(seas)
      # calculate naive forecast error
      nf <- sum(seas, na.rm = TRUE) * (dl$nyear/(dl$nyear-dl$P))
      # determine which values to na in our pred
      to_na <- which(is.na(my)==TRUE)
      # na the y_pred values
      fmat[,to_na] <- NA
      # subtract our y values
      ets <- sweep(fmat, 2, my)
      # take absolute value
      ets <- abs(ets)
      # sum error for every sample and divide
      MASE <- rowSums(ets, na.rm = TRUE) / nf
      return(MASE)
    }
  }
  
  pploss <- function(post = NULL, yobs = NULL){
    yobs <- as.numeric(yobs)
    a1 <- abs(sweep(post, 2, yobs))
    togo <- which(is.na(yobs))
    am <- apply(a1[,-togo], 2, mean)
    
    a2 <- apply(post[,-togo], 2, var)
    
    dsel <- sum(am) + sum(a2)
    return(dsel)
  }