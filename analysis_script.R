
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

#####################################################################
# Assuming the files are saved in your working directory and you have
### loaded the above functions
#####################################################################

# load the reshape package
package_load("reshape2")


# lets do the coyote
z <- df_2_array(read.table("z_matrix_sp10_sp13.txt", header = TRUE, sep = "\t"))[1,,]

# get number of sites occupied


soc <- colSums(z, na.rm = TRUE) / apply(z, 2, function(x) 100 - sum(is.na(x)))

# get beta a and b from mean and standard deviation of soc

betaABfromMeanSD = function( mean , sd ) {
  if ( mean <=0 | mean >= 1) stop("must have 0 < mean < 1")
  if ( sd <= 0 ) stop("sd must be > 0")
  kappa = mean*(1-mean)/sd^2 - 1
  if ( kappa <= 0 ) stop("invalid combination of mean and sd")
  a = mean * kappa
  b = ( 1.0 - mean ) * kappa
  return( list( a=a , b=b ) )
}


prior_for_occ <- betaABfromMeanSD(mean(soc), sd(soc))
# figure out how many sites there were sampled each season



# read in the species names. Note: The 1st dimension of the y and z array
### are ordered this way. This means that species specific covariate data
### should also be ordered identically.
species_names <- read.table("species_used_in_fa10_sp13_analysis.txt", header = TRUE)

# read in the site names. Note: the 2nd dimension of the y and z array
### are ordered this way. This means that covariate data should also
### be ordered identically.  Furthermore, the 1st dimension of the
### j matrix is ordered this way
site_names <- read.table("sites_used_in_sp10_sp13_analysis.txt", header = TRUE)

covdat <- read.csv("urban_pc.csv", header = TRUE)
# read in the y array

y_array <- df_2_array(read.table("y_matrix_sp10_sp13.txt", header = TRUE, sep = "\t"))

j_mat <- read.table("j_matrix_sp10_sp13.txt", header = TRUE, sep = "\t")



inits <- function(chain){
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

params <- c("lp", "ly", "g0", "p0", "gy", "py", "gy_sd",
            "py_sd", "ly_sd", "psinit")

# get coyote data from the y_array
data_list <- list(y = as.matrix(y_array[which(species_names$x=="Coyote"),,]), nyear = ncol(z), 
                  nsite = nrow(z), 
                  spa = prior_for_occ$a, spb = prior_for_occ$b,
                  jmat = as.matrix(j_mat),
                  cov = covdat$pc1)#,
                  #pi = 3.14159)

data_list$jmat[68,4] <- floor(median(data_list$jmat))

# run the jags model.


mod_mcmc <- run.jags( model= "ranef_year_jags.R" , 
                                   monitor=params , 
                                   data=data_list ,  
                                   inits=inits , 
                                   n.chains=7 ,
                                   adapt=3000,
                                   burnin=3000 , 
                                   sample=ceiling(10000/7) ,
                                   thin=5 ,
                                   summarise=FALSE ,
                                   plots=FALSE,
                                   method = "parallel")
k <- 100
tok <- which(!is.na(data_list$y))
tok <- tok[-which(tok<101)]

####
sim_ans <- rep(0,10)
for(i in 2:10){
  y <- data_list$y
tg <- sample(1:length(tok), k, replace = TRUE)
stored_val <- rep(0, length(tg))
stored_val <- y[tok[tg]]
y[tok[tg]] <- NA

new_list <- list( y = as.matrix(y), nyear = ncol(z), nsite = nrow(z),
                  spa = prior_for_occ$a, spb = prior_for_occ$b,
                  jmat = as.matrix(j_mat),
                  cov = covdat$pc1)
new_list$jmat[68,4] <- floor(median(data_list$jmat))
  dropped <- run.jags( model= "ranef_year_jags.R" , 
                        monitor="y" , 
                        data=new_list ,  
                        inits=inits , 
                        n.chains=7 ,
                        adapt=3000,
                        burnin=3000 , 
                        sample=ceiling(10000/7) ,
                        thin=5 ,
                        summarise=FALSE ,
                        plots=FALSE,
                        method = "parallel")
  dropped <- as.matrix(as.mcmc.list(dropped), chains = TRUE)[,-1]
  # get just the ones we predicted
  preds <- dropped[,tok[tg]]
  bfs <- (sweep(preds, 2, stored_val)^2)
  sim_ans[i] <- mean(bfs)
}

fy <- data_list$y
fy[,10:13] <- NA
forecast_list <- list( y = as.matrix(y), nyear = ncol(z), nsite = nrow(z),
                  spa = prior_for_occ$a, spb = prior_for_occ$b,
                  jmat = as.matrix(j_mat),
                  cov = covdat$pc1)

dropped <- run.jags( model= "ranef_year_jags.R", 
                     monitor="y_pred" , 
                     data=new_list ,  
                     inits=inits , 
                     n.chains=7 ,
                     adapt=3000,
                     burnin=3000 , 
                     sample=ceiling(10000/7) ,
                     thin=5 ,
                     summarise=FALSE ,
                     plots=FALSE,
                     method = "parallel")
dm <- as.matrix(as.mcmc.list(dropped), chains = TRUE)[,-1]

m1 <- mase(dm, data_list, type = "season")

# fix this function so we can use it for both types
# of effects
mase <- function(fmat = NULL, dl = NULL, type = "ns"){
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
  ets <- sweep(dm, 2, my)
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
      
  }
}
j1 <- mm[1,-c(1:100)]
j2 <- dm[2,-c(1:100)]
j1[to_na] <- NA
j2[to_na] <- NA
jj1 <- abs(j1 - my)
jj2 <- abs(j2 - my)
#
mm <- as.jags(mod_mcmc)

source("drop.k_travis.R")
test <- drop.k(mod_mcmc, dropvars = 'y[1:100, 1:13]', k = 100, simulations = 2, n.cores = 6,
               max.time = "5m", silent.jags = FALSE)


test <- drop.k()

# do it again with the trig functions

inits <- function(chain){
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

cs <- make_c_s_mat(1:12, 3)
y
data_list <- list(y = as.matrix(y_array[which(species_names$x=="Coyote"),,]), nyear = ncol(z), 
                  nsite = nrow(z), 
                  spa = prior_for_occ$a, spb = prior_for_occ$b,
                  jmat = as.matrix(j_mat),
                  pi = 3.14159, C = cs[[1]], S = cs[[2]],
                  cov = covdat$pc1)
data_list$jmat[68,4] <- 26

params <- c("lp", "ly", "p0", "gy", "py", "gy_sd",
            "py_sd", "ly_sd", "psinit", "theta",
            "dprobs", "dp", "g_mu", "y_pred")
j_mat[68,4] <- 5
mod_pulse2 <- run.jags( model= "pulse_year_jags_trig.R" , 
                                   monitor= params , 
                                   data=data_list ,  
                                   inits=inits , 
                                   n.chains=7 ,
                                   adapt=3000,
                                   burnin=3000 , 
                                   sample=ceiling(10000/7) ,
                                   thin=5 ,
                                   summarise=FALSE ,
                                   plots=FALSE,
                                   method = "parallel")
mm <- as.matrix(as.mcmc.list(mod_pulse2), chains = TRUE)[,-1]
mm <- mm[,grep("y_pred", colnames(mm))]

tt <- mase(mm, data_list, type = "ns")


theta <- mm[,45]
delta <- mm[,47]
gys <- mm[,grep("gy", colnames(mm))]

gys <- apply(gys, 2, median)

plot(colSums(z, na.rm = TRUE) /(100 - apply(z, 2, function(x) sum(is.na(x)))),
     type = "l")

zs <- mm[,grep("z", colnames(mm))]

ans <- matrix(0, ncol = 12, nrow = 100)
av <- matrix(0, nrow = 4000, ncol = 12)
for(j in 1:4000){
  nz <- matrix(zs[j,], ncol = 13, nrow = 100)
for(i in 2:13){
  ans[,i-1] <- nz[,i] - nz[,i-1]
  boop <- nz[,i] + nz[,i-1]
  av[j,i-1] <- sum(ans[,i-1]==1, na.rm = TRUE)/(100 - sum(boop==2))
}
}

mye <- apply(av, 2, quantile, probs= c(0.05,0.5,0.95))

for(i in 2:12){
  ans[,i-1] <- z[,i] - z[,i-1]
  av[i-1] <- sum(ans[,i-1]==1, na.rm = TRUE)/(100 - sum(is.na(ans[,i-1])))
}

############### opossum ###############

inits <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = z,
      g0 = runif(1, -3, 3),
      p0 = runif(1, -3, 3),
      py = runif(12, -3, 3),
      gy = runif(12, -3, 3),
      ly = runif(13, -3, 3),
      lp = runif(1, -3, 3),
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

params <- c("lp", "ly", "g0", "p0", "gy", "py", "gy_sd",
            "py_sd", "ly_sd", "psinit")

# get coyote data from the y_array
data_list <- list(y = as.matrix(y_array[which(species_names$x=="Opossum"),,]), nyear = ncol(z), 
                  nsite = nrow(z), 
                  spa = prior_for_occ$a, spb = prior_for_occ$b,
                  jmat = as.matrix(j_mat),
                  pi = 3.14159)

# run the jags model.

mod_mcmc <- as.mcmc.list(run.jags( model= "ranef_year_jags.R" , 
                                   monitor=params , 
                                   data=data_list ,  
                                   inits=inits , 
                                   n.chains=1 ,
                                   adapt=3000,
                                   burnin=3000 , 
                                   sample=10000 ,
                                   thin=5 ,
                                   summarise=FALSE ,
                                   plots=FALSE,
                                   method = "parallel"))

make_c_s <- function( ti = NULL, p = NULL){

  C <- cos((2 * pi * ti)/p)
  S <- cos((2 * pi * ti)/p)

  return(list(C=C, S=S))
}

cs <- make_c_s(1:12, p = 2)

data_list <- list(y = as.matrix(y_array[which(species_names$x=="Opossum"),,]), nyear = ncol(z), 
                  nsite = nrow(z), 
                  spa = prior_for_occ$a, spb = prior_for_occ$b,
                  jmat = as.matrix(j_mat),
                  pi = 3.14159, C = cs[[1]], S = cs[[2]], P = 2)

inits <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = z,
      g0 = runif(1, -3, 3),
      p0 = runif(1, -3, 3),
      py = runif(12, -3, 3),
      gy = runif(12, -3, 3),
      ly = runif(13, -3, 3),
      lp = runif(1, -3, 3),
      a1_phi = runif(1, 0.01, 3),
      a1_gam = runif(1, 0.01, 3),
      dprobs_phi = 1,
      dprobs_gam = 0,
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


params <- c("lp", "ly", "p0", "gy", "py", "gy_sd",
            "py_sd", "ly_sd", "psinit", "a1_phi", "a1_gam",
            "dprobs_gam", "dprobs_phi",
            "dp", "g_mu", "z")

mod_boom <- as.mcmc.list(run.jags( model= "boom_bust_jags_trig.R" , 
                                     monitor=params , 
                                     data=data_list ,  
                                     inits=inits , 
                                     n.chains=1 ,
                                     adapt=3000,
                                     burnin=3000 , 
                                     sample=4000 ,
                                     thin=5 ,
                                     summarise=FALSE ,
                                     plots=FALSE,
                                     method = "parallel"))


om <- as.matrix(mod_boom, chains = TRUE)



a1p <- om[,45]
a1g <- om[,46]

a2p <- om[,47]
a2g <- om[,48]
gys <- om[,grep("gy\\[", colnames(om))]
pys<- om[,grep("py\\[", colnames(om))]

gys <- apply(gys, 2, median)

zs <- om[,grep("z", colnames(om))]

ans <- matrix(0, ncol = 12, nrow = 100)
av <- matrix(0, nrow = 4000, ncol = 12)
for(j in 1:4000){
  nz <- matrix(zs[j,], ncol = 13, nrow = 100)
  for(i in 2:13){
    ans[,i-1] <- nz[,i] - nz[,i-1]
    boop <- nz[,i] + nz[,i-1]
    av[j,i-1] <- sum(ans[,i-1]==1, na.rm = TRUE)/(100 - sum(boop==2))
  }
}

mye <- apply(av, 2, quantile, probs= c(0.05,0.5,0.95))
