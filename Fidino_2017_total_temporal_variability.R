############################
#
# Total temporal variability script for Fidino and Magle 2017:
# Using Fourier series to estimate periodic
# trends in dynamic occupancy models.
#
#


### DESCRIPTION ###### DESCRIPTION ###### DESCRIPTION ###### DESCRIPTION ###
# This script uses the functions in Fidino_2017_periodic_utility_functions.R
# to fit a periodic time model to each species that also includes
# a random temporal component. The primary variable of interest
# from these models is the standard devation associated to the
# random temporal component, as it can be compared to a 
# the standard deviation from a model that only contains
# a random temporal component (i.e., no periodic element).
### DESCRIPTION ###### DESCRIPTION ###### DESCRIPTION ###### DESCRIPTION ###


# read in the utility functions, assuming it is in the
# working directory
source("Fidino_2017_periodic_utility_functions.R")

#####################################################################
# Assuming the files are saved in your working directory and you have
### sourced the utility script
#####################################################################

# load the necessary packages for this analysis
package_load(c("reshape2", "runjags", "rjags", "parallel", "data.table"))
# load glm module for JAGS
load.module('glm')


# the community incidence matrix
# z is a species * site * season array.
# Species order in z is
# 1. Coyote
# 2. Red Fox
# 3. Striped Skunk
# 4. Raccoon
# 5. Virginia Opossum

#
# COYOTE ANALYSIS, reading in just the first species
#
z <- df_2_array(read.table("Fidino_2017_community_incidence_matrix_sp11_sp13.txt", header = TRUE, sep = "\t"))[1,,]

# proportion of sites occupied each time step
prop_sites_occupied <- colSums(z, na.rm = TRUE) / apply(z, 2, function(x) 95 - sum(is.na(x)))

# get beta a and b from mean and standard deviation of soc
prior_for_occ <- betaABfromMeanSD(mean(prop_sites_occupied), sd(prop_sites_occupied))

covdat <- read.table("Fidino_2017_URB_covariate.txt", header = TRUE, sep ="\t")

# read in the y array (days a species was seen per site and season)
# ordered identically to community incidence matrix z
y_array <- df_2_array(read.table("Fidino_2017_chicago_detection_data_sp11_sp13.txt", header = TRUE, sep = "\t"))

# read in number of days a camera was active per site and season
j_mat <- read.table("Fidino_2017_Chicago_days_camera_active_sp11_sp13.txt", header = TRUE, sep = "\t")

# make temporally varying covariates for Fourier analysis
cs <- make_c_s_mat(1:8, 4)


# data to supply to model
data_list <- list(y = as.matrix(y_array[1,,]), nyear = ncol(z), 
  nsite = nrow(z), 
  spa = prior_for_occ$a, spb = prior_for_occ$b,
  jmat = as.matrix(j_mat),
  pi = 3.14159, C = cs[[1]], S = cs[[2]],
  x = covdat$x, P = 4,
  delta = 2)




# paramaters to track from each model
to_monitor <-  c("psinit", "Gt","Gt_sd", "d0","d1","m0",
  "m1", "A", "delta", "omega", "omega_sd",
  "epsilon", "epsilon_sd","f0", "f1", "Ut", "Ut_sd")

inits <- function(chain){
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
      A = rgamma(1,1,1),
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

coyote_out <- run.jags( model= "Fidino_2017_periodic_time_full_model_pulse.R" , 
  monitor= to_monitor, 
  data=data_list ,  
  inits=inits , 
  n.chains=detectCores()-2 ,
  adapt=1000,
  burnin=4000 , 
  sample=ceiling(200000/(detectCores()-2)) ,
  thin=2,
  summarise=FALSE ,
  plots=FALSE,
  method = "parallel")


coyote_out <- as.matrix(as.mcmc.list(coyote_out), chains = TRUE)
# save CPO scores for each model 
write.table(coyote_out, "coyote_full_model.txt", row.names = FALSE, sep = "\t")
rm(coyote_out)
#
# RED FOX ANALYSIS, reading in the second species
#
z <- df_2_array(read.table(
  "Fidino_2017_community_incidence_matrix_sp11_sp13.txt", 
  header = TRUE, sep = "\t"))[2,,]

prop_sites_occupied <- colSums(z, na.rm = TRUE) / 
  apply(z, 2, function(x) 95 - sum(is.na(x)))

prior_for_occ <- betaABfromMeanSD(mean(prop_sites_occupied), 
  sd(prop_sites_occupied))

data_list <- list(y = as.matrix(y_array[2,,]), nyear = ncol(z), 
  nsite = nrow(z), 
  spa = prior_for_occ$a, spb = prior_for_occ$b,
  jmat = as.matrix(j_mat),
  pi = 3.14159, C = cs[[1]], S = cs[[2]],
  x = covdat$x, P = 4,
  delta = 2)

fox_out <- run.jags( model= "Fidino_2017_periodic_time_full_model_pulse.R" , 
  monitor= to_monitor, 
  data=data_list ,  
  inits=inits , 
  n.chains=detectCores()-2 ,
  adapt=1000,
  burnin=4000 , 
  sample=ceiling(200000/(detectCores()-2)) ,
  thin=2,
  summarise=FALSE ,
  plots=FALSE,
  method = "parallel")

fox_out <- as.matrix(as.mcmc.list(fox_out), chains = TRUE)
# save CPO scores for each model 
write.table(fox_out, "fox_full_model.txt", row.names = FALSE, sep = "\t")
rm(fox_out)

#
# STRIPED SKUNK ANALYSIS, reading in the third species
#

z <- df_2_array(read.table(
  "Fidino_2017_community_incidence_matrix_sp11_sp13.txt", 
  header = TRUE, sep = "\t"))[3,,]

prop_sites_occupied <- colSums(z, na.rm = TRUE) / 
  apply(z, 2, function(x) 95 - sum(is.na(x)))

prior_for_occ <- betaABfromMeanSD(mean(prop_sites_occupied), 
  sd(prop_sites_occupied))

data_list <- list(y = as.matrix(y_array[3,,]), nyear = ncol(z), 
  nsite = nrow(z), 
  spa = prior_for_occ$a, spb = prior_for_occ$b,
  jmat = as.matrix(j_mat),
  pi = 3.14159, C = cs[[1]], S = cs[[2]],
  x = covdat$x, P = 4,
  delta = 2)

skunk_out <- run.jags( model= "Fidino_2017_periodic_time_full_model_pulse.R", 
  monitor= to_monitor, 
  data=data_list ,  
  inits=inits , 
  n.chains=detectCores()-2 ,
  adapt=1000,
  burnin=4000 , 
  sample=ceiling(200000/(detectCores()-2)),
  thin=2,
  summarise=FALSE ,
  plots=FALSE,
  method = "parallel")

skunk_out <- as.matrix(as.mcmc.list(skunk_out), chains = TRUE)
# save CPO scores for each model 
write.table(skunk_out, "skunk_full_model.txt", row.names = FALSE, sep = "\t")
rm(skunk_out)

#
# RACCOON ANALYSIS, reading in the fourth species
#

z <- df_2_array(read.table(
  "Fidino_2017_community_incidence_matrix_sp11_sp13.txt", 
  header = TRUE, sep = "\t"))[4,,]

prop_sites_occupied <- colSums(z, na.rm = TRUE) / 
  apply(z, 2, function(x) 95 - sum(is.na(x)))

prior_for_occ <- betaABfromMeanSD(mean(prop_sites_occupied), 
  sd(prop_sites_occupied))

# make C and S for the boom bust Fourier series
cs <- make_c_s(1:8, 2)


data_list <- list(y = as.matrix(y_array[4,,]), nyear = ncol(z), 
  nsite = nrow(z), 
  spa = prior_for_occ$a, spb = prior_for_occ$b,
  jmat = as.matrix(j_mat),
  pi = 3.14159, C = cs[[1]], S = cs[[2]],
  x = covdat$x, P = 2,
  delta = 0)

raccoon_out <- run.jags(
  model= "Fidino_2017_periodic_time_full_model_boom_bust.R", 
  monitor= to_monitor, 
  data=data_list ,  
  inits=inits , 
  n.chains=detectCores()-2 ,
  adapt=1000,
  burnin=4000 , 
  sample=ceiling(200000/(detectCores()-2)),
  thin=2,
  summarise=FALSE ,
  plots=FALSE,
  method = "parallel")

raccoon_out <- as.matrix(as.mcmc.list(raccoon_out), chains = TRUE)
# save CPO scores for each model 
write.table(raccoon_out, "raccoon_full_model.txt", 
  row.names = FALSE, sep = "\t")
rm(raccoon_out)

#
# VIRGINIA OPOSSUM ANALYSIS, reading in the fifth species
#

z <- df_2_array(read.table(
  "Fidino_2017_community_incidence_matrix_sp11_sp13.txt", 
  header = TRUE, sep = "\t"))[5,,]

prop_sites_occupied <- colSums(z, na.rm = TRUE) / 
  apply(z, 2, function(x) 95 - sum(is.na(x)))

prior_for_occ <- betaABfromMeanSD(mean(prop_sites_occupied), 
  sd(prop_sites_occupied))

data_list <- list(y = as.matrix(y_array[5,,]), nyear = ncol(z), 
  nsite = nrow(z), 
  spa = prior_for_occ$a, spb = prior_for_occ$b,
  jmat = as.matrix(j_mat),
  pi = 3.14159, C = cs[[1]], S = cs[[2]],
  x = covdat$x, P = 2, delta = 0)

opossum_out <- run.jags(
  model= "Fidino_2017_periodic_time_full_model_boom_bust.R" , 
  monitor= to_monitor, 
  data=data_list ,  
  inits=inits , 
  n.chains=detectCores()-2 ,
  adapt=1000,
  burnin=4000 , 
  sample=ceiling(200000/(detectCores()-2)) ,
  thin=2,
  summarise=FALSE ,
  plots=FALSE,
  method = "parallel")

opossum_out <- as.matrix(as.mcmc.list(opossum_out), chains = TRUE)
# save CPO scores for each model 
write.table(opossum_out, "opossum_full_model.txt",
  row.names = FALSE, sep = "\t")
rm(opossum_out)

### calc_proportion variance explained by periodic trend

# get standardard devation from 
m1_sd <- fread("coyote_full_model.txt", select = "Ut_sd", data.table = FALSE)
m2_sd <- fread("Fidino_2017_stochastic_time_model_coyote_post.txt", 
  select = "Ut_sd", data.table = FALSE)

m_sd <- apply(cbind(m1_sd, m2_sd), 2, median)

variance_explained <- rep(0, 5)
variance_explained[1] <- 1 - (m_sd[1]/m_sd[2])

m1_sd <- fread("fox_full_model.txt", select = "Ut_sd", data.table = FALSE)
m2_sd <- fread("Fidino_2017_stochastic_time_model_redfox_post.txt", 
  select = "Ut_sd", data.table = FALSE)

m_sd <- apply(cbind(m1_sd, m2_sd), 2, median)
variance_explained[2] <- 1 - (m_sd[1]/m_sd[2])

m1_sd <- fread("skunk_full_model.txt", select = "Ut_sd", data.table = FALSE)
m2_sd <- fread("Fidino_2017_stochastic_time_model_skunk_post.txt", 
  select = "Ut_sd", data.table = FALSE)

m_sd <- apply(cbind(m1_sd, m2_sd), 2, median)
variance_explained[3] <- 1 - (m_sd[1]/m_sd[2])

m1_sd <- fread("raccoon_full_model.txt", select = "Ut_sd", data.table = FALSE)
m2_sd <- fread("Fidino_2017_stochastic_time_model_raccoon_post.txt", 
  select = "Ut_sd", data.table = FALSE)

m_sd <- apply(cbind(m1_sd, m2_sd), 2, median)
variance_explained[4] <- 1 - (m_sd[1]/m_sd[2])

m1_sd <- fread("opossum_full_model.txt", select = "Ut_sd", data.table = FALSE)
m2_sd <- fread("Fidino_2017_stochastic_time_model_opossum_post.txt", 
  select = "Ut_sd", data.table = FALSE)

m_sd <- apply(cbind(m1_sd, m2_sd), 2, median)
variance_explained[5] <- 1 - (m_sd[1]/m_sd[2])


