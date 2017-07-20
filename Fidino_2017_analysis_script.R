############################
#
# Analysis script for Fidino and Magle 2017:
# Using Fourier series to estimate periodic
# trends in dynamic occupancy models.
#
#


### DESCRIPTION ###### DESCRIPTION ###### DESCRIPTION ###### DESCRIPTION ###
# This script uses the functions in Fidino_2017_periodic_utility_functions.R
# to fit periodic time, stochastic time, and homogeneous time dynamic
# occupancy models to 9 seasons of Chicago camera trap data
# for coyote, red fox, striped skunk, raccoon, and virginia opossum.
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
z <- df_2_array(read.table(
  "Fidino_2017_community_incidence_matrix_sp11_sp13.txt", 
  header = TRUE, sep = "\t"))[1,,]

# proportion of sites occupied each time step
prop_sites_occupied <- colSums(z, na.rm = TRUE) / 
  apply(z, 2, function(x) 95 - sum(is.na(x)))

# get beta a and b from mean and standard deviation of soc
prior_for_occ <- betaABfromMeanSD(mean(prop_sites_occupied), 
  sd(prop_sites_occupied))

covdat <- read.table("Fidino_2017_URB_covariate.txt", header = TRUE, sep ="\t")

# read in the y array (days a species was seen per site and season)
# ordered identically to community incidence matrix z
y_array <- df_2_array(read.table(
  "Fidino_2017_chicago_detection_data_sp11_sp13.txt", 
  header = TRUE, sep = "\t"))

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
                  x = covdat$x, P = 4)


# models to fit to the data
models <- c("Fidino_2017_periodic_time_model_pulse.R",
            "Fidino_2017_stochastic_time_model.R",
            "Fidino_2017_homogeneous_time_model.R")

# paramaters to track from each model
to_monitor <- list( c("psinit", "Gt","Gt_sd", "d0","d1","m0",
                      "m1", "A", "delta", "omega", "omega_sd",
                      "epsilon", "epsilon_sd","f0", "f1", "z", "p"),
                    c("psinit", "Gt","Gt_sd", "d0","d1","m0",
                      "m1", "Ut", "Ut_sd", "omega", "omega_sd",
                      "epsilon", "epsilon_sd","f0", "f1", "z", "p"),
                    c("psinit", "Gt","Gt_sd", "d0","d1","m0",
                      "m1", "omega", "omega_sd",
                      "epsilon", "epsilon_sd","f0", "f1", "z", "p"))

# functions that generate initial values for each model
init_functions <- list(inits_PTM_pulse, 
                       inits_STM,
                       inits_HTM)

# fit each model, calculate CPO for each
coyote_scores <- fit_models(models = models,
                            data_list = data_list,
                            init_function = init_functions,
                            to_monitor = to_monitor,
                            species = "coyote")
# save CPO scores for each model 
write.table(coyote_scores, "coyote_scores.txt", row.names = FALSE, sep = "\t")

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
                  x = covdat$x, P = 4)

fox_scores <- fit_models(models = models,
                            data_list = data_list,
                            init_function = init_functions,
                            to_monitor = to_monitor,
                            species = "redfox")
write.table(fox_scores, "fox_scores.txt", row.names = FALSE, sep = "\t")

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
                  x = covdat$x, P = 4)

skunk_scores <- fit_models(models = models,
                         data_list = data_list,
                         init_function = init_functions,
                         to_monitor = to_monitor,
                         species = "skunk")
write.table(skunk_scores, "striped_skunk_scores.txt", 
  row.names = FALSE, sep = "\t")

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

models <- c("Fidino_2017_periodic_time_model_boom_bust.R",
            "Fidino_2017_stochastic_time_model.R",
            "Fidino_2017_homogeneous_time_model.R")

init_functions <- list(inits_PTM_boom, 
                       inits_STM,
                       inits_HTM)

data_list <- list(y = as.matrix(y_array[4,,]), nyear = ncol(z), 
                  nsite = nrow(z), 
                  spa = prior_for_occ$a, spb = prior_for_occ$b,
                  jmat = as.matrix(j_mat),
                  pi = 3.14159, C = cs[[1]], S = cs[[2]],
                  x = covdat$x, P = 2)



raccoon_scores <- fit_models(models = models,
                             data_list = data_list,
                             init_function = init_functions,
                             to_monitor = to_monitor,
                             species = "raccoon")

write.table(raccoon_scores, "raccoon_scores.txt", 
  row.names = FALSE, sep = "\t")
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
                  x = covdat$x, P = 2)

opossum_scores <- fit_models(models = models,
                           data_list = data_list,
                           init_function = init_functions,
                           to_monitor = to_monitor,
                           species = "opossum")

write.table(opossum_scores, "opossum_scores.txt", row.names = FALSE, sep = "\t")


