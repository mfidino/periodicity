


source("periodicity_utility.R")

#####################################################################
# Assuming the files are saved in your working directory and you have
### sourced the utility script
#####################################################################

# load the reshape package
package_load(c("reshape2", "runjags", "rjags", "parallel"))


# lets do the coyote
z <- df_2_array(read.table("z_matrix_sp10_sp13.txt", header = TRUE, sep = "\t"))[1,,]

# get number of sites occupied


soc <- colSums(z, na.rm = TRUE) / apply(z, 2, function(x) 100 - sum(is.na(x)))

# get beta a and b from mean and standard deviation of soc

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
j_mat[68,3:4] <- 26




params <- c("lp", "ly", "g0", "p0", "gy", "py", "gy_sd",
            "py_sd", "ly_sd", "psinit")

# coyote analysis
cs <- make_c_s_mat(1:12, 4)

data_list_trig <- list(y = as.matrix(y_array[which(species_names$x=="Coyote"),,]), nyear = ncol(z), 
                  nsite = nrow(z), 
                  spa = prior_for_occ$a, spb = prior_for_occ$b,
                  jmat = as.matrix(j_mat),
                  pi = 3.14159, C = cs[[1]], S = cs[[2]],
                  cov = covdat$pc1, P = 4)


models <- c("ranef_year_jags.R", 
            "pulse_year_jags_trig.R",
            "only_pulse_year_jags_trig.R",
            "homog_time.R")

to_mon <- list(c("psinit", "gy", "py", "ly", "ls", "g_mu", "p_mu",
                     "lp", "g1", "p1", "l1", "gy_sd", "py_sd", "ly_sd",
                     "ls_sd", "y_pred", "z", "p"),
                   c("psinit", "ly", "ls", "g_mu", "p_mu",
                     "lp", "g1", "p1", "l1", "gy_sd", "py_sd", "ly_sd",
                     "ls_sd", "y_pred", "theta", "dp", "z", "p"),
                   c("psinit", "gy", "py", "ly", "ls", "g_mu", "p_mu",
                     "lp", "g1", "p1", "l1", "ly_sd",
                     "ls_sd", "y_pred", "theta", "dp", "z","p" ),
                   c("psinit", "py", "ly", "ls", "g_mu", "p_mu",
                     "lp", "g1", "p1", "l1", "py_sd", "ly_sd",
                     "ls_sd", "y_pred","z","p"))
species <- "coyote"
inl <- list(inits_ranef, inits_pulse, 
            inits_only_pulse, inits_homog)

# the model outputs will be saved
coyote_scores <- fit_models(models, data_list_trig, inl, to_mon, "coyote")

coyote_scores <- data.frame(loss_scores = c(1346.364, 1334.906, 1338.537, 1345.094),
                            species = "coyote",
                            model = c("ranef_year_jags", "pulse_year_jags_trig",
                                      "only_pulse_year_jags_trig", "homog_time"))
write.table(coyote_scores, "./model_outputs/coyote_scores.txt", row.names = FALSE, sep = "\t")

### red fox
z <- df_2_array(read.table("z_matrix_sp10_sp13.txt", header = TRUE, sep = "\t"))[5,,]

# get number of sites occupied


soc <- colSums(z, na.rm = TRUE) / apply(z, 2, function(x) 100 - sum(is.na(x)))

# get beta a and b from mean and standard deviation of soc

prior_for_occ <- betaABfromMeanSD(mean(soc), sd(soc))

data_list_trig <- list(y = as.matrix(y_array[which(species_names$x=="Redfox"),,]), nyear = ncol(z), 
                       nsite = nrow(z), 
                       spa = prior_for_occ$a, spb = prior_for_occ$b,
                       jmat = as.matrix(j_mat),
                       pi = 3.14159, C = cs[[1]], S = cs[[2]],
                       cov = covdat$pc1, P = 4)

# the model outputs will be saved
fox_scores <- fit_models(models, data_list_trig, inl, to_mon, "redfox")
write.table(fox_scores, "./model_outputs/fox_scores.txt", row.names = FALSE, sep = "\t")
# now striped skunk

z <- df_2_array(read.table("z_matrix_sp10_sp13.txt", header = TRUE, sep = "\t"))[6,,]

# get number of sites occupied


soc <- colSums(z, na.rm = TRUE) / apply(z, 2, function(x) 100 - sum(is.na(x)))

# get beta a and b from mean and standard deviation of soc

prior_for_occ <- betaABfromMeanSD(mean(soc), sd(soc))

data_list_trig <- list(y = as.matrix(y_array[which(species_names$x=="Skunk"),,]), nyear = ncol(z), 
                       nsite = nrow(z), 
                       spa = prior_for_occ$a, spb = prior_for_occ$b,
                       jmat = as.matrix(j_mat),
                       pi = 3.14159, C = cs[[1]], S = cs[[2]],
                       cov = covdat$pc1, P = 4)

# the model outputs will be saved
skunk_scores <- fit_models(models, data_list_trig, inl, to_mon, "skunk")
write.table(skunk_scores, "./model_outputs/skunk_scores.txt", row.names = FALSE, sep = "\t")
#######################################################
###### Opossum
#######################################################

z <- df_2_array(read.table("z_matrix_sp10_sp13.txt", header = TRUE, sep = "\t"))[3,,]

# get number of sites occupied


soc <- colSums(z, na.rm = TRUE) / apply(z, 2, function(x) 100 - sum(is.na(x)))

# get beta a and b from mean and standard deviation of soc

prior_for_occ <- betaABfromMeanSD(mean(soc), sd(soc))

cs <- make_c_s(1:12, 2)

data_list_boom <- list(y = as.matrix(y_array[which(species_names$x=="Opossum"),,]), nyear = ncol(z), 
                       nsite = nrow(z), 
                       spa = prior_for_occ$a, spb = prior_for_occ$b,
                       jmat = as.matrix(j_mat),
                       pi = 3.14159, C = cs[[1]], S = cs[[2]],
                       cov = covdat$pc1, P = 2)


models <- c("ranef_year_jags.R", 
            "boom_bust_jags_trig.R",
            "only_boom_bust_jags_trig.R",
            "homog_time.R")

to_mon <- list(c("psinit", "gy", "py", "ly", "ls", "g_mu", "p_mu",
                 "lp", "g1", "p1", "l1", "gy_sd", "py_sd", "ly_sd",
                 "ls_sd", "y_pred", "z", "p"),
               c("psinit", "ly", "ls", "g_mu", "p_mu",
                 "lp", "g1", "p1", "l1", "gy_sd", "py_sd", "ly_sd",
                 "ls_sd", "y_pred", "a1_gam", "a2_gam", "z", "p"),
               c("psinit", "gy", "py", "ly", "ls", "g_mu", "p_mu",
                 "lp", "g1", "p1", "l1", "ly_sd",
                 "ls_sd", "y_pred", "a1_gam", "a2_gam", "z", "p"),
               c("psinit", "py", "ly", "ls", "g_mu", "p_mu",
                 "lp", "g1", "p1", "l1", "py_sd", "ly_sd",
                 "ls_sd", "y_pred","z", "p"))
species <- "opossum"
inl <- list(inits_ranef, inits_boom, inits_only_boom, inits_homog)

# the model outputs will be saved
opossum_scores <- fit_models(models, data_list_boom, inl, to_mon, "opossum")
write.table(opossum_scores, "./model_outputs/opossum_scores.txt", row.names = FALSE, sep = "\t")
###############################
#raccoon
##################################

z <- df_2_array(read.table("z_matrix_sp10_sp13.txt", header = TRUE, sep = "\t"))[4,,]
prior_for_occ <- betaABfromMeanSD(mean(soc), sd(soc))


data_list_boom <- list(y = as.matrix(y_array[which(species_names$x=="Raccoon"),,]), nyear = ncol(z), 
                       nsite = nrow(z), 
                       spa = prior_for_occ$a, spb = prior_for_occ$b,
                       jmat = as.matrix(j_mat),
                       pi = 3.14159, C = cs[[1]], S = cs[[2]],
                       cov = covdat$pc1, P = 2)

raccoon_scores <- fit_models(models, data_list_boom, inl, to_mon, "raccoon")
write.table(raccoon_scores, "./model_outputs/raccoon_scores.txt", row.names = FALSE, sep = "\t")
all_scores <- rbind(coyote_scores, fox_scores, skunk_scores, opossum_scores,raccoon_scores)
write.table(all_scores, "./model_outputs/all_scores.txt", row.names = FALSE, sep = "\t")
# swtich things around to do the opossum and raccoon

py1 <- mase(puly, data_list, type = "naive")
py2 <- mase(puly, data_list, type = "season")

all_ans <- cbind(ry1, py1, ry2, py2)
apply(all_ans, 2, quantile)
pz1 <- mase(pulsez, data_list, type = "naive")
pz2 <- mase(pulsez, data_list, type = "season")
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
