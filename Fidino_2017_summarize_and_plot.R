############################
#
# Summary and Plotting script for Fidino and Magle 2017:
# Using Fourier series to estimate periodic
# trends in dynamic occupancy models.
#
#

### DESCRIPTION ###### DESCRIPTION ###### DESCRIPTION ###### DESCRIPTION ###
# This script uses the functions in Fidino_2017_periodic_utility_functions.R
# to summarize periodic time, stochastic time, and homogeneous time dynamic
# occupancy models to 9 seasons of Chicago camera trap data
# for coyote, red fox, striped skunk, virginia opossum, and raccoon.
### DESCRIPTION ###### DESCRIPTION ###### DESCRIPTION ###### DESCRIPTION ###


# coyote

# PTM
species_posterior <- fread("Fidino_2017_periodic_time_model_pulse_coyote_post.txt", data.table = FALSE) 

# read in observed data
z <- df_2_array(read.table("Fidino_2017_community_incidence_matrix_sp11_sp13.txt", header = TRUE, sep = "\t"))[1,,]

# get estimated z from species_posterior
zest <- species_posterior[,grep("z", colnames(species_posterior))]

# calculate most likely value for delta
delta <- species_posterior[,grep("delta", colnames(species_posterior))]
# save delta
write.table(prop.table(table(delta)), "Fidino_2017_coyote_when_pulse.txt", row.names = FALSE)

# summarise and save posterior
write.table(round(t(apply(species_posterior, 2, quantile, 
                          probs = c(0.025,0.5,0.975))),2),
            "Fidino_2017_coyote_PTM_summary.txt", sep = "\t")

# get observed colonization, 
observed_colonization <- matrix(0, ncol = 8, nrow = 5)
observed_colonization[1,] <- zest_observed(z)
 
# convert zest to 3d array, steps by site by session
zest_array <- array(as.matrix(zest), dim = c(nrow(zest), 95, 9))

# array for zest_posterior
zest_PTM <- array(0, dim = c(3,8,5))
zest_PTM[,,1] <- zest_posterior(zest_array)

#make colonization predictions
col_predict_PTM <-  array(0, dim=c(3, 8, 5))
col_predict_PTM[,,1] <- make_pred_PTM_pulse(species_posterior, 2, 4,1:8)

rm(zest_array, species_posterior)

# STM model
species_posterior <- fread("Fidino_2017_stochastic_time_model_coyote_post.txt", data.table = FALSE)

# get estimated z from species_posterior
zest <- species_posterior[,grep("z", colnames(species_posterior))]

# summarise and save posterior
write.table(round(t(apply(species_posterior, 2, quantile, 
                          probs = c(0.025,0.5,0.975))),2),
            "Fidino_2017_coyote_STM_summary.txt", sep = "\t")

# convert zest to 3d array, steps by site by session
zest_array <- array(as.matrix(zest), dim = c(nrow(zest), 95, 9))

# array for zest_posterior ranef
zest_STM <- array(0, dim = c(3,8,5))
zest_STM[,,1] <- zest_posterior(zest_array)

col_predict_STM <-  array(0, dim = c(5,8,5))
col_predict_STM[,,1] <- make_pred_STM(species_posterior)

rm(species_posterior, zest_array)

# HTM
species_posterior <- fread("Fidino_2017_homogeneous_time_model_coyote_post.txt", data.table = FALSE)

# get estimated z from species_posterior
zest <- species_posterior[,grep("z", colnames(species_posterior))]

# summarise and save posterior
write.table(round(t(apply(species_posterior, 2, quantile, 
                          probs = c(0.025,0.5,0.975))),2),
            "Fidino_2017_coyote_HTM_summary.txt", sep = "\t")

# convert zest to 3d array, steps by site by session
zest_array <- array(as.matrix(zest), dim = c(nrow(zest), 95, 9))

# array for zest_posterior ranef
zest_HTM <- array(0, dim = c(3,8,5))
zest_HTM[,,1] <- zest_posterior(zest_array)

col_predict_HTM <-  array(0, dim = c(3,8,5))
col_predict_HTM[,,1] <- make_pred_HTM(species_posterior, 8)

rm(species_posterior, zest_array)

#######################red fox

# PTM
species_posterior <- fread("Fidino_2017_periodic_time_model_pulse_redfox_post.txt", data.table = FALSE) 

# read in observed data
z <- df_2_array(read.table("Fidino_2017_community_incidence_matrix_sp11_sp13.txt", header = TRUE, sep = "\t"))[2,,]

# get estimated z from species_posterior
zest <- species_posterior[,grep("z", colnames(species_posterior))]

# calculate most likely value for delta
delta <- species_posterior[,grep("delta", colnames(species_posterior))]
# save delta
write.table(prop.table(table(delta)), "Fidino_2017_redfox_when_pulse.txt", row.names = FALSE)

# summarise and save posterior
write.table(round(t(apply(species_posterior, 2, quantile, 
                          probs = c(0.025,0.5,0.975))),2),
            "Fidino_2017_redfox_PTM_summary.txt", sep = "\t")

# get observed colonization, 
observed_colonization[2,] <- zest_observed(z)

# convert zest to 3d array, steps by site by session
zest_array <- array(as.matrix(zest), dim = c(nrow(zest), 95, 9))

# array for zest_posterior
zest_PTM[,,2] <- zest_posterior(zest_array)

#make colonization predictions
col_predict_PTM[,,2] <- make_pred_PTM_pulse(species_posterior, 2, 4,1:8)

rm(zest_array, species_posterior)

# STM model
species_posterior <- fread("Fidino_2017_stochastic_time_model_redfox_post.txt", data.table = FALSE)

# get estimated z from species_posterior
zest <- species_posterior[,grep("z", colnames(species_posterior))]

# summarise and save posterior
write.table(round(t(apply(species_posterior, 2, quantile, 
                          probs = c(0.025,0.5,0.975))),2),
            "Fidino_2017_redfox_STM_summary.txt", sep = "\t")

# convert zest to 3d array, steps by site by session
zest_array <- array(as.matrix(zest), dim = c(nrow(zest), 95, 9))

# array for zest_posterior ranef
zest_STM[,,2] <- zest_posterior(zest_array)

col_predict_STM[,,2] <- make_pred_STM(species_posterior)

rm(species_posterior, zest_array)

# HTM
species_posterior <- fread("Fidino_2017_homogeneous_time_model_redfox_post.txt", data.table = FALSE)

# get estimated z from species_posterior
zest <- species_posterior[,grep("z", colnames(species_posterior))]

# summarise and save posterior
write.table(round(t(apply(species_posterior, 2, quantile, 
                          probs = c(0.025,0.5,0.975))),2),
            "Fidino_2017_redfox_HTM_summary.txt", sep = "\t")

# convert zest to 3d array, steps by site by session
zest_array <- array(as.matrix(zest), dim = c(nrow(zest), 95, 9))

# array for zest_posterior ranef
zest_HTM[,,2] <- zest_posterior(zest_array)

# make predictions
col_predict_HTM[,,2] <- make_pred_HTM(species_posterior, 8)

rm(species_posterior, zest_array)

######################### skunk

# PTM
species_posterior <- fread("Fidino_2017_periodic_time_model_pulse_skunk_post.txt", data.table = FALSE) 

# read in observed data
z <- df_2_array(read.table("Fidino_2017_community_incidence_matrix_sp11_sp13.txt", header = TRUE, sep = "\t"))[3,,]

# get estimated z from species_posterior
zest <- species_posterior[,grep("z", colnames(species_posterior))]

# calculate most likely value for delta
delta <- species_posterior[,grep("delta", colnames(species_posterior))]
# save delta
write.table(prop.table(table(delta)), "Fidino_2017_skunk_when_pulse.txt", row.names = FALSE)

# summarise and save posterior
write.table(round(t(apply(species_posterior, 2, quantile, 
                          probs = c(0.025,0.5,0.975))),2),
            "Fidino_2017_skunk_PTM_summary.txt", sep = "\t")

# get observed colonization, 
observed_colonization[3,] <- zest_observed(z)

# convert zest to 3d array, steps by site by session
zest_array <- array(as.matrix(zest), dim = c(nrow(zest), 95, 9))

# array for zest_posterior
zest_PTM[,,3] <- zest_posterior(zest_array)

#make colonization predictions
col_predict_PTM[,,3] <- make_pred_PTM_pulse(species_posterior, 2, 4,1:8)

rm(zest_array, species_posterior)

# STM model
species_posterior <- fread("Fidino_2017_stochastic_time_model_skunk_post.txt", data.table = FALSE)

# get estimated z from species_posterior
zest <- species_posterior[,grep("z", colnames(species_posterior))]

# summarise and save posterior
write.table(round(t(apply(species_posterior, 2, quantile, 
                          probs = c(0.025,0.5,0.975))),2),
            "Fidino_2017_skunk_STM_summary.txt", sep = "\t")

# convert zest to 3d array, steps by site by session
zest_array <- array(as.matrix(zest), dim = c(nrow(zest), 95, 9))

# array for zest_posterior ranef
zest_STM[,,3] <- zest_posterior(zest_array)

col_predict_STM[,,3] <- make_pred_STM(species_posterior)

rm(species_posterior, zest_array)

# HTM
species_posterior <- fread("Fidino_2017_homogeneous_time_model_skunk_post.txt", data.table = FALSE)

# get estimated z from species_posterior
zest <- species_posterior[,grep("z", colnames(species_posterior))]

# summarise and save posterior
write.table(round(t(apply(species_posterior, 2, quantile, 
                          probs = c(0.025,0.5,0.975))),2),
            "Fidino_2017_skunk_HTM_summary.txt", sep = "\t")

# convert zest to 3d array, steps by site by session
zest_array <- array(as.matrix(zest), dim = c(nrow(zest), 95, 9))

# array for zest_posterior ranef
zest_HTM[,,3] <- zest_posterior(zest_array)

# make predictions
col_predict_HTM[,,3] <- make_pred_HTM(species_posterior, 8)

rm(species_posterior, zest_array)

################# raccoon

# PTM
species_posterior <- fread("Fidino_2017_periodic_time_model_boom_bust_raccoon_post.txt", data.table = FALSE) 

# read in observed data
z <- df_2_array(read.table("Fidino_2017_community_incidence_matrix_sp11_sp13.txt", header = TRUE, sep = "\t"))[4,,]

# get estimated z from species_posterior
zest <- species_posterior[,grep("z", colnames(species_posterior))]

# calculate most likely value for delta
delta <- species_posterior[,grep("delta", colnames(species_posterior))]
# save delta
write.table(prop.table(table(delta)), "Fidino_2017_raccoon_when_pulse.txt", row.names = FALSE)

# summarise and save posterior
write.table(round(t(apply(species_posterior, 2, quantile, 
                          probs = c(0.025,0.5,0.975))),2),
            "Fidino_2017_raccoon_PTM_summary.txt", sep = "\t")

# get observed colonization, 
observed_colonization[4,] <- zest_observed(z)

# convert zest to 3d array, steps by site by session
zest_array <- array(as.matrix(zest), dim = c(nrow(zest), 95, 9))

# array for zest_posterior
zest_PTM[,,4] <- zest_posterior(zest_array)

#make colonization predictions
col_predict_PTM[,,4] <- make_pred_PTM_boom(species_posterior, 0, 2,1:8)

rm(zest_array, species_posterior)

# STM model
species_posterior <- fread("Fidino_2017_stochastic_time_model_raccoon_post.txt", data.table = FALSE)


# get estimated z from species_posterior
zest <- species_posterior[,grep("z", colnames(species_posterior))]

# summarise and save posterior
write.table(round(t(apply(species_posterior, 2, quantile, 
                          probs = c(0.025,0.5,0.975))),2),
            "Fidino_2017_raccoon_STM_summary.txt", sep = "\t")

# convert zest to 3d array, steps by site by session
zest_array <- array(as.matrix(zest), dim = c(nrow(zest), 95, 9))

# array for zest_posterior ranef
zest_STM[,,4] <- zest_posterior(zest_array)

col_predict_STM[,,4] <- make_pred_STM(species_posterior)

rm(species_posterior, zest_array)

# HTM
species_posterior <- fread("Fidino_2017_homogeneous_time_model_raccoon_post.txt", data.table = FALSE)

# get estimated z from species_posterior
zest <- species_posterior[,grep("z", colnames(species_posterior))]

# summarise and save posterior
write.table(round(t(apply(species_posterior, 2, quantile, 
                          probs = c(0.025,0.5,0.975))),2),
            "Fidino_2017_raccoon_HTM_summary.txt", sep = "\t")

# convert zest to 3d array, steps by site by session
zest_array <- array(as.matrix(zest), dim = c(nrow(zest), 95, 9))

# array for zest_posterior ranef
zest_HTM[,,4] <- zest_posterior(zest_array)

# make predictions
col_predict_HTM[,,4] <- make_pred_HTM(species_posterior, 8)

rm(species_posterior, zest_array)

################# opossum

# PTM
species_posterior <- fread("Fidino_2017_periodic_time_model_boom_bust_opossum_post.txt", data.table = FALSE) 

# read in observed data
z <- df_2_array(read.table("Fidino_2017_community_incidence_matrix_sp11_sp13.txt", header = TRUE, sep = "\t"))[5,,]

# get estimated z from species_posterior
zest <- species_posterior[,grep("z", colnames(species_posterior))]

# calculate most likely value for delta
delta <- species_posterior[,grep("delta", colnames(species_posterior))]
# save delta
write.table(prop.table(table(delta)), "Fidino_2017_opossum_when_pulse.txt", row.names = FALSE)

# summarise and save posterior
write.table(round(t(apply(species_posterior, 2, quantile, 
                          probs = c(0.025,0.5,0.975))),2),
            "Fidino_2017_opossum_PTM_summary.txt", sep = "\t")

# get observed colonization, 
observed_colonization[5,] <- zest_observed(z)

# convert zest to 3d array, steps by site by session
zest_array <- array(as.matrix(zest), dim = c(nrow(zest), 95, 9))

# array for zest_posterior
zest_PTM[,,5] <- zest_posterior(zest_array)

#make colonization predictions
col_predict_PTM[,,5] <- make_pred_PTM_boom(species_posterior, 0, 2,1:8)

rm(zest_array, species_posterior)

# STM model
species_posterior <- fread("Fidino_2017_stochastic_time_model_opossum_post.txt", data.table = FALSE)

# get estimated z from species_posterior
zest <- species_posterior[,grep("z", colnames(species_posterior))]

# summarise and save posterior
write.table(round(t(apply(species_posterior, 2, quantile, 
                          probs = c(0.025,0.5,0.975))),2),
            "Fidino_2017_opossum_STM_summary.txt", sep = "\t")

# convert zest to 3d array, steps by site by session
zest_array <- array(as.matrix(zest), dim = c(nrow(zest), 95, 9))

# array for zest_posterior ranef
zest_STM[,,5] <- zest_posterior(zest_array)

col_predict_STM[,,5] <- make_pred_STM(species_posterior)

rm(species_posterior, zest_array)

# HTM
species_posterior <- fread("Fidino_2017_homogeneous_time_model_opossum_post.txt", data.table = FALSE)

# get estimated z from species_posterior
zest <- species_posterior[,grep("z", colnames(species_posterior))]

# summarise and save posterior
write.table(round(t(apply(species_posterior, 2, quantile, 
                          probs = c(0.025,0.5,0.975))),2),
            "Fidino_2017_opossum_HTM_summary.txt", sep = "\t")

# convert zest to 3d array, steps by site by session
zest_array <- array(as.matrix(zest), dim = c(nrow(zest), 95, 9))

# array for zest_posterior ranef
zest_HTM[,,5] <- zest_posterior(zest_array)

# make predictions
col_predict_HTM[,,5] <- make_pred_HTM(species_posterior, 8)

rm(species_posterior, zest_array)

# write the remaining files
saveRDS(col_predict_PTM,"./model_summary/col_predict_PTM.RDS")
saveRDS(col_predict_STM,"./model_summary/col_predict_STM.RDS")
saveRDS(col_predict_HTM,"./model_summary/col_predict_HTM.RDS")

saveRDS(zest_PTM, "./model_summary/zest_PTM.RDS")
saveRDS(zest_HTM, "./model_summary/zest_HTM.RDS")
saveRDS(zest_STM, "./model_summary/zest_STM.RDS")

saveRDS(observed_colonization, "./model_summary/observed_colonization.RDS")
#####################
### actual plotting #
#####################

col_predict_HTM <- readRDS("./model_summary/col_predict_HTM.RDS")
col_predict_PTM <- readRDS("./model_summary/col_predict_PTM.RDS")
col_predict_STM <- readRDS("./model_summary/col_predict_STM.RDS")

zest_HTM <- readRDS("./model_summary/zest_HTM.RDS")
zest_PTM <- readRDS("./model_summary/zest_PTM.RDS")
zest_STM <- readRDS("./model_summary/zest_STM.RDS")

observed_colonization <- readRDS("./model_summary/observed_colonization.RDS")

for(pup in 1:1){
tiff("all_col5.tiff", height = 14, width = 12, units = "in",
     compression = "lzw", res = 600)
#windows(12, 14)
#tiff("periodic_plot_3p2.tiff", height = 6, width = 7.5,
 #    units = "in", compression = "lzw", res = 600)
m <- matrix(c(1:15), ncol = 3, nrow = 5 )
#m <- matrix(c(1:6), ncol = 2, nrow = 3 )
layout(m)
par( mar = c(2,2.9, 1.8,0.5),
     oma = c(2,2.9,1.8,0.5) + 0.1)

#lims <- c(1,0.75,0.75, 1, 0.75)
lims <- c(1,1,1, 1, 1)

nmsa <- c(1.2, 1.25, 1.85, 1.3, 1.35)

pls_letters <- paste0(LETTERS[seq(1, 13, 3)], ")")
for(i in 1:5){
  plot(1~1, type = "n", xlim = c(1,8), ylim = c(0,lims[i]), xlab = "",
       ylab = "", xaxt = "n", yaxt="n", bty = "n")
  
  axis(1, at= c(1:8), labels = F, tck = -.05)
  se <- c("SU", "FA", "WI", "SP")
  mtext("Year:", 1, line = 1.2, at = 0.1, cex = 1.2)
  mtext("Season:", 1, line = 2.6, at = -.18, cex = 1.2)
  mtext(text = c("'11", "'12", "'13"), 1, line = 1.2, at = c(1,4,8), cex = 1.2)
  mtext(text = rep(se, 2),1, line = 2.6, at = c(1:8), cex = 1.2)
  
  
  axis(2, at=seq(0,lims[i], 0.25), labels=F, tck=-.05)
  axis(2, at=seq(0,lims[i], 0.125), labels=F, tck=-.03)
  mtext(text = as.character(seq(0,lims[i],0.25)), 2, line = 1.2, at = seq(0,lims[i],0.25),
        las = 1, cex = 1.1)
  
  mtext(bquote(paste("Pr(", italic(gamma), ")", sep = "")),
        2, at = lims[i]/2, line = 3.9, cex = 1.4)
  #mtext(nms[i],3, at = nmsa[i], line = 0.4, cex = 1)
  
  #mtext("B) Red fox",3, at = 1.25, line = 0.1, cex = 1.2)
  
  y <- col_predict_PTM[1,,i]
  y2 <- rev(col_predict_PTM[3,,i])
  x <- 1:8
  x2 <- rev(x)
  polygon(c(x, x2), c(y, y2), col = "azure3", border = FALSE)
  lines(col_predict_PTM[2,,i], lwd = 2)
  #lines(mry[2,,i], lty = 3, lwd = 0.75)
  #lines(ary[1,,i], lwd = 2, lty = 2)
  #lines(ary[3,,i], lwd = 2, lty = 2)
  points(1:8, observed_colonization[i,], pch = 21, bg = "white", col = "black", cex = 1.5)
  for(j in 1:8){
    lines(c(j,j), zest_PTM[-2,j,i], lwd = 1.25)
  }
  points(zest_PTM[2,,i], pch = 21, bg = "black",col = "black", cex = 1.5)
  legend("topleft", pls_letters[i], bty = "n", cex = 1.6, inset = c(-0.05, 0))
}
ran_letters <- paste0(LETTERS[seq(2, 16, 3)], ")")
# plotting for ranef

for(i in 1:5){
  plot(1~1, type = "n", xlim = c(1,8), ylim = c(0,lims[i]), xlab = "",
       ylab = "", xaxt = "n", yaxt="n", bty = "n")
  
  axis(1, at= c(1:8), labels = F, tck = -.05)
  se <- c("SU", "FA", "WI", "SP")
  mtext(text = c("'11", "'12", "'13"), 1, line = 1.2, at = c(1,4,8), cex = 1.2)
  mtext(text = rep(se, 2),1, line = 2.6, at = c(1:8), cex = 1.2)
  
  
  axis(2, at=seq(0,lims[i], 0.25), labels=F, tck=-.05)
  axis(2, at=seq(0,lims[i], 0.125), labels=F, tck=-.03)
  mtext(text = as.character(seq(0,lims[i],0.25)), 2, line = 1.2, at = seq(0,lims[i],0.25),
        las = 1, cex = 1.1)
  
  
  y <- col_predict_STM[2,,i]
  y2 <- rev(col_predict_STM[4,,i])
  x <- 1:8
  x2 <- rev(x)
  polygon(c(x, x2), c(y, y2), col = "azure3", border = FALSE)
  lines(col_predict_STM[3,,i], lwd = 2)
  lines(col_predict_STM[1,,i], lty = 2)
  lines(col_predict_STM[5,,i], lty = 2)

  points(1:8, observed_colonization[i,], pch = 21, bg = "white", col = "black", cex = 1.5)
  for(j in 1:8){
    lines(c(j,j), zest_STM[-2,j,i], lwd = 1.25)
  }
  points(zest_STM[2,,i], pch = 21, bg = "black",col = "black", cex = 1.5)
  legend("topleft", ran_letters[i], bty = "n", cex = 1.6, inset = c(-0.05, 0))
}
homog_letters <- paste0(LETTERS[seq(3, 17, 3)], ")")
for(i in 1:5){
  plot(1~1, type = "n", xlim = c(1,8), ylim = c(0,lims[i]), xlab = "",
       ylab = "", xaxt = "n", yaxt="n", bty = "n")
  
  axis(1, at= c(1:8), labels = F, tck = -.05)
  se <- c("SU", "FA", "WI", "SP")
  #mtext("Year:", 1, line = 0.55, at = 0, cex = 0.6)
  #mtext("Season:", 1, line = 1.45, at = -.275, cex = 0.6)
  mtext(text = c("'11", "'12", "'13"), 1, line = 1.2, at = c(1,4,8), cex = 1.2)
  mtext(text = rep(se, 2),1, line = 2.6, at = c(1:8), cex = 1.2)
  
  
  axis(2, at=seq(0,lims[i], 0.25), labels=F, tck=-.05)
  axis(2, at=seq(0,lims[i], 0.125), labels=F, tck=-.03)
  mtext(text = as.character(seq(0,lims[i],0.25)), 2, line = 1.2, at = seq(0,lims[i],0.25),
        las = 1, cex = 1.1)
  

  
  y <- col_predict_HTM[1,,i]
  y2 <- rev(col_predict_HTM[3,,i])
  x <- 1:8
  x2 <- rev(x)
  polygon(c(x, x2), c(y, y2), col = "azure3", border = FALSE)
  lines(col_predict_HTM[2,,i], lwd = 2)

  points(1:8, observed_colonization[i,], pch = 21, bg = "white", col = "black", cex = 1.5)
  for(j in 1:8){
    lines(c(j,j), zest_HTM[-2,j,i], lwd = 1.25)
  }
  points(zest_HTM[2,,i], pch = 21, bg = "black",col = "black", cex = 1.5)
  legend("topleft", homog_letters[i], bty = "n", cex = 1.6, inset = c(-0.05, 0))
}
dev.off()
}


text()






  tiff("coyote_test2.tiff", height = 3, width = 3.5, units = "in",
       compression = "lzw", res = 600)
  windows(7.5, 6)
  tiff("periodic_plot_3p2.tiff", height = 6, width = 7.5,
       units = "in", compression = "lzw", res = 600)
  m <- matrix(c(1:6), ncol = 2, nrow = 3 )
  layout(m)
  par( mar = c(3,4, 2.25,1.25),
       oma = c(3,4,2.25,1.25) + 0.1)
  
  lims <- c(0.75,0.5,0.5, 1, 0.75)
  nms <- c("A) Coyote", "B) Red fox", "C) Striped skunk",
           "D) Raccoon","E) Opossum")
  nmsa <- c(1.2, 1.25, 1.85, 1.3, 1.35)
  
  for(i in 1:5){
  plot(1~1, type = "n", xlim = c(1,8), ylim = c(0,lims[i]), xlab = "",
       ylab = "", xaxt = "n", yaxt="n", bty = "n")
  
  axis(1, at= c(1:8), labels = F, tck = -.05)
  se <- c("SU", "FA", "WI", "SP")
  mtext("Year:", 1, line = 0.55, at = 0, cex = 0.6)
  mtext("Season:", 1, line = 1.45, at = -.275, cex = 0.6)
  mtext(text = c("'11", "'12", "'13"), 1, line = 0.6, at = c(1,4,8), cex = 0.75)
  mtext(text = rep(se, 2),1, line = 1.5, at = c(1:8), cex = 0.75)
  
  
  axis(2, at=seq(0,lims[i], 0.25), labels=F, tck=-.05)
  axis(2, at=seq(0,lims[i], 0.125), labels=F, tck=-.03)
  mtext(text = as.character(seq(0,lims[i],0.25)), 2, line = 0.6, at = seq(0,lims[i],0.25),
        las = 1)
  
  mtext(bquote(paste("Pr(", italic(gamma), ")", sep = "")),
        2, at = lims[i]/2, line = 3.2, cex = 1)
  mtext(nms[i],3, at = nmsa[i], line = 0.4, cex = 1)
  
  #mtext("B) Red fox",3, at = 1.25, line = 0.1, cex = 1.2)

  y <- ary[1,,i]
  y2 <- rev(ary[3,,i])
  x <- 1:8
  x2 <- rev(x)
  polygon(c(x, x2), c(y, y2), col = "azure3", border = FALSE)
  lines(ary[2,,i], lwd = 2)
  #lines(mry[2,,i], lty = 3, lwd = 0.75)
  #lines(ary[1,,i], lwd = 2, lty = 2)
  #lines(ary[3,,i], lwd = 2, lty = 2)
  points(1:8, kry[i,], pch = 21, bg = "white", col = "black", cex = 1.25)
    for(j in 1:8){
    lines(c(j,j), mry[-2,j,i], lwd = 1.25)
  }
  points(mry[2,,i], pch = 21, bg = "black",col = "black", cex = 1.25)
  }

  
  #########################
  dev.off()


abline(v = 2, lty = 2)
abline(h = -0.2, lty = 2)
lines(x = c(2,3.2), y = c(0.653, 0.653), lty = 2)
arrows(3,-0.2,y1 = 0.653, length = 0.1,
       code = 3)
text(x = 2.3, y = -0.75, 
     labels = bquote(paste(italic(alpha)['2'], "= 2", sep = "")),
     cex = .8)
text(x = 1.58, y = -0.33,
     labels = bquote(paste(italic(m)['0'], "= -0.2", sep = "")),
     cex = .8)
text(x = 3.4, y = 0.25,
     labels = bquote(paste(italic(alpha)['1'], "= 2.1", sep = "")),
     cex = .8)
plot(1~1, type = "n", xlim = c(1,4), ylim = c(-1,1), xlab = "",
     ylab = "", xaxt = "n", yaxt="n", bty = "n")

axis(1, at= c(1:4), labels = F, tck = -.05)
mtext(text = c(1:4),1, line = 0.3, at = c(1:4))

axis(2, at=seq(-1,1, 0.5), labels=F, tck=-.05)
axis(2, at=seq(-1,1, 0.25), labels=F, tck=-.03)
mtext(text = as.character(seq(-1,1,0.5)), 2, line = 0.4, at = seq(-1,1,0.5),
      las = 1)
mtext("B)",3, at = c(0.2), line = -0.1, cex = 1.2)
lines(p2(.1, 0.8, 2, 1:4, 2), lwd = 3)
abline(v = 2, lty = 2)
abline(h = 0.1, lty = 2)
lines(x = c(2,3), y = c(0.9, 0.9), lty = 2)
arrows(2.7,0.1,y1 = 0.9, length = 0.1,
       code = 3)
text(x = 2.3, y = -0.4, 
     labels = bquote(paste(italic(alpha)['2'], "= 2", sep = "")),
     cex = .8)
text(x = 3, y = -0.05,
     labels = bquote(paste(italic(m)['0'], "= 0.1", sep = "")),
     cex = .8)
text(x = 3.1, y = 0.5,
     labels = bquote(paste(italic(alpha)['1'], "= 0.8", sep = "")),
     cex = .8)
mtext(bquote(paste("logit(", italic(gamma), ")", sep = "")),
      2, at = 0, line = 1.7, cex = 1.2)
mtext("Time", 1, at = 2.5, line = 1.55, cex = 1.2)
dev.off()

uc <- read.csv('urban_covs.csv', header = TRUE)
up <- read.csv('urban_pc.csv', header = TRUE)

# get only covariates
uc <- uc[,-c(1:2)]

ucp <- prcomp(scale(uc))
loadings(ucp)
