get_things <- function(sp = NULL, tg = tg, dotg = TRUE){
  # get the site stuff
  sm <- sp[,grep("ls\\.", colnames(sp))]
  
  # get the mu
  lp <- sp[,grep("lp", colnames(sp))]
  for(i in 1:ncol(sm)){
    sm[,i] <- sm[,i] + lp
  }
  
  # remove the blank sites
  if(dotg) sm <- sm[,-tg]
  
  zs <- sp[,grep("z\\.", colnames(sp))]
  if(dotg){
  zs <- matrix(apply(zs, 2, median), nrow = 100, ncol = 9)
  zs <- zs[-tg,]
  }else{
    zs <- matrix(apply(zs, 2, median), nrow = 95, ncol = 9)

  }
  pres <- rep(0, 8)
  for(i in 2:ncol(zs)){
    boop <- zs[,i] + zs[,i-1]
    pres[i-1] <- sum(boop==2)/ (sum(zs[,i-1]==1 & zs[,i]==0) + sum(boop==2))
  }
  zs <- apply(zs, 2, mean)
  
  
  # get quantiles
  sm <- apply(sm, 2, quantile, probs = c(0.025,0.5,0.975))
  
  
  
  # get mu and sd
  pmu <- sp[,grep("lp|ls_sd", colnames(sp))]
  
  plow <- median(qnorm(0.025, pmu[,1], pmu[,2]) )
  phigh <- median(qnorm(0.975, pmu[,1], pmu[,2]))
  
  return(list(sm = sm, plow = plow,
              phigh = phigh, 
              pmu = quantile(pmu[,1], probs = c(0.025,0.5,0.975)),
              zs = zs, pres = pres))
         
}
prpl <- function(op = NULL){
  plow <- op$plow
  phigh <- op$phigh
  sm <- op$sm
  pmu <- op$pmu
  plot(1~1, type = "n", xlim = c(0,95), ylim = c(-9,3), xlab = "",
       ylab = "", xaxt = "n", yaxt="n", bty = "n")
  
  #
  
  axis(2, at= seq(-9,3, 1), labels = F, tck = -.05)
  #axis(1, at= seq(0,0.5, 0.05), labels = F, tck = -.03)
  mtext(text = as.character(seq(-9,3,3)), 2, line = 0.6, at = seq(-9,3,3),
        las = 1, cex = 0.75)
  #
  rect( 0, plow, 96, phigh, col = "gray80")
  rect(0, pmu[1], 96, pmu[3], col = "gray50", border = FALSE)
  abline(h = pmu[2], lwd = 1)
  for(i in 1:ncol(sm)){
    lines(c(i,i), sm[-2,sor[i]])
  }
  points(c(1:95),sm[2,sor], pch = 19, cex = 0.8)
  
}
# coyote
sp <- read.table("./model_outputs/coyote4_only_pulse_year_jags_trig.txt", header = TRUE,
                 sep = "\t")
z <- df_2_array(read.table("z_matrix_sp10_sp13.txt", header = TRUE, sep = "\t"))[1,,-c(1:4)]

tg <- which(rowSums(is.na(z))==9)

# coyote
coy <- get_things(sp, tg)
# fox
sp <- read.table("./model_outputs/redfox3_only_pulse_year_jags_trig.txt", header = TRUE,
                 sep = "\t")
fox <- get_things(sp, tg)
# skunk
sp <- read.table("./model_outputs/skunk3_only_pulse_year_jags_trig.txt", header = TRUE,
                 sep = "\t")
sku <- get_things(sp, tg)
# raccoon
sp <- read.table("./model_outputs/raccoon4_ranef_year_jags.txt", header = TRUE,
                 sep = "\t")
rac <- get_things(sp, tg)
# possum
sp <- read.table("./model_outputs/opossum4_boom_bust_jags_trig.txt", header = TRUE,
                 sep = "\t")
pos <- get_things(sp, tg, FALSE)


mns <- data.frame(rac$zs, coy$zs, pos$zs, sku$zs, fox$zs)
mns <- data.frame(rac$pres, coy$pres, pos$pres, sku$pres, fox$pres)
windows(5,3)
tiff("sp_occ.tiff", height = 3, width = 5, units = "in",
     res = 600, compression = 'lzw')
par( mar = c(1,2, .6,0),
     oma = c(1,2,.6,0) + 0.1)

plot(1~1, type = "n", xlim = c(.5,5.5), ylim = c(0,1), xlab = "",
     ylab = "", xaxt = "n", yaxt="n", bty = "n")

#

axis(1, at= seq(1,5, 1), labels = F, tck = -.05)
#axis(1, at= seq(0,0.5, 0.05), labels = F, tck = -.03)
mtext(text = c("Raccoon", "Coyote", "Opossum",
               "Skunk", "Red fox"), 1, line = 0.6, at = seq(1,5,1),
      las = 1, cex = 1)
axis(2, at = seq(0, 1, 0.1), labels = F, tck = -0.05)
axis(2, at = seq(0, 1, 0.05), labels = F, tck = -0.03)
mtext(text = seq(0,1,0.2), 2, line = 0.6, at = seq(0,1,0.2),
      las = 2)
mtext("Proportion of sites occupied", 2, line =2)

boxplot(mns, add = TRUE, outline = FALSE, names = "")
dev.off()
# order them
sor <- order(sm[2,])

lims <- c(0.75,0.5,0.5, 1, 0.75)
nms <- c("A) Coyote", "B) Red fox", "C) Striped skunk",
         "D) Raccoon","E) Opossum")
nmsa <- c(1.2, 1.25, 1.85, 1.3, 1.35)

  tiff("detection_plot.tiff", height = 8, width = 6.25, units = "in",
       res = 600, compression = "lzw")
  m <- matrix(c(1:5), ncol = 1, nrow = 5 )
  layout(m)

  par( mar = c(1,2, 1.25,0),
       oma = c(1,2,1.25,0) + 0.1)
  prpl(coy)
  mtext(nms[1], 3, at = 2.5)
  #abline(v = -2.7)
  prpl(fox)
  #abline(v = -2.7)
  mtext(nms[2], 3, at = 3)
  prpl(sku)
  #abline(v = -2.7)
  mtext(nms[3], 3, at = 7.5)
  mtext("logit ( probability of detection )", 2, line = 2.4, at = -3)
  prpl(rac)
  #abline(v = -2.7)
  mtext(nms[4], 3, at = 4)
  prpl(pos)
  #abline(v = -2.7)
  mtext(nms[5], 3, at = 4.5)
  mtext("Site estimates", side = 1, line = 0)
  dev.off()
 


  abline(h = median(pmu[,1]), cex = 2)
  abline(h = quantile(pmu[,1], probs = 0.025), cex = 2, lty = 2)
  abline(h = quantile(pmu[,1], probs = 0.975), cex = 2, lty = 2)
  
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




plot(1~)
s2_coy <- s2
dp <- sp[,grep("dp", colnames(sp))]
rm(sp)

do_flp <- function(x, y){
  a <- x[,2,y]
  b <- x[,1,y]
  return(a/b)
  
}
t(apply(s2_coy, 2, quantile, probs = c(0.025,0.5,0.975)))
t(apply(s2_fox, 2, quantile, probs = c(0.025,0.5,0.975)))
t(apply(s2_skunk, 2, quantile, probs = c(0.025,0.5,0.975)))
t(apply(s2_possum, 2, quantile, probs = c(0.025,0.5,0.975)))
t(apply(s2_raccoon, 2, quantile, probs = c(0.025,0.5,0.975)))
ex(quantile(s2_coy[,3], pro))


ob_per <- rep(0, 8)
for(i in 2:9){
  #a2[,i-1] <- z[,i] - z[,i-1]
  boop <- z[,i] + z[,i-1]
  ob_per[i-1] <- sum(boop==2, na.rm = TRUE)/
    (100 - sum(is.na(boop)))
}

kry <- matrix(0, ncol = length(kk), nrow = 5)
kry[1,] <- kk


# get the proportion of sites colonized by sp
# get the proportion of sites colonized by sp
ans <- matrix(0, ncol = 8, nrow = 95)
av <- matrix(0, nrow = nrow(spz), ncol = 8)
for(j in 1:nrow(spz)){
  nz <- matrix(as.numeric(spz[j,]), ncol = 9, nrow = 100)
  nz <- nz[-tg,]
  for(i in 2:9){
    ans[,i-1] <- nz[,i] - nz[,i-1]
    #a2[,i-1] <- z[-tg,i] - z[-tg,i-1]
    #nna <- which(is.na(a2[,i-1]))
    boop <- nz[,i] + nz[,i-1]
    #b2 <- z[,i] + z[,i-1]
    #nna2 <- which(is.na(b2))
    av[j,i-1] <- sum(ans[,i-1]==1, na.rm = TRUE)/
      (95 - sum(boop==2))
  }
}

mye <- apply(av, 2, quantile, probs= c(0.05,0.5,0.95))

mry <- array(0, dim = c(nrow(mye), ncol(mye), 5))
mry[,,1] <- mye

# make logit_preds
thet <- s2[,grep("theta", colnames(s2))]
gmu <- s2[,grep("g_mu", colnames(s2))]
#hist(thet)
ack <- matrix(0, , ncol = 8, nrow = length(thet))
# for the pulse species
for(i in 1:nrow(ack)){
  ack[i,] <- plsg(thet[i], 2, 4, 1:8, gmu[i] )
}
ex <- function(x) 1 / (1 + exp(-x))
am <- ex(apply(ack, 2, quantile, probs = c(0.025,0.5,0.975)))
ary <- array(0, dim=c(nrow(am), ncol(am), 5))
ary[,,1] <- am
rm(spz)
rm(ack)
rm(av)
rm(s2)
################# red fox

sp <- read.table("./model_outputs/redfox3_only_pulse_year_jags_trig.txt", header = TRUE,
                 sep = "\t")
z <- df_2_array(read.table("z_matrix_sp10_sp13.txt", header = TRUE, sep = "\t"))[5,,-c(1:4)]

# get the z from sp
spz <- sp[,grep("z", colnames(sp))]

s2 <- sp[,-grep("pred", colnames(sp))]
s2 <- s2[,-grep("p\\.", colnames(s2))]
s2 <- s2[,-grep("z\\.", colnames(s2))]
s2 <- s2[,-grep("ls\\.", colnames(s2))]
s2_fox <- s2
dp <- sp[,grep("dp", colnames(sp))]
rm(sp)

#t(apply(s2, 2, quantile, probs = c(0.025,0.5,0.975)))


a2 <-  matrix(0, ncol = 8, nrow = 100)
kk <- rep(0, 8)
for(i in 2:9){
  a2[,i-1] <- z[,i] - z[,i-1]
  boop <- z[,i] + z[,i-1]
  kk[i-1] <- sum(a2[,i-1]==1, na.rm = TRUE)/
    (100 - sum(is.na(boop))- sum(boop==2, na.rm = TRUE))
}

kry[2,] <- kk


# get the proportion of sites colonized by sp
# get the proportion of sites colonized by sp
ans <- matrix(0, ncol = 8, nrow = 95)
av <- matrix(0, nrow = nrow(spz), ncol = 8)
for(j in 1:nrow(spz)){
  nz <- matrix(as.numeric(spz[j,]), ncol = 9, nrow = 100)
  nz <- nz[-tg,]
  for(i in 2:9){
    ans[,i-1] <- nz[,i] - nz[,i-1]
    # a2[,i-1] <- z[-tg,i] - z[-tg,i-1]
    #nna <- which(is.na(a2[,i-1]))
    boop <- nz[,i] + nz[,i-1]
    #b2 <- z[,i] + z[,i-1]
    #nna2 <- which(is.na(b2))
    av[j,i-1] <- sum(ans[,i-1]==1, na.rm = TRUE)/
      (95 - sum(boop==2))
  }
}

mye <- apply(av, 2, quantile, probs= c(0.05,0.5,0.95))

mry[,,2] <- mye

# make logit_preds
thet <- s2[,grep("theta", colnames(s2))]
gmu <- s2[,grep("g_mu", colnames(s2))]
#hist(thet)
ack <- matrix(0, , ncol = 8, nrow = length(thet))
# for the pulse species
for(i in 1:nrow(ack)){
  ack[i,] <- plsg(thet[i], 2, 4, 1:8, gmu[i] )
}
ex <- function(x) 1 / (1 + exp(-x))
am <- ex(apply(ack, 2, quantile, probs = c(0.025,0.5,0.975)))
ary[,,2] <- am

rm(spz)
rm(ack)
rm(av)
rm(s2)

################# skunk

sp <- read.table("./model_outputs/skunk3_only_pulse_year_jags_trig.txt", header = TRUE,
                 sep = "\t")
z <- df_2_array(read.table("z_matrix_sp10_sp13.txt", header = TRUE, sep = "\t"))[6,,-c(1:4)]

# get the z from sp
spz <- sp[,grep("z", colnames(sp))]

s2 <- sp[,-grep("pred", colnames(sp))]
s2 <- s2[,-grep("p\\.", colnames(s2))]
s2 <- s2[,-grep("z\\.", colnames(s2))]
s2 <- s2[,-grep("ls\\.", colnames(s2))]
s2_skunk <- s2
dp <- sp[,grep("dp", colnames(sp))]
rm(sp)

t(apply(s2, 2, quantile, probs = c(0.025,0.5,0.975)))


a2 <-  matrix(0, ncol = 8, nrow = 100)
kk <- rep(0, 8)
for(i in 2:9){
  a2[,i-1] <- z[,i] - z[,i-1]
  boop <- z[,i] + z[,i-1]
  kk[i-1] <- sum(a2[,i-1]==1, na.rm = TRUE)/
    (100 - sum(is.na(boop))- sum(boop==2, na.rm = TRUE))
}

kry[3,] <- kk


# get the proportion of sites colonized by sp
# get the proportion of sites colonized by sp
ans <- matrix(0, ncol = 8, nrow = 95)
av <- matrix(0, nrow = nrow(spz), ncol = 8)
for(j in 1:nrow(spz)){
  nz <- matrix(as.numeric(spz[j,]), ncol = 9, nrow = 100)
  nz <- nz[-tg,]
  for(i in 2:9){
    ans[,i-1] <- nz[,i] - nz[,i-1]
    #a2[,i-1] <- z[-tg,i] - z[-tg,i-1]
    #nna <- which(is.na(a2[,i-1]))
    boop <- nz[,i] + nz[,i-1]
    #b2 <- z[,i] + z[,i-1]
    #nna2 <- which(is.na(b2))
    av[j,i-1] <- sum(ans[,i-1]==1, na.rm = TRUE)/
      (95 - sum(boop==2))
  }
}

mye <- apply(av, 2, quantile, probs= c(0.05,0.5,0.95))

mry[,,3] <- mye

# make logit_preds
thet <- s2[,grep("theta", colnames(s2))]
gmu <- s2[,grep("g_mu", colnames(s2))]
#hist(thet)
ack <- matrix(0, , ncol = 8, nrow = length(thet))
# for the pulse species
for(i in 1:nrow(ack)){
  ack[i,] <- plsg(thet[i], 2, 4, 1:8, gmu[i] )
}
ex <- function(x) 1 / (1 + exp(-x))
am <- ex(apply(ack, 2, quantile, probs = c(0.025,0.5,0.975)))
ary[,,3] <- am

rm(spz)
rm(ack)
rm(av)
rm(s2)

################# raccoon

sp <- read.table("./model_outputs/raccoon4_ranef_year_jags.txt", header = TRUE,
                 sep = "\t")
z <- df_2_array(read.table("z_matrix_sp10_sp13.txt", header = TRUE, sep = "\t"))[4,,-c(1:4)]

# get the z from sp
spz <- sp[,grep("z", colnames(sp))]

s2 <- sp[,-grep("pred", colnames(sp))]
s2 <- s2[,-grep("p\\.", colnames(s2))]
s2 <- s2[,-grep("z\\.", colnames(s2))]
s2 <- s2[,-grep("ls\\.", colnames(s2))]
s2_raccoon <- s2
dp <- sp[,grep("a2_gam", colnames(sp))]
rm(sp)

t(apply(s2, 2, quantile, probs = c(0.025,0.5,0.975)))


a2 <-  matrix(0, ncol = 8, nrow = 100)
kk <- rep(0, 8)
for(i in 2:9){
  a2[,i-1] <- z[,i] - z[,i-1]
  boop <- z[,i] + z[,i-1]
  kk[i-1] <- sum(a2[,i-1]==1, na.rm = TRUE)/
    (100 - sum(is.na(boop))- sum(boop==2, na.rm = TRUE))
}

kry[4,] <- kk


# get the proportion of sites colonized by sp
# get the proportion of sites colonized by sp
ans <- matrix(0, ncol = 8, nrow = 95)
av <- matrix(0, nrow = nrow(spz), ncol = 8)
for(j in 1:nrow(spz)){
  nz <- matrix(as.numeric(spz[j,]), ncol = 9, nrow = 100)
  nz <- nz[-tg,]
  for(i in 2:9){
    ans[,i-1] <- nz[,i] - nz[,i-1]
    #a2[,i-1] <- z[-tg,i] - z[-tg,i-1]
    #nna <- which(is.na(a2[,i-1]))
    boop <- nz[,i] + nz[,i-1]
    #b2 <- z[,i] + z[,i-1]
    #nna2 <- which(is.na(b2))
    av[j,i-1] <- sum(ans[,i-1]==1, na.rm = TRUE)/
      (95 - sum(boop==2))
  }
}

mye <- apply(av, 2, quantile, probs= c(0.05,0.5,0.95))

mry[,,4] <- mye

# make logit_preds
#thet <- s2[,grep("a1_gam", colnames(s2))]
gmu <- s2[,grep("g_mu", colnames(s2))]
gy <- s2[,grep("gy", colnames(s2))]
# remove sd
gy <- gy[,-9]
gy <- apply(gy, 2, quantile, probs = c(0.025,0.5,0.975))
#hist(thet)
#ack <- matrix(0, , ncol = 8, nrow = length(thet))
# for the pulse species
#for(j in 1:nrow(ack)){
# ack[j,] <- p2(0, thet[j], 0, 1:8, 2 )
#}
# add gy to ack
#ack <- ack + gy
ex <- function(x) 1 / (1 + exp(-x))
am <- ex(gy)


ary[,,4] <- am



rm(spz)
rm(ack)
rm(av)
rm(s2)

################# opossum

#tg <- which(rowSums(is.na(z))==9)
sp <- read.table("./model_outputs/opossum4_boom_bust_jags_trig.txt", header = TRUE,
                 sep = "\t")
z <- df_2_array(read.table("z_matrix_sp10_sp13.txt", header = TRUE, sep = "\t"))[3,,-c(1:4)]

# get the z from sp
spz <- sp[,grep("z", colnames(sp))]

s2 <- sp[,-grep("pred", colnames(sp))]
s2 <- s2[,-grep("p\\.", colnames(s2))]
s2 <- s2[,-grep("z\\.", colnames(s2))]
s2 <- s2[,-grep("ls\\.", colnames(s2))]
s2_possum <- s2
dp <- sp[,grep("a2_gam", colnames(sp))]
rm(sp)

t(apply(s2, 2, quantile, probs = c(0.025,0.5,0.975)))


# with tg
a2 <-  matrix(0, ncol = 8, nrow = 95)
kk <- rep(0, 8)
for(i in 2:9){
  a2[,i-1] <- z[-tg,i] - z[-tg,i-1]
  boop <- z[-tg,i] + z[-tg,i-1]
  kk[i-1] <- sum(a2[-tg,i-1]==1, na.rm = TRUE)/
    (95 - sum(is.na(boop))- sum(boop==2, na.rm = TRUE))
}


kry[5,] <- kk
write.table(kry, "obs_z_col.txt", sep = "\t", row.names = FALSE)

# get the proportion of sites colonized by sp
ans <- matrix(0, ncol = 8, nrow = 95)
av <- matrix(0, nrow = nrow(spz), ncol = 8)
for(j in 1:nrow(spz)){
  nz <- matrix(as.numeric(spz[j,]), ncol = 9, nrow = 95)
  #nz <- nz[-tg,]
  for(i in 2:9){
    ans[,i-1] <- nz[,i] - nz[,i-1]
    #a2[,i-1] <- z[-tg,i] - z[-tg,i-1]
    #nna <- which(is.na(a2[,i-1]))
    boop <- nz[,i] + nz[,i-1]
    #b2 <- z[,i] + z[,i-1]
    #nna2 <- which(is.na(b2))
    av[j,i-1] <- sum(ans[,i-1]==1, na.rm = TRUE)/
      (95 - sum(boop==2))
  }
}

mye <- apply(av, 2, quantile, probs= c(0.05,0.5,0.95))

mry[,,5] <- mye
saveRDS(mry, "z_preds_col.RDS")

# make logit_preds
thet <- s2[,grep("a1_gam", colnames(s2))]
gmu <- s2[,grep("g_mu", colnames(s2))]
gmu <- quantile(gmu, probs = c(0.025,0.5,0.975))
gmu <- matrix(gmu, ncol = 8, nrow = 3)

gy <- s2[,grep("gy", colnames(s2))]
# remove sd
gy <- gy[,-9]
#gy <- apply(gy, 2, quantile, probs = c(0.025,0.5,0.975))
#gy <- s2[,grep("gy", colnames(s2))]
# remove sd
#gy <- gy[,-9]
#hist(thet)
ack <- matrix(0, , ncol = 8, nrow = length(thet))
# for the pulse species
for(j in 1:nrow(ack)){
  ack[j,] <- p2(0, thet[j], 0, 1:8, 2 )
}
# add gy to ack
ack <- ack + gy
ex <- function(x) 1 / (1 + exp(-x))
am <- ex(apply(ack, 2, quantile, probs = c(0.025,0.5,0.975)))
ary[,,5] <- am
saveRDS(ary, "col_preds.RDS")


rm(spz)
rm(ack)
rm(av)
rm(s2)

plot(am[2,], type = 'l', ylim = c(0,1))
lines(am[1,], type = 'l', lty = 2)
lines(am[3,], type = 'l', lty = 2)
points(kk, col = 'red')
points(mye[2,], pch = 21, col = 'black')
# species plots, modify for 5 species
# should probably rescale the y-axes

