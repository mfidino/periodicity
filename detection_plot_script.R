


get_things <- function(species_posterior = NULL){
  # get the site stuff
  site_ests <- species_posterior[,grep("epsilon", 
    colnames(species_posterior))]
  
  # get the mu
  f0 <- species_posterior[,grep("f0", colnames(species_posterior))]
  for(i in 1:(ncol(site_ests)-1)){
    site_ests[,i] <- site_ests[,i] + f0
  }
  
  
  # get quantiles
  sm <- apply(site_ests, 2, quantile, probs = c(0.025,0.5,0.975))
  
  
  
  # get mu and sd
  pmu <- species_posterior[,grep("lp|ls_sd", colnames(species_posterior))]
  
  plow <- median(qnorm(0.025, median(f0), site_ests[,96]) )
  phigh <- median(qnorm(0.975, median(f0), site_ests[,96]))
  
  return(list(sm = sm[,-96], plow = plow,
              phigh = phigh, 
              pmu = quantile(f0, probs = c(0.025,0.5,0.975))))
         
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
species_posterior <- fread("Fidino_2017_periodic_time_model_pulse_coyote_post#.txt", data.table = FALSE) 


# coyote
coy <- get_things(species_posterior)

# fox
species_posterior <- fread("Fidino_2017_periodic_time_model_pulse_redfox_post.txt", data.table = FALSE) 
fox <- get_things(species_posterior)
# skunk
species_posterior <- fread("Fidino_2017_stochastic_time_model_skunk_post.txt", data.table = FALSE)
sku <- get_things(species_posterior)
# raccoon
species_posterior <- fread("Fidino_2017_homogeneous_time_model_raccoon_post.txt", data.table = FALSE)
rac <- get_things(species_posterior)
# possum
species_posterior <- fread("Fidino_2017_periodic_time_model_boom_bust_opossum_post.txt", data.table = FALSE) 
pos <- get_things(species_posterior)









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


### initial occupancy plots

do_init <- function(x, a, b, addy = TRUE){
  prior <- density(rbeta(1e6, a, b))
  dens <- density(x)
  plot(1~1, type = "n", xlim = c(0,1), ylim = c(0,14), xlab = "",
    ylab = "", xaxt = "n", yaxt="n", bty = "n")
  axis(2, at= seq(0,14, 1), labels = F, tck = -.05)
  #axis(1, at= seq(0,0.5, 0.05), labels = F, tck = -.03)
  if(addy){
  mtext(text = as.character(seq(0,14,2)), 2, line = 0.8, at = seq(0,14,2),
    las = 1, cex = 0.75)
  }
  axis(1, at= seq(0,1, 0.25), labels = F, tck = -.05)
  axis(1, at= seq(0,1, 0.125), labels = F, tck = -.03)
  #axis(1, at= seq(0,0.5, 0.05), labels = F, tck = -.03)
  mtext(text = as.character(seq(0,1,.25)), 1, line = 0.8, at = seq(0,1,0.25),
    las = 1, cex = 0.75)
  lines(dens, lwd = 2)
  lines(prior, lwd = 2, lty = 2)
}


windows(6.5,6.5)


coypsi <- fread("Fidino_2017_periodic_time_model_pulse_coyote_post.txt", data.table = FALSE, select = "psinit")$psinit 



# fox
foxpsi <- fread("Fidino_2017_periodic_time_model_pulse_redfox_post.txt", 
  data.table = FALSE, select = "psinit")$psinit


# skunk
skupsi <- fread("Fidino_2017_stochastic_time_model_skunk_post.txt", data.table = FALSE, select = "psinit")$psinit


# raccoon
racpsi <- fread("Fidino_2017_homogeneous_time_model_raccoon_post.txt", data.table = FALSE, select = "psinit")$psinit

nms <- c("A) Coyote", "B) Red fox", "C) Striped skunk",
  "D) Raccoon","E) Opossum")
nmsa <- c(1.2, 1.25, 1.85, 1.3, 1.35)
# possum
opopsi <- fread("Fidino_2017_periodic_time_model_boom_bust_opossum_post.txt", data.table = FALSE, select = "psinit")$psinit 

tiff("initial_occupany2.tiff", height = 6.5, width = 6.5, units = "in",
  res = 600, compression = "lzw")
m <- matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2, byrow = TRUE)
layout(m)
par( mar = c(2,2.5, 1.25,0),
  oma = c(2,2.5,1.25,0) + 0.1)
do_init(coypsi, 24.35, 31.08)
mtext(nms[1], 3, at = 0.12)
legend("topright", c("Posterior", "Prior"), lty = c(1,2), lwd = c(2,2),
  col = c("black", "black"), box.col = "white", box.lwd = 0, bg = "white")
mtext("Density", 2, line = 2.4, at = 7)
do_init(foxpsi, 2.13, 26.69, FALSE)
mtext(nms[2], 3, at = 0.12)
legend("topright", c("Posterior", "Prior"), lty = c(1,2), lwd = c(2,2),
  col = c("black", "black"), box.col = "white", box.lwd = 0, bg = "white")
do_init(skupsi, 5.70, 30.78)
legend("topright", c("Posterior", "Prior"), lty = c(1,2), lwd = c(2,2),
  col = c("black", "black"), box.col = "white", box.lwd = 0, bg = "white")
mtext("Density", 2, line = 2.4, at = 7)
mtext(nms[3], 3, at = 0.205)
do_init(racpsi, 24.35, 31.08, FALSE)
legend("topright", c("Posterior", "Prior"), lty = c(1,2), lwd = c(2,2),
  col = c("black", "black"), box.col = "white", box.lwd = 0, bg = "white")
mtext(nms[4], 3, at = 0.13)
do_init(opopsi, 36.03, 29.78)
legend("topright", c("Posterior", "Prior"), lty = c(1,2), lwd = c(2,2),
  col = c("black", "black"), box.col = "white", box.lwd = 0, bg = "white")
mtext("Density", 2, line = 2.4, at = 7)
mtext(nms[5], 3, at = .145)
mtext("Probability of Initial Occupancy", side = 1, line = 2.5,at = 0.5)
dev.off()