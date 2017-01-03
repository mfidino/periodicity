
pls <- function(theta = NULL, delta = NULL, t = NULL, a0 = NULL, av = NULL){

  n <- 1:4  
  a1 <- (theta/(1*pi))*sin((pi*1)/4)*cos((pi*1 * (t - delta))/2)
  a2 <- (theta/(2*pi))*sin((pi*2)/4)*cos((pi*2 * (t - delta))/2)
  a3 <- (theta/(3*pi))*sin((pi*3)/4)*cos((pi*3 * (t - delta))/2)
  a4 <- (theta/(4*pi))*sin((pi*4)/4)*cos((pi*4 * (t - delta))/2)
  ans <- a0 + a1 +a2 + a3 + a4 + av
  return(ans)
}

pp <- matrix(0, nrow = length(theta), ncol = length(1:12))
for(i in 1:length(theta)){
  pp[i,] <- pls(theta[i], 2, seq(1,12,1), 0, gys[1:12])
}

pq <- apply(pp, 2, quantile, probs = c(0.025,0.5,0.975))

pqp <- 1 / (1 + exp(-pq))
windows(width = 3, height = 3)
tiff("tiny.tiff", height = 2.5, width = 3, units = "in", res = 600,
     compression = "lzw")
par( mar = c(1,1.5, .35,0.2),
     oma = c(1,1.5,0.6,0.5) + 0.1)
plot(pqp[2,], type = "l", ylim = c(0,0.5), bty = "n", xlab = "Season",
     ylab = "Probability of colonization", xaxt = "n", yaxt = "n", lwd = 2)

axis(1, at= c(1:12), labels = F, tck = -.02)
a <- c("SP", "SU", "FA","WI")
my_text <- paste0(paste0(a, c(rep(10,3), rep(11,4), rep(12,4), rep(13,1))),"\n-\n",
       c("SU10", "FA10", "WI11", "SP11", "SU11", "FA11", "WI12", "SP12", "SU12",
         "FA12", "WI13", "SP13"))
  mtext(text = 1:12, 1, line = 0.1, at = 1:12, cex = 0.8)
  mtext("Season",1, at = 6.5, line = 0.9, cex = 1)

axis(2, at=seq(0, 0.5, 0.1), labels=F, tck=-.02)
mtext(text = seq(0,0.5, 0.1),2, line = 0.4, at = seq(0,0.5,0.1), las = 2, cex = 0.8)
mtext(text = "Probability of colonization", 2, line = 1.8, at = 0.25)
mtext(text = "A) Coyote", 3, at = 1.4, line = 0.2)


lines(pqp[1,], lty = 2)

lines(pqp[3,], lty = 2)

for(i in 1:12){
  lines( c(i,i), c(mye[1,i], mye[3,i]))
}
points(mye[2,], pch = 21, bg = "azure3")

dev.off()

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
p2ranef <- function(gs, a1, a2, ti, p){
  p2m <- rep(0, length(ti))

  return(p2m)
} 


# doing opossum

opo <- matrix(0, nrow = length(a1p), ncol = length(1:12))
gys <- apply(gys, 2, median)
pys <- apply(pys, 2, median)

for(i in 1:length(a1p)){
  opo[i,] <- p2(gys, a1g[i],0, seq(1,12,1), 2)
}


opoq <- apply(opo, 2, quantile, probs = c(0.025,0.5,0.975))

opoq <- 1 / (1 + exp(-opoq))
opoqe <- 1 - opoq

for(i in 1:length(a1p)){
  opo[i,] <- p2(gys, a1g[i],0, seq(1,12,1), 2)
}


opoq <- apply(opo, 2, quantile, probs = c(0.025,0.5,0.975))

opoq <- 1 / (1 + exp(-opoq))



windows(width = 3, height = 3)
tiff("tiny_pos.tiff", height = 2.5, width = 3, units = "in", res = 600,
     compression = "lzw")
par( mar = c(1,1.5, .35,0.2),
     oma = c(1,1.5,0.6,0.5) + 0.1)
plot(opoq[2,], type = "l", ylim = c(0,0.5), bty = "n", xlab = "Season",
     ylab = "Probability", xaxt = "n", yaxt = "n", lwd = 2)

axis(1, at= c(1:12), labels = F, tck = -.02)

mtext(text = 1:12, 1, line = 0.1, at = 1:12, cex = 0.8)
mtext("Season",1, at = 6.5, line = 0.9, cex = 1)

axis(2, at=seq(0, 0.5, 0.1), labels=F, tck=-.02)
mtext(text = seq(0,0.5, 0.1),2, line = 0.4, at = seq(0,0.5,0.1), las = 2, cex = 0.8)
mtext(text = "Probability of colonizatoin", 2, line = 1.8, at = 0.25)
mtext(text = "B) Opossum", 3, at = 1.8, line = 0.2)


lines(opoq[1,], lty = 2)

lines(opoq[3,], lty = 2)

#lines(opoqe[2,], col = "grey40", lwd = 2)
#lines(opoqe[1,], col = "grey40", lty = 2)
#lines(opoqe[3,], col = "grey40", lty = 2)

for(i in 1:12){
  lines( c(i,i), c(mye[1,i], mye[3,i]))
}
points(mye[2,], pch = 21, bg = "azure3")

dev.off()


plot(p2(.1, 0.8, 2, 1:12, 2), type = "l")

test <- plsg(2.1, 0, 2, 1:8, -.2)
plot(test, type = "l")

windows(height = 4, width = 3)
tiff("fig_1_periodicity2.tiff", height = 4, width = 3, units = "in",
     compression = "lzw", res = 600)
par(mfrow = c(2,1))
par( mar = c(1.4,2, .75,0.1),
     oma = c(1.5,1.5,1,0.1) + 0.1)
plot(1~1, type = "n", xlim = c(1,4), ylim = c(-1,1), xlab = "",
     ylab = "", xaxt = "n", yaxt="n", bty = "n")

axis(1, at= c(1:4), labels = F, tck = -.05)
mtext(text = c(1:4),1, line = 0.3, at = c(1:4))

axis(2, at=seq(-1,1, 0.5), labels=F, tck=-.05)
axis(2, at=seq(-1,1, 0.25), labels=F, tck=-.03)
mtext(text = as.character(seq(-1,1,0.5)), 2, line = 0.4, at = seq(-1,1,0.5),
      las = 1)

mtext(bquote(paste("logit(", italic(gamma), ")", sep = "")),
      2, at = 0, line = 1.7, cex = 1.2)
mtext("A)",3, at = c(0.2), line = -0.1, cex = 1.2)
lines(plsg(2.1, 2, 4, 1:4, -.2), lwd = 3)
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


pulser <- function(tau, p, n, ti){
  ans <- matrix(0, ncol = length(ti), nrow = length(n))
  nn <- n
  for(i in 1:length(nn)){
    ans[i,] <- tau/p + ((2/(n[i]*pi))*sin((pi*n[i]*tau)/p)*cos((2*pi*n[i]*(ti-(tau/2))/p)))
     #ans[i,] <- tau/p + (p/(pi*(2*i + tau))) * sin(2*pi * (2*i+tau) * ti)
    #ans[i,] <- (4/pi) * sin((2*pi*ti * i)/p)
    }
  return(colSums(ans))
}

plot(pulser(1, 4, seq(1,45,2), 1:8) , type = "l")
tst <- matrix(0, nrow = nrow(mm), ncol = length(1:12))
for(i in 1:nrow(mm)){
  tst[i,] <- pls(mm[i,45], 2, 1:12,median(mm[,48]),0)
}

tp <- apply(tst, 2, quantile, probs = c(0.025, 0.5, 0.975))
tp <- 1 / (1 + exp(-tp))
plot(tp[2,], ylim = c(0, 0.5), type = "l", col = "black", xlab = "season",
     ylab = "probability of colonization")
lines(tp[1,], col = "grey")
lines(tp[3,], col = "grey")

points(av, pch = 15, col = "red")
for(i in 2:nrow(tst)){
  lines(tst[i,], col = "grey")
}
test <- lapply(as.matrix(mm[,45]), pls, delta = 2, t = 1:12 )
plot(-4.1 +mt[18:29], type = "l" )
ans <- pls(0.42, 2, 1:12, -4.1, mt[18:29])
pls2 <- function(theta = NULL, delta = NULL, ti = NULL){
a <- 1:3
aa <- matrix(0, nrow = length(ti), ncol = length(a))
for(i in a){
  aa[,i] <- theta * cos((pi*i*delta)/2) * (cos((pi*i*ti)/2) * sin((pi*i)/4))/(pi*i) + 
            theta * sin((pi*i*delta)/2) * (sin((pi*i*ti)/2) * sin((pi*i)/4))/(pi*i)

}

ans <- rowSums(aa)
return(ans)

}


make_c_s <- function( ti = NULL, n = NULL){
  C <- matrix(0, nrow = length(ti), ncol = n)
  S <- matrix(0, nrow = length(ti), ncol = n)
  for(i in 1:n){
    C[,i] <- (cos((pi*i*ti)/2) * sin((pi*i)/4))/(pi*i)
    S[,i] <- (sin((pi*i*ti)/2) * sin((pi*i)/4))/(pi*i)
  }
  return(list(C=C, S=S))
}

pls3 <- function(theta = NULL, delta = NULL, ti = NULL, cs = NULL){
  a <- 1:3
  aa <- matrix(0, nrow = length(ti), ncol = length(a))
  for(i in a){
    aa[,i] <- theta * cos((pi*i*delta)/2) * cs[[1]][,i] + 
      theta * sin((pi*i*delta)/2) * cs[[2]][,i]
    
  }
  ans <- rowSums(aa)
  return(ans)
}

longshot <- pls3(0.49, 2, 1:12, cs)
  
ans <- pls2(0.49, 2, 1:13)
