
# some stuff for drop-k

# drop the last year of data and predict it

forecast_list <- list(y = as.matrix(y_array[which(species_names$x=="Coyote"),,]), nyear = ncol(z), 
                  nsite = nrow(z), 
                  spa = prior_for_occ$a, spb = prior_for_occ$b,
                  jmat = as.matrix(j_mat),
                  pi = 3.14159, C = cs[[1]], S = cs[[2]],
                  cov = covdat$pc1, P = 4)

forecast_list$y[,10:13] <- NA

ran_cast <- run.jags( model= "ranef_year_jags.R" , 
                      monitor= c("z" , "y") , 
                      data=forecast_list ,  
                      inits=inits_ranef , 
                      n.chains=detectCores()-1 ,
                      adapt=3000,
                      burnin=3000 , 
                      sample=ceiling(10000/7) ,
                      thin=5 ,
                      summarise=FALSE ,
                      plots=FALSE,
                      method = "parallel")

rr <- as.matrix(as.mcmc.list(ran_cast), chains = TRUE)[,-1]


rany <- rr[,grep("y", colnames(rr))]
ranz <- rr[,grep("z", colnames(rr))]
rm(rr)

ryf1 <- mase(rany, data_list, type = "naive", TRUE)
ryf2 <- mase(rany, data_list, type = "season", TRUE)
rzf1 <- mase(ranz, data_list, type = "naive", TRUE)
rzf2 <- mase(ranz, data_list, type = "season", TRUE)


mod_pulse_cast <- run.jags( model= "pulse_year_jags_trig.R" , 
                       monitor= c("y","Z") , 
                       data=forecast_list ,  
                       inits=inits_pulse , 
                       n.chains=detectCores()-1 ,
                       adapt=3000,
                       burnin=3000 , 
                       sample=ceiling(10000/7) ,
                       thin=5 ,
                       summarise=FALSE ,
                       plots=FALSE,
                       method = "parallel")

pp <- as.matrix(as.mcmc.list(mod_pulse_cast), chains = TRUE)[,-1]


puly <- pp[,grep("y", colnames(pp))]
pulz <- pp[,grep("z", colnames(pp))]
pyf1 <- mase(puly, data_list, type = "naive", TRUE)
pyf2 <- mase(puly, data_list, type = "season", TRUE)
rzf1 <- mase(ranz, data_list, type = "naive", TRUE)
rzf2 <- mase(ranz, data_list, type = "season", TRUE)


k <- 100
tok <- which(!is.na(data_list$y))
tok <- tok[-which(tok<101)]

####
sim_ans <- vector("list", length = 10)
for(i in 2:10){
  y <- data_list$y
  tg <- sample(1:length(tok), k, replace = FALSE)
  stored_val <- rep(0, length(tg))
  stored_val <- y[tok[tg]]
  old_val <- y[tok[tg]-100]
  y[tok[tg]] <- NA
  
  new_list <- list( y = as.matrix(y), nyear = ncol(z), nsite = nrow(z),
                    spa = prior_for_occ$a, spb = prior_for_occ$b,
                    jmat = as.matrix(j_mat),
                    cov = covdat$pc1, P = 4, pi = 3.141593)

  dropped <- run.jags( model= "ranef_year_jags.R" , 
                       monitor="y" , 
                       data=new_list ,  
                       inits=inits_ranef , 
                       n.chains=detectCores()-1 ,
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
  ets <- (sweep(preds, 2, stored_val))
  ets <- abs(ets)
  nf <- abs(stored_val - old_val)
  togo <- which(is.na(nf))
  et <- rowSums(ets[,-togo])/ (sum(nf, na.rm=TRUE) * 13/12)
  sim_ans[[i]] <- et
}

sim_ans_pulse <- vector("list", length = 10)
for(i in 2:10){
  y <- data_list$y
  tg <- sample(1:length(tok), k, replace = FALSE)
  stored_val <- rep(0, length(tg))
  stored_val <- y[tok[tg]]
  old_val <- y[tok[tg]-100]
  y[tok[tg]] <- NA
  
  new_list <- list( y = as.matrix(y), nyear = ncol(z), nsite = nrow(z),
                    spa = prior_for_occ$a, spb = prior_for_occ$b,
                    jmat = as.matrix(j_mat),
                    cov = covdat$pc1, P = 4, pi = 3.141593,
                    C = cs[[1]], S = cs[[2]])
  
  dropped <- run.jags( model= "pulse_year_jags_trig.R" , 
                       monitor="y" , 
                       data=new_list ,  
                       inits=inits_pulse , 
                       n.chains=detectCores()-1 ,
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
  ets <- (sweep(preds, 2, stored_val))
  ets <- abs(ets)
  nf <- abs(stored_val - old_val)
  togo <- which(is.na(nf))
  et <- rowSums(ets[,-togo])/ (sum(nf, na.rm=TRUE) * 13/12)
  sim_ans_pulse[[i]] <- et
}

mean(sapply(sim_ans_pulse, mean))
