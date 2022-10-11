### simulate 100 times ###
cpara.WLS <- cpara.Hom <- cpara.Hom.OGP <- cpara.Het <- cpara.Het.OGP <- 
  ll.Hom <- ll.Hom.OGP <- ll.Het <- ll.Het.OGP <- 
  tm.Hom <- tm.Hom.OGP <- tm.Het <- tm.Het.OGP <- 
  rmse.WLS <- rmse.Hom <- rmse.Hom.OGP <- rmse.Het <- rmse.Het.OGP <- 
  score.Hom <- score.Hom.OGP <- score.Het <- score.Het.OGP <- coverage <- rep(0, 100)

for(ii in 1:100){
  print(ii)
  set.seed(ii)
  
  # simulate X and Z
  X <- matrix(rep(X0, n.rep), ncol = 1)
  Z <- rep(0, sum(n.rep))
  for(i in 1:length(X0)) {
    Z[(ifelse(i==1,0,sum(n.rep[1:(i-1)]))+1):sum(n.rep[1:i])] <- pmean[i] + rnorm(n.rep[i], 0, sd = sqrt(var.y[i]))
  }
  
  ## WLS estimate ------------------------------------------------------------ 
  Z0 <- hetGP::find_reps(X, Z)$Z0
  Sigma_inv <- diag(1/sapply(hetGP::find_reps(X, Z)$Zlist, var))
  min.index <- which.min(apply(matrix(seq(cpara_min, cpara_max, length.out = 11),ncol=1), 1, function(x) optim(x, fn = function(g) t(Z0 - f.sim(X0, g)) %*% Sigma_inv %*% (Z0 - f.sim(X0, g)), 
                                                                                                               lower = cpara_min, upper = cpara_max, method = "L-BFGS-B")$value))
  
  cpara.WLS[ii] <- optim(seq(cpara_min, cpara_max, length.out = 11)[min.index], 
                         fn = function(g) t(Z0 - f.sim(X0, g)) %*% Sigma_inv %*% (Z0 - f.sim(X0, g)), 
                         lower = cpara_min, upper = cpara_max, method = "L-BFGS-B")$par
  xgrid <- matrix(seq(min(X0), max(X0), length.out = 101), ncol = 1) 
  rmse.WLS[ii] <- sqrt(mean((p.fun(xgrid) - f.sim(xgrid, cpara.WLS[ii]))^2))
  
  ## Hom without orthogonality ------------------------------------------------------------
  model <- vector("list", 3)
  jj <- 0
  for(cpara.init in cpara_init.vt){
    jj <- jj + 1
    model[[jj]] <- mleHomCalibrate(X = X, Z = Z, cpara_max = cpara_max, cpara_min = cpara_min,
                                   lower = lower, upper = upper,
                                   init = list("cpara" = cpara.init),
                                   covtype = "Matern5_2", orthogonal = FALSE, f.sim = f.sim, df.sim = df.sim)
  }
  tm.Hom[ii] <- sum(sapply(model, function(x) x$time))
  
  llmax.index <- which.max(sapply(model, function(x) x$ll))
  model <- model[[llmax.index]]
  cpara.Hom[ii] <- model$cpara
  
  # Create a prediction grid and obtain predictions
  xgrid <- matrix(seq(min(X0), max(X0), length.out = 101), ncol = 1) 
  predictions <- predict(x = xgrid, object =  model)
  rmse.Hom[ii] <- sqrt(mean((p.fun(xgrid) - predictions$mean)^2))
  var.pred <- predictions$sd2 + predictions$nugs
  score.Hom[ii] <- -mean((p.fun(xgrid) - predictions$mean)^2/var.pred) - mean(var.f(xgrid)/var.pred) - mean(log(var.pred))
  
  ## Home with orthogonality ------------------------------------------------------------
  model <- vector("list", 3)
  jj <- 0
  for(cpara.init in cpara_init.vt){
    jj <- jj + 1
    model[[jj]] <- mleHomCalibrate(X = X, Z = Z, cpara_max = cpara_max, cpara_min = cpara_min, 
                                   lower = lower, upper = upper,
                                   init = list("cpara" = cpara.init),
                                   covtype = "Matern5_2", orthogonal = TRUE, f.sim = f.sim, df.sim = df.sim)
  }
  tm.Hom.OGP[ii] <- sum(sapply(model, function(x) x$time))
  
  llmax.index <- which.max(sapply(model, function(x) x$ll))
  model <- model[[llmax.index]]
  cpara.Hom.OGP[ii] <- model$cpara
  
  # Create a prediction grid and obtain predictions
  xgrid <- matrix(seq(min(X0), max(X0), length.out = 101), ncol = 1) 
  predictions <- predict(x = xgrid, object =  model)
  rmse.Hom.OGP[ii] <- sqrt(mean((p.fun(xgrid) - predictions$mean)^2))
  var.pred <- predictions$sd2 + predictions$nugs
  score.Hom.OGP[ii] <- -mean((p.fun(xgrid) - predictions$mean)^2/var.pred) - mean(var.f(xgrid)/var.pred) - mean(log(var.pred))
  
  ## Het without orthogonality ------------------------------------------------------------
  model <- vector("list", 3)
  jj <- 0
  for(cpara.init in cpara_init.vt){
    jj <- jj + 1
    model[[jj]] <- mleHetCalibrate(X = X, Z = Z, cpara_max = cpara_max, cpara_min = cpara_min,
                                   lower = lower, upper = upper,
                                   init = list("cpara" = cpara.init),
                                   settings = list(checkHom = FALSE, linkThetas = "none"),
                                   covtype = "Matern5_2", orthogonal = FALSE, f.sim = f.sim, df.sim = df.sim)
    
  }
  tm.Het[ii] <- sum(sapply(model, function(x) x$time))
  
  
  llmax.index <- which.max(sapply(model, function(x) x$ll))
  model <- model[[llmax.index]]
  cpara.Het[ii] <- model$cpara
  
  ## Create a prediction grid and obtain predictions
  xgrid <- matrix(seq(min(X0), max(X0), length.out = 101), ncol = 1) 
  predictions <- predict(x = xgrid, object =  model)
  rmse.Het[ii] <- sqrt(mean((p.fun(xgrid) - predictions$mean)^2))
  var.pred <- predictions$sd2 + predictions$nugs
  score.Het[ii] <- -mean((p.fun(xgrid) - predictions$mean)^2/var.pred) - mean(var.f(xgrid)/var.pred) - mean(log(var.pred))
  
  
  ## Het with orthogonality ------------------------------------------------------------
  model <- vector("list", 3)
  jj <- 0
  for(cpara.init in cpara_init.vt){
    jj <- jj + 1
    model[[jj]] <- mleHetCalibrate(X = X, Z = Z, cpara_max = cpara_max, cpara_min = cpara_min,
                                   lower = lower, upper = upper,
                                   init = list("cpara" = cpara.init),
                                   settings = list(checkHom = FALSE, linkThetas = "none"),
                                   covtype = "Matern5_2", orthogonal = TRUE, f.sim = f.sim, df.sim = df.sim)
    
  }
  tm.Het.OGP[ii] <- sum(sapply(model, function(x) x$time))
  
  llmax.index <- which.max(sapply(model, function(x) x$ll))
  model <- model[[llmax.index]]
  cpara.Het.OGP[ii] <- model$cpara
  
  ### contruct confidence interval for the calibration parameter
  Info.mx <- computeInfo(model) # computer the information matrix
  sd.cpara <- sqrt(diag(Info.mx)[1])
  LCL <- qnorm(0.025, model$cpara, sd.cpara)
  UCL <- qnorm(0.975, model$cpara, sd.cpara)
  coverage[ii] <- UCL > true.cpara & LCL < true.cpara
  

  ## Create a prediction grid and obtain predictions
  xgrid <- matrix(seq(min(X0), max(X0), length.out = 101), ncol = 1) 
  predictions <- predict(x = xgrid, object =  model)
  rmse.Het.OGP[ii] <- sqrt(mean((p.fun(xgrid) - predictions$mean)^2))
  var.pred <- predictions$sd2 + predictions$nugs
  score.Het.OGP[ii] <- -mean((p.fun(xgrid) - predictions$mean)^2/var.pred) - mean(var.f(xgrid)/var.pred) - mean(log(var.pred))
}
print(max(tm.Het.OGP)) # print the maximum computational time
print(sum(coverage)/100) # print the coverage rate of the confidence intervals

# plot results
par(mar=c(5,4,0.5,0.5))
par(mfrow = c(1,3))
boxplot(list(cpara.WLS-true.cpara, cpara.Hom-true.cpara, 
             cpara.Hom.OGP-true.cpara, cpara.Het-true.cpara, cpara.Het.OGP-true.cpara), 
        ylab = "estimation bias", xaxt="n", ylim = c(-0.1,0.3)) 
axis(1, 1:5, labels = c("WLS", "HomGP", "HomOGP", "HetGP", "HetOGP"), las = 2)
abline(h = 0, col = 2)
boxplot(list(rmse.WLS, rmse.Hom, rmse.Hom.OGP, rmse.Het, rmse.Het.OGP), ylab = "RMSE", xaxt="n") 
axis(1, 1:5, labels = c("WLS", "HomGP", "HomOGP", "HetGP", "HetOGP"), las = 2)
boxplot(list(NULL, score.Hom, score.Hom.OGP, score.Het, score.Het.OGP), ylab = "score", xaxt="n") 
axis(1, 1:5, labels = c("WLS", "HomGP", "HomOGP", "HetGP", "HetOGP"), las = 2)
# 9X3







