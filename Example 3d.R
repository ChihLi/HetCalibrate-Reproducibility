### setting ###
n.timestep <- 30

f.sim <- function(x, cpara) {
  
  if(is.null(dim(x))){
    x1 <- x[1]
    x2 <- x[2]
  }else{
    x1 <- x[,1]
    x2 <- x[,2]
  }  
  
  out <- cpara[1] + x1*cpara[2] + x2*cpara[3]
  return(c(out))
}
df.sim <- function(x, cpara) {
  if(is.null(dim(x))){
    return(c(1,x[,1],x[,2]))
  }else{
    return(cbind(1,x))
  }
}


# variance process
var.f <- function(x) 0.01*exp(-10*sin(x[,1]*pi) * cos(x[,2]*pi))

# physical process
p.fun <- function(x) {
  if(is.null(dim(x))){
    x1 <- x[1]
    x2 <- x[2]
  }else{
    x1 <- x[,1]
    x2 <- x[,2]
  }  
  4*x1+x1*sin(5*x2)
}


# numerical solution
# true.cpara <- optim(c(0,0,0), fn = function(g) {
#   x.grid <- as.matrix(expand.grid(seq(0,1,0.01),seq(0,1,0.01)))
#   mean((p.fun(x.grid) - f.sim(x.grid, g))^2)
# }, 
# lower = -5, upper = 5, method = "L-BFGS-B")$par

# true parameters
A <- matrix(0,3,3)
b <- rep(0,3)
A[1,] <- c(10,5,5)
b[1] <- 21-cos(5)
A[2,] <- c(300,200,150)
b[2] <- 840-40*cos(5)
A[3,] <- c(300,150,200)
b[3] <- 600-60*cos(5)+12*sin(5)
true.cpara <- solve(A,b)




# setting for lower and upper bounds of parameters
lower <- 0.01
upper <- 2.5
cpara_min <- rep(-5, 3)  
cpara_max <- rep(5, 3)
cpara_init.mx <- as.matrix(expand.grid(c(-2.5,0,2.5),c(-2.5,0,2.5),c(-2.5,0,2.5)))
colnames(cpara_init.mx) <- NULL

### simulate 100 times ###
cpara.WLS <- cpara.Hom <- cpara.Hom.OGP <- cpara.Het <- cpara.Het.OGP <- matrix(0, ncol = 3, nrow = 100) 

ll.Hom <- ll.Hom.OGP <- ll.Het <- ll.Het.OGP <- 
  tm.Hom <- tm.Hom.OGP <- tm.Het <- tm.Het.OGP <- 
  rmse.WLS <- rmse.Hom <- rmse.Hom.OGP <- rmse.Het <- rmse.Het.OGP <- 
  score.Hom <- score.Hom.OGP <- score.Het <- score.Het.OGP <- rep(0, 100)

for(n.rep.each in c(5, 10, 2)){
  n.rep <- rep(n.rep.each, n.timestep)
  cat("n.rep.each:", n.rep.each, "\n")
  for(ii in 1:100){
    print(ii)
    set.seed(ii)
    
    # observed input
    X0 <- maximinSLHD(1, n.timestep, 2)$StandDesign
    # mean process
    pmean <- p.fun(X0)
    # variance process
    var.y <- var.f(X0)
    
    # simulate X and Z
    X <- matrix(rep(X0, rep(n.rep, ncol(X0))), ncol = 2)
    Z <- rep(0, sum(n.rep))
    for(i in 1:nrow(X0)) {
      Z[(ifelse(i==1,0,sum(n.rep[1:(i-1)]))+1):sum(n.rep[1:i])] <- pmean[i] + rnorm(n.rep[i], 0, sd = sqrt(var.y[i]))
    }
    
    ## WLS estimate ------------------------------------------------------------ 
    Z0 <- hetGP::find_reps(X, Z)$Z0
    Sigma_inv <- diag(1/sapply(hetGP::find_reps(X, Z)$Zlist, var))
    min.index <- which.min(apply(cpara_init.mx, 1, function(x) optim(x, fn = function(g) t(Z0 - f.sim(X0, g)) %*% Sigma_inv %*% (Z0 - f.sim(X0, g)), 
                                                                     lower = cpara_min, upper = cpara_max, method = "L-BFGS-B")$value))
    
    cpara.WLS[ii,] <- optim(cpara_init.mx[min.index,], 
                            fn = function(g) t(Z0 - f.sim(X0, g)) %*% Sigma_inv %*% (Z0 - f.sim(X0, g)), 
                            lower = cpara_min, upper = cpara_max, method = "L-BFGS-B")$par
    
    ## Hom without orthogonality ------------------------------------------------------------
    model <- vector("list", nrow(cpara_init.mx))
    for(jj in 1:nrow(cpara_init.mx)){
      cpara.init.vt <- cpara_init.mx[jj,]
      model[[jj]] <- mleHomCalibrate(X = X, Z = Z, cpara_max = cpara_max, cpara_min = cpara_min,
                                     lower = lower, upper = upper,
                                     inputBounds = matrix(c(0,0,1,1), ncol = 2, byrow = TRUE),
                                     init = list("cpara" = cpara.init.vt),
                                     covtype = "Matern5_2", orthogonal = FALSE, f.sim = f.sim, df.sim = df.sim)
    }
    tm.Hom[ii] <- sum(sapply(model, function(x) x$time))
    
    llmax.index <- which.max(sapply(model, function(x) x$ll))
    model <- model[[llmax.index]]
    cpara.Hom[ii,] <- model$cpara
    
    ## Home with orthogonality ------------------------------------------------------------
    model <- vector("list", nrow(cpara_init.mx))
    for(jj in 1:nrow(cpara_init.mx)){
      cpara.init.vt <- cpara_init.mx[jj,]
      model[[jj]] <- mleHomCalibrate(X = X, Z = Z, cpara_max = cpara_max, cpara_min = cpara_min, 
                                     lower = lower, upper = upper,
                                     inputBounds = matrix(c(0,0,1,1), ncol = 2, byrow = TRUE),
                                     init = list("cpara" = cpara.init.vt),
                                     covtype = "Matern5_2", orthogonal = TRUE, f.sim = f.sim, df.sim = df.sim)
    }
    tm.Hom.OGP[ii] <- sum(sapply(model, function(x) x$time))
    
    llmax.index <- which.max(sapply(model, function(x) x$ll))
    model <- model[[llmax.index]]
    cpara.Hom.OGP[ii,] <- model$cpara
    
    ## Het without orthogonality ------------------------------------------------------------
    model <- vector("list", nrow(cpara_init.mx))
    for(jj in 1:nrow(cpara_init.mx)){
      cpara.init.vt <- cpara_init.mx[jj,]
      model[[jj]] <- mleHetCalibrate(X = X, Z = Z, cpara_max = cpara_max, cpara_min = cpara_min,
                                     lower = lower, upper = upper,
                                     inputBounds = matrix(c(0,0,1,1), ncol = 2, byrow = TRUE),
                                     init = list("cpara" = cpara.init.vt),
                                     settings = list(checkHom = FALSE, linkThetas = "none"),
                                     covtype = "Matern5_2", orthogonal = FALSE, f.sim = f.sim, df.sim = df.sim)
      
    }
    tm.Het[ii] <- sum(sapply(model, function(x) x$time))
    
    llmax.index <- which.max(sapply(model, function(x) x$ll))
    model <- model[[llmax.index]]
    cpara.Het[ii,] <- model$cpara
    
    ## Het with orthogonality ------------------------------------------------------------
    model <- vector("list", 3)
    for(jj in 1:nrow(cpara_init.mx)){
      cpara.init.vt <- cpara_init.mx[jj,]
      model[[jj]] <- mleHetCalibrate(X = X, Z = Z, cpara_max = cpara_max, cpara_min = cpara_min,
                                     lower=lower, upper = upper, 
                                     inputBounds = matrix(c(0,0,1,1), ncol = 2, byrow = TRUE),
                                     init = list("cpara" = cpara.init.vt),
                                     settings = list(checkHom = FALSE, linkThetas = "none"),
                                     covtype = "Matern5_2", orthogonal = TRUE, f.sim = f.sim, df.sim = df.sim)
      
    }
    tm.Het.OGP[ii] <- sum(sapply(model, function(x) x$time))
    
    llmax.index <- which.max(sapply(model, function(x) x$ll))
    model <- model[[llmax.index]]
    cpara.Het.OGP[ii,] <- model$cpara
    
    # plot results
    par(mar=c(5,4,0.5,0.5))
    par(mfrow = c(1,3))
    for(i in 1:3){
      boxplot(list(cpara.WLS[1:ii,i]-true.cpara[i], 
                   cpara.Hom[1:ii,i]-true.cpara[i], 
                   cpara.Hom.OGP[1:ii,i]-true.cpara[i], 
                   cpara.Het[1:ii,i]-true.cpara[i], 
                   cpara.Het.OGP[1:ii,i]-true.cpara[i]), 
              ylab = "estimation bias", xaxt="n") 
      axis(1, 1:5, labels = c("WLS", "HomGP", "HomOGP", "HetGP", "HetOGP"), las = 2)
      abline(h = 0, col = 2)
    }
  }
  par(mar=c(5,4,0.5,0.5))
  par(mfrow = c(1,3))
  for(i in 1:3){
    boxplot(list(cpara.WLS[1:ii,i]-true.cpara[i], 
                 cpara.Hom[1:ii,i]-true.cpara[i], 
                 cpara.Hom.OGP[1:ii,i]-true.cpara[i], 
                 cpara.Het[1:ii,i]-true.cpara[i], 
                 cpara.Het.OGP[1:ii,i]-true.cpara[i]), 
            ylab = "estimation bias", xaxt="n", ylim=c(-1,1)) 
    axis(1, 1:5, labels = c("WLS", "HomGP", "HomOGP", "HetGP", "HetOGP"), las = 2)
    abline(h = 0, col = 2)
  }
  save(cpara.WLS, cpara.Hom, cpara.Hom.OGP, cpara.Het, cpara.Het.OGP, file = paste0("example_3d_n_", n.timestep, "_rep_", n.rep.each, ".Rdata"))
}



load("example_3d_n_30_rep_2.Rdata")
list.rep2 <- list(cpara.WLS, cpara.Hom, cpara.Hom.OGP, cpara.Het, cpara.Het.OGP)
load("example_3d_n_30_rep_5.Rdata")
list.rep5 <- list(cpara.WLS, cpara.Hom, cpara.Hom.OGP, cpara.Het, cpara.Het.OGP)
load("example_3d_n_30_rep_10.Rdata")
list.rep10 <- list(cpara.WLS, cpara.Hom, cpara.Hom.OGP, cpara.Het, cpara.Het.OGP)

par(mfrow=c(1,3), oma=c(7, 4, 1, 2.5), mar=c(0.5,0.5,1.5,0.5))
for(i in 1:3){
  theta_bias <- c(lapply(list.rep2, function(x.ls) x.ls[,i] - true.cpara[i]),
                  lapply(list.rep5, function(x.ls) x.ls[,i] - true.cpara[i]),
                  lapply(list.rep10, function(x.ls) x.ls[,i] - true.cpara[i]))
  if(i == 1) boxplot(theta_bias, xaxt="n", ylim = c(-2,2)) else boxplot(theta_bias, yaxt = "n", xaxt="n", ylim = c(-2,2))
  if(i == 1) mtext("estimation bias", 2, 3, las = 0)
  abline(h = 0, col = 2)
  abline(v = c(5.5, 10.5), lty = 3)
  axis(1, at = c(5.5, 10.5), labels = NA, tck = -0.4, lty = 3)
  axis(1, 1:15, labels = rep(c("WLS", "HomGP", "HomOGP", "HetGP", "HetOGP"), 3), las = 2, cex = 0.7)
  axis(1, c(3,8, 13), labels = c("rep.=2", "rep.=5", "rep.=8"), las = 1, tick = FALSE, padj=7)
  if(i == 1) mtext(expression(theta[1]), 3, las = 0, padj = -0.5)
  if(i == 2) mtext(expression(theta[2]), 3, las = 0, padj = -0.5)
  if(i == 3) mtext(expression(theta[3]), 3, las = 0, padj = -0.5)
}
#9X4
