### setting ###
n.timestep <- 8
n.rep.each <- 5

# computer model
f.sim <- function(x, cpara) {
  return(c(exp(x/10)*sin(x) - sqrt(cpara^2 - cpara + 1) * (sin(cpara*x)+cos(cpara*x))))
}
df.sim <- function(x, cpara) {
  return(c(-sqrt(cpara^2-cpara+1)*(x*cos(x*cpara)-x*sin(x*cpara))-((2*cpara-1)*(sin(x*cpara)+cos(x*cpara)))/(2*sqrt(cpara^2-cpara+1))))
}

# physical process
p.fun <- function(x) exp(x/10)*sin(x)

# true parameter
true.cpara <- optim(0, fn = function(g) {
  x.grid <- seq(0,2*pi,0.01)
  mean((p.fun(x.grid) - f.sim(x.grid, g))^2)
},
lower = -0.3, upper = 0.3, method = "L-BFGS-B")$par


# observed input
X0 <- seq(0,2*pi,length.out = n.timestep)
# mean process
pmean <- p.fun(X0)
# number of replicates
n.rep <- rep(n.rep.each,length(X0))

# setting for lower and upper bounds of parameters
lower <- 0.01*max(X0)
upper <- 2.5*max(X0)
cpara_min <- -0.3
cpara_max <- 0.3
cpara_init.vt <- c(-0.2, 0, 0.2)

empirical.power <- rep(0,3)
cpara.WLS.ls <- cpara.Hom.ls <- cpara.Hom.OGP.ls <- cpara.Het.ls <- cpara.Het.OGP.ls <- vector("list", 3)

for(kk in 1:3){
  # variance process
  if(kk == 1){
    var.f <- function(x) rep(1, length(x))
  }else if(kk == 2){
    var.f <- function(x) 3*exp(-3*(x-pi/2)^2)+2*exp(-3*(x-3*pi/2)^2)+0.01
  }else{
    var.f <- function(x) 6*exp(-6*(x-pi/2)^2)+0.01
  }

  # variance process
  var.y <- var.f(X0)

  ### simulate 100 times ###
  cpara.WLS <- cpara.Hom <- cpara.Hom.OGP <- cpara.Het <- cpara.Het.OGP <- chitest <- rep(0, 100)

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

    llmax.index <- which.max(sapply(model, function(x) x$ll))
    model <- model[[llmax.index]]
    cpara.Hom[ii] <- model$cpara

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

    llmax.index <- which.max(sapply(model, function(x) x$ll))
    model <- model[[llmax.index]]
    cpara.Hom.OGP[ii] <- model$cpara

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


    llmax.index <- which.max(sapply(model, function(x) x$ll))
    model <- model[[llmax.index]]
    cpara.Het[ii] <- model$cpara

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

    llmax.index <- which.max(sapply(model, function(x) x$ll))
    model <- model[[llmax.index]]
    cpara.Het.OGP[ii] <- model$cpara

    Info.mx <- computeInfo(model) # compute the information matrix

    nmean <- drop(rowSums(model$Kgi) %*% model$Delta / sum(model$Kgi))
    sKgi <- sum(model$Kgi)
    R <- diag(1,n.timestep) - matrix(1,n.timestep,n.timestep) %*% model$Kgi / sKgi
    chi_stat <- drop((model$Delta - nmean) %*% solve(R %*% Info.mx[(nrow(Info.mx)-length(model$Delta)+1):ncol(Info.mx), (nrow(Info.mx)-length(model$Delta)+1):ncol(Info.mx)]%*% t(R) + diag(1e-8, n.timestep))  %*% (model$Delta - nmean))
    chitest[ii] <- 1-pchisq(chi_stat, df = length(model$Delta))
  }
  empirical.power[kk] <- sum(chitest < 0.01)
  cpara.WLS.ls[[kk]] <- cpara.WLS
  cpara.Hom.ls[[kk]] <- cpara.Hom
  cpara.Hom.OGP.ls[[kk]] <- cpara.Hom.OGP
  cpara.Het.ls[[kk]] <- cpara.Het
  cpara.Het.OGP.ls[[kk]] <- cpara.Het.OGP
}


# plot results
par(mfrow=c(2,3), oma=c(0, 4, 3, 2.5), mar=c(5,0.5,0,0.5))
for(kk in 1:3){
  if(kk == 1){
    var.f <- function(x) rep(1, length(x))
    curve(var.f, 0, 2*pi, ylim = c(0,6), col = 1, lty = 2, lwd = 1, xlab = "x", ylab = "variance")
    mtext("variance", 2, 3, las = 0)
    mtext("(i)", 3, 1, las = 0)
  }else if(kk == 2){
    var.f <- function(x) 3*exp(-3*(x-pi/2)^2)+3*exp(-3*(x-3*pi/2)^2)+0.01
    curve(var.f,0, 2*pi, ylim = c(0,6), col = 1, lty = 2, lwd = 1, yaxt="n")
    mtext("(ii)", 3, 1, las = 0)
  }else{
    var.f <- function(x) 6*exp(-6*(x-pi/2)^2)+0.01
    curve(var.f,0, 2*pi, ylim = c(0,6), col = 1, lty = 2, lwd = 1, yaxt="n")
    mtext("(iii)", 3, 1, las = 0)
  }
}


for(kk in 1:3){
  if(kk == 1){
    boxplot(list(cpara.WLS.ls[[kk]]-true.cpara,
                 cpara.Hom.ls[[kk]]-true.cpara,
                 cpara.Hom.OGP.ls[[kk]]-true.cpara,
                 cpara.Het.ls[[kk]]-true.cpara,
                 cpara.Het.OGP.ls[[kk]]-true.cpara),
            ylab = "estimation bias", xaxt="n", ylim = c(-0.05,0.2))
    mtext("estimation bias", 2, 3, las = 0)
  }else{
    boxplot(list(cpara.WLS.ls[[kk]]-true.cpara,
                 cpara.Hom.ls[[kk]]-true.cpara,
                 cpara.Hom.OGP.ls[[kk]]-true.cpara,
                 cpara.Het.ls[[kk]]-true.cpara,
                 cpara.Het.OGP.ls[[kk]]-true.cpara), xaxt="n", yaxt="n", ylim = c(-0.05,0.2))
  }
  axis(1, 1:5, labels = c("WLS", "HomGP", "HomOGP", "HetGP", "HetOGP"), las = 2)
  abline(h = 0, col = 2)

}

## MSE of theta hat
print(sapply(cpara.Hom.OGP.ls, FUN = function(x.vt) mean((x.vt-true.cpara)^2))*10000)
print(sapply(cpara.Het.OGP.ls, FUN = function(x.vt) mean((x.vt-true.cpara)^2))*10000)
## empirical power
print(empirical.power)
#7 X 4.5

