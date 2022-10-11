##### heteroscedastic test #####
empirical.power.mx <- matrix(0, ncol=3, nrow=4)
samplesize.mx <- matrix(c(8,8,16,16,5,10,10,20), ncol=2)

for(aa in 1:nrow(samplesize.mx)){
  n.timestep <- samplesize.mx[aa,1]
  n.rep.each <- samplesize.mx[aa,2]

  cat("n.timestep=", n.timestep, ", n.rep.each=", n.rep.each,"\n")

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
    cpara.WLS <- cpara.Hom <- cpara.Hom.OGP <- cpara.Het <- cpara.Het.OGP <- het.test <- rep(0, 1000)

    for(ii in 1:1000){
      print(ii)
      set.seed(ii)

      # simulate X and Z
      X <- matrix(rep(X0, n.rep), ncol = 1)
      Z <- rep(0, sum(n.rep))
      for(i in 1:length(X0)) {
        Z[(ifelse(i==1,0,sum(n.rep[1:(i-1)]))+1):sum(n.rep[1:i])] <- pmean[i] + rnorm(n.rep[i], 0, sd = sqrt(var.y[i]))
      }

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
      chi_stat <- drop((model$Delta - nmean) %*% solve(R %*% Info.mx[(nrow(Info.mx)-n.timestep+1):ncol(Info.mx), (nrow(Info.mx)-n.timestep+1):ncol(Info.mx)]%*% t(R) + diag(1e-8, n.timestep))  %*% (model$Delta - nmean))
      het_stat <- (chi_stat - n.timestep)/sqrt(2*n.timestep)
      het.test[ii] <- 1-pnorm(het_stat)
    }
    empirical.power[kk] <- sum(het.test < 0.05)
  }
  empirical.power.mx[aa,] <- empirical.power
  print(empirical.power)
}
