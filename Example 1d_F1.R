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

# variance process
var.f <- function(x) (0.01+0.2*(x-pi)^2)^2

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
# variance process
var.y <- var.f(X0)
# number of replicates
n.rep <- rep(n.rep.each,length(X0))

# setting for lower and upper bounds of parameters
lower <- 0.01*max(X0)
upper <- 2.5*max(X0)
cpara_min <- -0.3  
cpara_max <- 0.3
cpara_init.vt <- c(-0.2, 0, 0.2)


### plot one example 
set.seed(1)

# simulate X and Z
X <- matrix(rep(X0, n.rep), ncol = 1)
Z <- rep(0, sum(n.rep))
for(i in 1:length(X0)) {
  Z[(ifelse(i==1,0,sum(n.rep[1:(i-1)]))+1):sum(n.rep[1:i])] <- pmean[i] + rnorm(n.rep[i], 0, sd = sqrt(var.y[i]))
}

## WLS
X0 <- hetGP::find_reps(X, Z)$X0
Z0 <- hetGP::find_reps(X, Z)$Z0
Sigma_inv <- diag(1/sapply(hetGP::find_reps(X, Z)$Zlist, var))
min.index <- which.min(apply(matrix(seq(cpara_min, cpara_max, length.out = 11),ncol=1), 1, function(x) optim(x, fn = function(g) t(Z0 - f.sim(X0, g)) %*% Sigma_inv %*% (Z0 - f.sim(X0, g)), 
                                                                                                             lower = cpara_min, upper = cpara_max, method = "L-BFGS-B")$value))

cpara.WLS <- optim(seq(cpara_min, cpara_max, length.out = 11)[min.index], 
                   fn = function(g) t(Z0 - f.sim(X0, g)) %*% Sigma_inv %*% (Z0 - f.sim(X0, g)), 
                   lower = cpara_min, upper = cpara_max, method = "L-BFGS-B")$par
xgrid <- matrix(seq(min(X0), max(X0), length.out = 101), ncol = 1) 
predictions.WLS <- f.sim(xgrid, cpara.WLS)
cat("============= WLS parameter:", cpara.WLS, "=============\n")

## KO
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
cpara.Hom <- model$cpara
predictions.Hom <- predict(x = xgrid, object =  model)
cat("============= HomGP parameter:", cpara.Hom, "=============\n")

## our model
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
cpara.Het.OGP <- model$cpara
predictions.Het <- predict(x = xgrid, object =  model)
cat("============= HetOGP parameter:", cpara.Het.OGP, "=============\n")

par(mfrow=c(2,3), oma=c(4.5, 4, 0, 2.5), mar=c(0.5,0.5,1.5,0.5))
## mean process
# WLS
plot(X, Z, ylab = 'y', xlab = "", xaxt="n")
mtext("y", 2, 3, las = 0)
axis(1, seq(0,6,1), labels = rep("", 7))
points(X0, Z0, pch = 20, cex = 1.2)
lines(xgrid, predictions.WLS, col = 'red', lwd = 2)
curve(p.fun, min(X0), max(X0), add = TRUE, col = 1, lty = 2, lwd = 1)
lines(xgrid, f.sim(xgrid, cpara.WLS), col = 4, lty = 2, lwd = 2)

# KO
plot(X, Z, ylab = 'y', xlab = "", xaxt="n", yaxt="n")
axis(1, seq(0,6,1), labels = rep("", 7))
points(X0, Z0, pch = 20, cex = 1.2)
lines(xgrid, predictions.Hom$mean, col = 'red', lwd = 2)
curve(p.fun, min(X0), max(X0), add = TRUE, col = 1, lty = 2, lwd = 1)
lines(xgrid, qnorm(0.025, predictions.Hom$mean, sqrt(predictions.Hom$sd2 + predictions.Hom$nugs)), 
      col = 3, lty = 3, lwd = 2)
lines(xgrid, qnorm(0.975, predictions.Hom$mean, sqrt(predictions.Hom$sd2 + predictions.Hom$nugs)), 
      col = 3, lty = 3, lwd = 2)
lines(xgrid, f.sim(xgrid, cpara.Hom), col = 4, lty = 2, lwd = 2)

# our model
plot(X, Z, ylab = 'y', xlab = "", xaxt="n", yaxt="n")
axis(1, seq(0,6,1), labels = rep("", 7))
points(X0, Z0, pch = 20, cex = 1.2)
lines(xgrid, predictions.Het$mean, col = 'red', lwd = 2)
curve(p.fun, min(X0), max(X0), add = TRUE, col = 1, lty = 2, lwd = 1)
lines(xgrid, qnorm(0.025, predictions.Het$mean, sqrt(predictions.Het$sd2 + predictions.Het$nugs)), 
      col = 3, lty = 3, lwd = 2)
lines(xgrid, qnorm(0.975, predictions.Het$mean, sqrt(predictions.Het$sd2 + predictions.Het$nugs)), 
      col = 3, lty = 3, lwd = 2)
lines(xgrid, f.sim(xgrid, cpara.Het.OGP), col = 4, lty = 2, lwd = 2)

## variance process
# WLS
curve(var.f, min(X0), max(X0), col = 1, lty = 2, lwd = 1, xlab = "x", ylab = "variance")
mtext("variance", 2, 3, las = 0)
points(X0, sapply(unstack(data.frame(Z,rep(1:length(X0),n.rep))), var), col = 1, pch = 20, cex = 1.2)
mtext("x", 1, 3, las = 0)
# KO
curve(var.f, min(X0), max(X0), col = 1, lty = 2, lwd = 1, xlab = "x", yaxt="n", ylab = "variance")
points(X0, sapply(unstack(data.frame(Z,rep(1:length(X0),n.rep))), var), col = 1, pch = 20, cex = 1.2)
lines(xgrid, predictions.Hom$nugs, col = 2, lwd = 2)
mtext("x", 1, 3, las = 0)
# our model
curve(var.f, min(X0), max(X0), col = 1, lty = 2, lwd = 1, xlab = "x", yaxt="n", ylab = "variance")
points(X0, sapply(unstack(data.frame(Z,rep(1:length(X0),n.rep))), var), col = 1, pch = 20, cex = 1.2)
lines(xgrid, predictions.Het$nugs, col = 2, lwd = 2)
mtext("x", 1, 3, las = 0)

