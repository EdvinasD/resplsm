rm(list=ls())
library(resplsm)
library(ggplot2)
library(reshape2)
source('~/GitHub/resplsm/R/test_script_functions.R')
# parameters for function:
# AR1 ARCH1 y_t = 0.5y_{t-1}+sqrt{0.05+0.1y^2_{t-1}}\varepsilon_t
burnperiod <- 10
N <- 1000+burnperiod
theta1 <-  0.5
theta2 <-  0.05
theta3 <- 0.1
theta_true <- c(theta1,theta2,theta3)



# with gaussian inovations
set.seed(95)
Y <- rep(0,N)
e <- rnorm(N)
# e[runif(N)<0.025]=5
Y[1] <- e[1]
for (i in 2:N){
  Y[i] <- theta1*Y[i-1] + sqrt(theta2 + theta3*Y[i-1]^2)*e[i]

}
Y <- Y[(burnperiod+1):N]
N <- length(Y)
X <- Y[-N]
Y <- Y[-1]

plot(Y,type='l')


n.points.X <- 10
X.theta.points=sort(X)[round(seq(0,1,length.out=n.points.X+2)[-c(1,n.points.X+2)]*N)]
X.theta.points=seq(from=X.theta.points[1], to=X.theta.points[n.points.X], length.out = n.points.X )

semi.parametric <- espls(Yt=Y, St=X,
                       s=X.theta.points,
                       initial.values=c(0.5,0.5),
                       bandwidth=0.15)


robust.semi.parametric <- respls(theta=semi.parametric,
                           Y=Y,
                           X=X,
                           c_bound=4,
                           iterations=5,
                           bindwidths=0.15)

getfunctionform(X.theta.points,theta_true)



# Estimation using parametric texniques
mu_thetaAR1 <- function(x,theta){
  theta[1]*x
}

sigma_thetaAR1 <- function(x,theta){
  if(theta[2]<0){
    theta[2]=NA
  }
  if(theta[3]<0){
    theta[3]=NA
  }
  theta[2]+theta[3]*x^2
}

d_mu_thetaAR1 <- function(x,theta){
  rbind(x,
        rep(0,length(x)),
        rep(0,length(x)))
}

d_sigma_thetaAR1 <- function(x,theta){
  theta[2] <- max(0,theta[2])
  theta[3] <- max(0,theta[3])

  # sigma_simp <- 1/(2*sqrt(theta[2]+theta[3]*x^2))
  #
  # rbind(rep(0,length(x)),
  #       sigma_simp,
  #       x^2*sigma_simp)
  rbind(rep(0,length(x)),
        rep(1,length(x)),
        x^2)

}

initial.theta <- c(0.25,0.25,0.25)


theta.pls <- pls(initial.theta=initial.theta,
                 Y=Y, X=X,
                 func_s=sigma_thetaAR1,
                 func_m=mu_thetaAR1)


theta.rls <- rls(theta = theta.pls,
                 Y = Y,
                 X = X,
                 c_bound = 7,
                 iterations = 20,
                 func_s = sigma_thetaAR1,
                 d_func_s = d_sigma_thetaAR1,
                 func_m = mu_thetaAR1,
                 d_func_m = d_mu_thetaAR1,return.all = T,
                 tolerance = 1)

plot.respls(rbind.thetas(list(RSMPML=robust.semi.parametric,
                              SPPML=semi.parametric,
                              PML=funform(X.theta.points,theta.pls),
                              RML=funform(X.theta.points,theta.rls[[20]]))
                         ),
            getfunctionform(X.theta.points,theta_true))

