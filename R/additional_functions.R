#' Estimating function
#'
#' @param yt A number.
#' @param thetas A vector of lengths 2.
#' @return A value of score function.
semi_est_func <- function(yt,thetas){
  thetas[2] <- max(0,thetas[2]^2)
  -0.5*log(2*pi)-0.5*(log(thetas[2]))-0.5*(yt-thetas[1])^2/thetas[2]
}

#' Residuals
#'
#' @param y A number.
#' @param x A number/vectors.
#' @param theta A vector of lengths 2 of data.frame, depends on func_mu and
#' func_sigma.
#' @param func_mu location function
#' @param func_sigma scale function
#' @return residual
u_resids <- function(y,x,theta,func_mu,func_sigma){
  (y-func_mu(x,theta))/sqrt(func_sigma(x,theta))
}

#' k1
#'
#' \deqn{k_1=\left.\frac{1}{2v_0(y_{t-1})^2}\frac{\partial v_0(y_{t-1})^2}{
#' \partial\theta(y_{t-1})}\right|_{\theta=\theta_0}}
#' @param x A number/vectors.
#' @param theta A data.frame
#' @param func scale function
#' @return A value of score function.
k1 <- function(theta,x,func){
  rbind(0,(1/sqrt(func(x,theta))))
}

#' k2
#'
#' @param x A number/vectors.
#' @param theta A data.frame
#' @param func scale function
#' @return A value of score function.
k2 <- function(theta,x,func){
  rbind((1/sqrt(func(x,theta))),0)
}

#' qn_function
#'
#' \eqn{q_n}
#'
#' @param u A number
#' @param parameters A list with given parameters to function:
#' \code{k1m, k2m, A, cb, tau, func_mu, func_sigma, x, theta0}
#' @return value of \eqn{q_n} function
qn_function <- function(u,parameters = list()){
  k1m <- parameters$k1m
  k2m <- parameters$k2m
  A <- parameters$A
  cb <- parameters$cb
  tau <- parameters$tau
  func_mu <- parameters$func_mu
  func_sigma <- parameters$func_sigma
  x <- parameters$x
  theta0 <- parameters$theta0


  Ymus <- func_mu(x,theta0)+sqrt(func_sigma(x,theta0))*u
  resids <- u_resids(Ymus,x,theta0,func_mu,func_sigma)
  score_v <- -k1m+k2m*resids+k1m*resids^2
  (-k1m+k2m*u+k1m*u^2)*cb/norm(A%*%(score_v-tau),type="2")
}


#'qn_function_z
#'
#' \eqn{\underline{q_n}}
#' @inheritParams qn_function
#' @param z A number
#' @return value of q_n function
qn_function_z <- function(z,u,parameters){
  qn_function(u+z,parameters)*exp(-.5*z^2)
}


Laplace_approx2 <- function(u,parameters,h=0.0001){
  q0 <- qn_function_z(0,u,parameters)
  q0ph <- qn_function_z(0+h,u,parameters)
  q0mh <- qn_function_z(0-h,u,parameters)
  q1 <-(q0ph-q0mh)/(2*h)
  q2 <-(q0ph-2*q0+q0mh)/(h^2)
  Laporxx <- q0+q1/u+q2/(u^2)
  return(Laporxx)
}

qn_function_den <- function(u,parameters = list()){
  k1m <- parameters$k1m
  k2m <- parameters$k2m
  A <- parameters$A
  cb <- parameters$cb
  tau <- parameters$tau
  func_mu <- parameters$func_mu
  func_sigma <- parameters$func_sigma
  x <- parameters$x
  theta0 <- parameters$theta0


  Ymus <- func_mu(x,theta0)+sqrt(func_sigma(x,theta0))*u
  resids <- u_resids(Ymus,x,theta0,func_mu=func_mu,func_sigma = func_sigma)
  score_v <- -k1m+k2m*resids+k1m*resids^2
  cb/norm(A%*%(score_v-tau),type="2")
}

#' Get \eqn{\theta(X)}
#'
#' @param dat data.frame which contains X and value of thate for that X
#' @param X vector for which new values of \eqn{\theta(X)} should be returned
#' @return Returns \eqn{\theta(X)}
interpolate_theta <- function(dat,X){
  splined <- spline(x=dat$X,y=dat$value,xout = X,method = 'natural')$y
  return(splined)
}

#'qn_function_z_den
#'
#' \eqn{\underline{q_n}}
#' @inheritParams qn_function
#' @param z A number
#' @return value of q_n function
qn_function_z_den <- function(z,u,parameters){
  qn_function_den(u+z,parameters)*exp(-.5*z^2)
}

Laplace_approx2_den <- function(u,parameters,h=0.0001){
  q0 <- qn_function_z_den(0,u,parameters)
  q0ph <- qn_function_z_den(0+h,u,parameters)
  q0mh <- qn_function_z_den(0-h,u,parameters)
  q1 <-(q0ph-q0mh)/(2*h)
  q2 <-(q0ph-2*q0+q0mh)/(h^2)
  Laporxx <- q0+q1/u+q2/(u^2)
  return(Laporxx)
}

sp_robust_est_function <- function(theta,params){
  func_mu=params$func_mu
  func_sigma=params$func_sigma
  tau_vector=params$tau_vector
  A_0=params$A_matrix
  Y=params$Y
  X=params$X
  cb=params$cb
  bindwidths=params$bindwidths
  S=params$S
  theta.vec.initial <- list()
  for(j in 1:length(theta))
    theta.vec.initial[[j]] <- data.frame(X=S,value=as.matrix(theta)[j,])

  theta=theta.vec.initial

  k1_0 <- k1(theta=theta,x=X,func=func_sigma)
  k2_0 <- k2(theta=theta,x=X,func=func_sigma)
  u_resids_0 <- u_resids(Y,X,theta,func_mu,func_sigma)
  score_0 <- -t(k1_0)+t(k2_0)*u_resids_0+t(k1_0)*u_resids_0^2
  score_tau <- score_0-t(tau_vector)

  weights <-t(matrix(dnorm((X-S)/bindwidths)))[c(1,1),]

  resultigvalues <- rowSums(weights*apply(t(score_tau),2,function(x){
    x*min(1,cb/norm(A_0%*%x,type="2"))}))
  return(resultigvalues)
}

sp_func_to_minimize  <- function(theta.local,params){
  sum(abs(sp_robust_est_function(theta.local,params)))
}

robust_est_function <- function(theta,params){
  func_mu=params$func_mu
  func_sigma=params$func_sigma
  d_func_mu=params$d_func_mu
  d_func_sigma=params$d_func_sigma
  tau_vector=params$tau_vector
  A_0=params$A_matrix
  Y=params$Y
  X=params$X
  cb=params$cb
  k1m = params$k1
  k2m = params$k2
  k1_0 <- k1m(theta=theta,x=X,func=func_sigma,d_func_sigma)
  k2_0 <- k2m(theta=theta,x=X,func=func_sigma,d_func=d_func_mu)



  u_resids_0 <- u_resids(Y,X,theta,func_mu,func_sigma)
  score_0 <- -t(k1_0)+t(k2_0)*u_resids_0+t(k1_0)*u_resids_0^2
  score_tau <- score_0-t(tau_vector)
  resultigvalues <- rowSums(apply(t(score_tau),2,function(x){
    x*min(1,cb/norm(A_0%*%x,type="2"))}))
  return(resultigvalues)
}

func_to_minimize <- function(theta,params){
  sum(abs(robust_est_function(theta,params)))
}

tau_for_2_roots <- function(qn_params,i_root_dw,i_root_up){
  # calculating tau_num
  int_u_up_inf <- 1/sqrt(2*pi)*exp(-.5*i_root_up^2)/i_root_up*
    Laplace_approx2(i_root_up,qn_params,h=0.0001)
  n_u <- list(dnorm_up=dnorm(i_root_up),
              dnorm_dw=dnorm(i_root_dw),
              pnorm_up=pnorm(i_root_up),
              pnorm_dw=pnorm(i_root_dw))
  M1m=n_u$dnorm_dw-n_u$dnorm_up
  M2m=i_root_dw*n_u$dnorm_dw-
    i_root_up*n_u$dnorm_up+
    n_u$pnorm_up-n_u$pnorm_dw
  int_u_down_u_up <- -qn_params$k1m*(n_u$pnorm_up-n_u$pnorm_dw)+
    qn_params$k2m*M1m+
    qn_params$k1m*M2m

  int_inf_u_dw <- -1/sqrt(2*pi)*exp(-.5*i_root_dw^2)/i_root_dw*
    Laplace_approx2(i_root_dw,qn_params,h=0.0001)
  tau_num=int_u_up_inf+int_u_down_u_up+int_inf_u_dw

  # calculating tau_den

  int_u_up_inf <- 1/sqrt(2*pi)*exp(-.5*i_root_up^2)/i_root_up*
    Laplace_approx2_den(i_root_up,qn_params,h=0.0001)
  int_u_down_u_up <- n_u$pnorm_up-n_u$pnorm_dw

  int_inf_u_dw <- -1/sqrt(2*pi)*exp(-.5*i_root_dw^2)/i_root_dw*
    Laplace_approx2_den(i_root_dw,qn_params,h=0.0001)
  tau_den=int_u_up_inf+int_u_down_u_up+int_inf_u_dw

  tau=tau_num/tau_den
  return(tau)
}


tau_for_4_roots_approximation <- function(qn_params,u_roots_i){
  # bad result for approximation not used
  i_root_up2 <- u_roots_i[4]
  i_root_up1 <- u_roots_i[3]
  i_root_dw1 <- u_roots_i[2]
  i_root_dw2 <- u_roots_i[1]


  # calculating tau_num
  int_u_up2_inf <- 1/sqrt(2*pi)*exp(-.5*i_root_up2^2)/i_root_up2*
    Laplace_approx2(i_root_up2,qn_params,h=0.001)

  n_u_up <- list(dnorm_up=dnorm(i_root_up2),
                 dnorm_dw=dnorm(i_root_up1),
                 pnorm_up=pnorm(i_root_up2),
                 pnorm_dw=pnorm(i_root_up1))
  M1m=n_u_up$dnorm_dw-n_u_up$dnorm_up
  M2m=i_root_up1*n_u_up$dnorm_dw-
    i_root_up2*n_u_up$dnorm_up+
    n_u_up$pnorm_up-n_u_up$pnorm_dw
  int_u_up1_u_up2 <- -qn_params$k1m*(n_u_up$pnorm_up-n_u_up$pnorm_dw)+
    qn_params$k2m*M1m+
    qn_params$k1m*M2m

  n_u_dw <- list(dnorm_up=dnorm(i_root_dw1),
                 dnorm_dw=dnorm(i_root_dw2),
                 pnorm_up=pnorm(i_root_dw1),
                 pnorm_dw=pnorm(i_root_dw2))
  M1m=n_u_dw$dnorm_dw-n_u_dw$dnorm_up
  M2m=i_root_dw2*n_u_dw$dnorm_dw-
    i_root_dw1*n_u_dw$dnorm_up+
    n_u_dw$pnorm_up-n_u_dw$pnorm_dw
  int_u_dw2_u_dw1 <- -qn_params$k1m*(n_u_dw$pnorm_up-n_u_dw$pnorm_dw)+
    qn_params$k2m*M1m+
    qn_params$k1m*M2m


  int_inf_u_dw2 <- -1/sqrt(2*pi)*exp(-.5*i_root_dw2^2)/i_root_dw2*
    Laplace_approx2(i_root_dw2,qn_params,h=0.0001)

  # qn_function1 <- function(u,parameters){qn_function(u,parameters)[1]*dnorm(u)}
  # qn_function2 <- function(u,parameters){qn_function(u,parameters)[2]*dnorm(u)}
  # c(integrate(Vectorize(qn_function1,vectorize.args='u'), lower = i_root_up2, upper = Inf,parameters = qn_params)$value,
  #   integrate(Vectorize(qn_function2,vectorize.args='u'), lower = i_root_up2, upper = Inf,parameters = qn_params)$value)
  #
  #
  int_u_dw1_u_up1 <- 0

  tau_num=int_inf_u_dw2+int_u_dw2_u_dw1+
    int_u_dw1_u_up1+
    int_u_up1_u_up2+int_u_up2_inf

  # calculating tau_den

  int_u_up2_inf <- 1/sqrt(2*pi)*exp(-.5*i_root_up2^2)/i_root_up2*
    Laplace_approx2_den(i_root_up2,qn_params,h=0.0001)

  int_u_up1_u_up2 <- n_u_up$pnorm_up-n_u_up$pnorm_dw

  int_u_dw1_u_up1 <- 0

  int_u_dw2_u_dw1 <- n_u_dw$pnorm_up-n_u_dw$pnorm_dw

  int_inf_u_dw2 <- -1/sqrt(2*pi)*exp(-.5*i_root_dw2^2)/i_root_dw2*
    Laplace_approx2_den(i_root_dw2,qn_params,h=0.0001)
  tau_den=int_inf_u_dw2+int_u_dw2_u_dw1+
    int_u_dw1_u_up1+
    int_u_up1_u_up2+int_u_up2_inf


  tau=tau_num/tau_den
  return(tau)
}



tau_for_4_roots <- function(qn_params,u_roots_i){

  i_root_up2 <- u_roots_i[4]
  i_root_up1 <- u_roots_i[3]
  i_root_dw1 <- u_roots_i[2]
  i_root_dw2 <- u_roots_i[1]

  qn_function.ind <- function(u,parameters,ind){
    (-parameters$k1m[ind]+parameters$k2m[ind]*u+parameters$k1m[ind]*u^2)*
      min(1,qn_function_den(u,parameters))*dnorm(u)}

  tau_num <- c()
  for(ind in 1:length(qn_params$theta0)){
    tau_num <- c(tau_num,
                 integrate(Vectorize(qn_function.ind,vectorize.args='u'),
                           lower = -Inf, upper = Inf,parameters = qn_params,ind=ind)$value)
  }


  # calculating tau_den

  qn_function1 <- function(u,parameters){min(1,qn_function_den(u,parameters))*dnorm(u)}
  tau_den <- c(integrate(Vectorize(qn_function1,vectorize.args='u'), lower = -Inf, upper = Inf,parameters = qn_params)$value)



  tau=tau_num/tau_den
  return(tau)
}







