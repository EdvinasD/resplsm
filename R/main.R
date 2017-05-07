.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Thank you for using resplsm package")
}

espls <- function(Yt, Xt,St=Yt,s,initial.values, FUN = function(...){-semi_est_func(...)},
                       bandwidth=1.06*sqrt(var(St))*length(St)^(-1/5),
                       int.of.par = c(0,1),print=F){
  # estimates parameters using newton-Raphson algorithm of semi-parametric
  # Pseudo Liklyhood.
  # input:
  #   FUN           - function for which derivative needs to be calculated
  #   Yt            - parmeter of a function which is not to be optimized,
  #                   usually Y_t
  #   Xt            - regresor parameter can be X's or lag(Y_t)
  #   initial.value - initial value of optimisible parameter might be a vector
  #   h             - limit parmerter
  #   max.iter      - number of maximum iterations until convergens
  #   delta         - stop if diferance less then it
  #   print         - self explanatory
  #   St            - at which points kernel should be applied, usually it is
  #                   Y_t but it can be done on errors too, if we estimating
  #                   ARCH processe or something similar
  #  s              - points at which kernel should be estimated


  theta.n <- length(initial.values)
  theta.s <- data.frame()
  for(si in s){
    theta0 <- initial.values


    if(print){
      print(paste("iteration: ", si))
    }
    weights <-t(matrix(dnorm((St-si)/bandwidth)))
    FUN.to.optim <- function(x){
      sum(weights*FUN(Yt,Xt,x))
    }
    if(length(initial.values)==1){
      answer <- optimize(f=FUN.to.optim,interval = int.of.par)
      theta0 <- answer$minimum
    }else{
      answer <- optim(par=initial.values,fn=FUN.to.optim)
      theta0 <- answer$par
    }


    theta.s <- rbind(theta.s,t(theta0))
  }
  return(theta.s)
}



respls <- function(theta,Y,X,
                   c_bound,iterations=5,n_par=length(theta),N=length(Y),bindwidths){
  ################################################################################
  ####                 STEP 1 initialize depending on model                   ####
  ################################################################################
  func_s=function(x,theta0){
    thetaOfX(theta0[[2]],x)^2
  }
  func_m=function(x,theta0){
    thetaOfX(theta0[[1]],x)
  }
  # inintial tau
  tau_0=matrix(0,nrow=n_par,ncol=length(Y))


  # calculate score for A matrix
  k1_0 <- k1(theta=theta,x=X,func=func_s)
  k2_0 <- k2(theta=theta,x=X,func=func_s)
  u_resids_0 <- u_resids(Y,X,theta,func_mu=func_m,func_sigma = func_s)
  score_0 <- -t(k1_0)+t(k2_0)*u_resids_0+t(k1_0)*u_resids_0^2

  # A matrixes
  A_0TA_0_inv <- t(score_0)%*%score_0/N
  A_0 <- t(solve(chol(A_0TA_0_inv)))
  A_0TA_0 <- t(A_0)%*%A_0
  A_0TA_0rep <- A_0TA_0[,rep(c(1:n_par),length(Y))]


  Thetass <- list()
  Thetass[[1]] <- theta
  for(ij in 1:iterations){
    ################################################################################
    ####                 STEP 2 Compute Tau and new A matrix                   ####
    ################################################################################

    k1_0 <- k1(theta=theta,x=X,func=func_s)
    k2_0 <- k2(theta=theta,x=X,func=func_s)
    u_resids_0 <- u_resids(Y,X,theta,func_m,func_s)
    score_0 <- -t(k1_0)+t(k2_0)*u_resids_0+t(k1_0)*u_resids_0^2

    # New Tau part -----------------------------------------------------------------
    a4=rowSums(t(k1_0)%*%A_0TA_0*t(k1_0))
    a3=2*rowSums(t(k1_0)%*%A_0TA_0*t(k2_0))
    a2=rowSums(t(k2_0)%*%A_0TA_0*t(k2_0))-2*a4-2*rowSums(t(k1_0)%*%A_0TA_0*t(tau_0))
    a1=-a3-2*rowSums(t(k2_0)%*%A_0TA_0*t(tau_0))
    a0=a4+2*rowSums(t(k1_0)%*%A_0TA_0*t(tau_0))+rowSums(t(tau_0)%*%A_0TA_0*t(tau_0))

    a00 <- a0-c_bound^2

    a_matrix <- rbind(a00,a1,a2,a3,a4)
    a=Sys.time()
    u_roots <- apply(a_matrix,2,
                     function(x){roots <- polyroot(x)
                     rotss <- sort(Re(roots)[abs(Im(roots)) < 1e-6])
                     if(length(rotss)==4){
                       warning("four roots detected, picked just 2, middle ones")
                     }
                     if(length(rotss)==0){
                       warning("0 roots detected")
                       return(rep(0,n_par))
                     }else{
                       return(rotss)
                     }
                     })
    if(!is.list(u_roots)){
      u_roots <- lapply(seq_len(ncol(u_roots)), function(i) u_roots[,i])
    }


    new_tau <- matrix(0,nrow=n_par,ncol=length(Y))

    for(i in 1:length(Y)){

      qn_params <- list(x=X[i],
                        theta0=theta,
                        k1m = k1_0[,i],
                        k2m = k2_0[,i],
                        A=A_0,
                        cb=c_bound,
                        tau=tau_0[,i],
                        func_mu=func_m,
                        func_sigma=func_s)

      u_roots_i <- u_roots[[i]]

      if(length(u_roots_i)==2){
        i_root_up <- u_roots_i[2]
        i_root_dw <- u_roots_i[1]

        new_tau[,i]=tau_for_2_roots(qn_params,i_root_dw,i_root_up)
      }
      if(length(u_roots_i)==4){
        stop('Four roots detected')
        new_tau[,i]=tau_for_4_roots(qn_params,u_roots_i)

      }
    }
    if(any(is.nan(new_tau)))warning("some of the tau vectors are not numbers")
    new_tau[is.nan(new_tau)] <- 0

    # qn_function1 <- function(u,parameters){qn_function(u,parameters)[1]*dnorm(u)}
    # qn_function2 <- function(u,parameters){qn_function(u,parameters)[2]*dnorm(u)}
    # c(integrate(Vectorize(qn_function1,vectorize.args='u'), lower = i_root_up, upper = Inf,parameters = qn_params)$value,
    #   integrate(Vectorize(qn_function2,vectorize.args='u'), lower = i_root_up, upper = Inf,parameters = qn_params)$value)

    # New A matrix part ------------------------------------------------------------

    score_tau <- score_0-t(tau_0)
    # A matrixes RECODE WITH SUM
    ATA_inv <- matrix(rowSums(apply(t(score_tau),2,function(x){
      x%*%t(x)*min(1,c_bound/norm(A_0%*%x,type="2"))^2})),ncol=n_par,nrow=n_par)/N
    new_A <- t(solve(chol(ATA_inv)))
    ATA <- t(new_A)%*%new_A
    ATArep <- ATA[,rep(c(1:n_par),length(Y))]

    params=list(func_mu=func_m,
                func_sigma=func_s,
                tau_vector=new_tau,
                A_matrix=new_A,
                Y=Y,
                X=X,
                cb=c_bound,
                bindwidths=bindwidths)

    # ToMinimizeforRobust(theta,params)
    if(ij==iterations){
      tolerance=0.01
    }else{
      tolerance=0.01*(iterations-ij)
    }


    new_theta <- theta
    for(Si in 1:length(theta[[1]]$X)){
      params$S <- theta[[1]]$X[Si]
      theta.local <- as.matrix(unlist(lapply(theta,function(x)x[Si,"value"])))
      new_theta_value=optim(par=theta.local,sp_func_to_minimize,params=params,control=list(abstol=tolerance))$par
      for(hh in 1:n_par){
        new_theta[[hh]][Si,"value"] <- new_theta_value[hh,]
      }
    }
    Thetass[[(ij+1)]]=new_theta
    theta=new_theta
    tau_0=new_tau
    A_0TA_0=ATA
    A_0=new_A
  }
  return(Thetass)
}

rls <- function(theta,Y,X,c_bound,
                                func_s,d_func_s,func_m,d_func_m,
                                iterations=5,n_par=length(theta),
                                N=length(Y)){
  ################################################################################
  ####                 STEP 1 initialize depending on model                   ####
  ################################################################################

  # inintial tau
  tau_0=matrix(0,nrow=n_par,ncol=length(Y))

  rls_k1 <- function(theta0,x,func,d_func){
    0.5/func(x,theta0)*d_func(x,theta0)
  }

  rls_k2 <- function(theta0,x,func,d_func){
    numrow <- length(theta0)
    matrix(rep(1/sqrt(func(x,theta0)),numrow),
           nrow=numrow,byrow=T)*d_func(x,theta0)
  }


  # calculate score for A matrix
  k1_0 <- rls_k1(theta=theta,x=X,func=func_s,d_func_s)
  k2_0 <- rls_k2(theta=theta,x=X,func=func_s,d_func=d_func_m)
  u_resids_0 <- u_resids(Y,X,theta,func_mu=func_m,func_sigma = func_s)
  score_0 <- -t(k1_0)+t(k2_0)*u_resids_0+t(k1_0)*u_resids_0^2

  # A matrixes
  A_0TA_0_inv <- t(score_0)%*%score_0/N
  A_0 <- t(solve(chol(A_0TA_0_inv)))
  A_0TA_0 <- t(A_0)%*%A_0
  A_0TA_0rep <- A_0TA_0[,rep(c(1:n_par),length(Y))]


  Thetass <- list()

  for(ij in 1:iterations){
    ################################################################################
    ####                 STEP 2 Compute Tau and new A matrix                   ####
    ################################################################################
    k1_0 <- rls_k1(theta=theta,x=X,func=func_s,d_func_s)
    k2_0 <- rls_k2(theta=theta,x=X,func=func_s,d_func=d_func_m)
    u_resids_0 <- u_resids(Y,X,theta,func_m,func_s)
    score_0 <- -t(k1_0)+t(k2_0)*u_resids_0+t(k1_0)*u_resids_0^2

    # New Tau part -----------------------------------------------------------------
    a4=rowSums(t(k1_0)%*%A_0TA_0*t(k1_0))
    a3=2*rowSums(t(k1_0)%*%A_0TA_0*t(k2_0))
    a2=rowSums(t(k2_0)%*%A_0TA_0*t(k2_0))-2*a4-2*rowSums(t(k1_0)%*%A_0TA_0*t(tau_0))
    a1=-a3-2*rowSums(t(k2_0)%*%A_0TA_0*t(tau_0))
    a0=a4+2*rowSums(t(k1_0)%*%A_0TA_0*t(tau_0))+rowSums(t(tau_0)%*%A_0TA_0*t(tau_0))

    a00 <- a0-c_bound^2

    a_matrix <- rbind(a00,a1,a2,a3,a4)

    u_roots <- apply(a_matrix,2,
                     function(x){roots <- polyroot(x)
                     rotss <- sort(Re(roots)[abs(Im(roots)) < 1e-6])

                     if(length(rotss)==4){
                       warning(paste("it:",ij,"four roots detected, careful"))
                     }
                     if(length(rotss)==0){
                       warning("0 roots detected")
                       return(rep(0,n_par))
                     }else{
                       return(rotss)
                     }
                     })
    if(!is.list(u_roots)){
      u_roots <- lapply(seq_len(ncol(u_roots)), function(i) u_roots[,i])
    }


    new_tau <- matrix(0,nrow=n_par,ncol=length(Y))

    for(i in 1:length(Y)){

      qn_params <- list(x=X[i],
                        theta0=theta,
                        k1m = k1_0[,i],
                        k2m = k2_0[,i],
                        A=A_0,
                        cb=c_bound,
                        tau=tau_0[,i],
                        func_mu=func_m,
                        func_sigma=func_s)

      u_roots_i <- u_roots[[i]]

      if(length(u_roots_i)==2){
        i_root_up <- u_roots_i[2]
        i_root_dw <- u_roots_i[1]

        new_tau[,i]=tau_for_2_roots(qn_params,i_root_dw,i_root_up)
      }
      if(length(u_roots_i)==4){
        # used integrate function in R because wierd results from approximation
        new_tau[,i]=tau_for_4_roots(qn_params,u_roots_i)

      }
    }
    if(any(is.nan(new_tau)))warning("some of the tau vectors are not numbers")
    new_tau[is.nan(new_tau)] <- 0


    # New A matrix part ------------------------------------------------------------

    k1_0 <- rls_k1(theta=theta,x=X,func=func_s,d_func_s)
    k2_0 <- rls_k2(theta=theta,x=X,func=func_s,d_func=d_func_m)
    u_resids_0 <- u_resids(Y,X,theta,func_mu=func_m,func_sigma = func_s)
    score_0 <- -t(k1_0)+t(k2_0)*u_resids_0+t(k1_0)*u_resids_0^2

    score_tau <- score_0-t(tau_0)
    # A matrixes RECODE WITH SUM
    ATA_inv <- matrix(rowSums(apply(t(score_tau),2,function(x){
      x%*%t(x)*min(1,c_bound/norm(A_0%*%x,type="2"))^2})),ncol=n_par,nrow=n_par)/N
    new_A <- t(solve(chol(ATA_inv)))
    ATA <- t(new_A)%*%new_A
    # ATArep <- ATA[,rep(c(1:n_par),length(Y))]

    params=list(func_mu=func_m,
                func_sigma=func_s,
                d_func_mu=d_func_m,
                d_func_sigma=d_func_s,
                tau_vector=new_tau,
                A_matrix=new_A,
                Y=Y,
                X=X,
                cb=c_bound,
                k1=k1,
                k2=k2)


    new_theta=optim(par=theta,func_to_minimize,params=params)$par
    Thetass[[ij]]=new_theta
    theta=new_theta
    tau_0=new_tau
    A_0TA_0=ATA
    A_0=new_A

  }
  return(Thetass)
}
