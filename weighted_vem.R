#########Variational EM by sample size weights;###############
#############################################################
ELBO_alpha_weighted <- function(alpha,a1,a2,S,B,pi,Xs,WB,ytrt,ysc,Ub=NULL,sigmabj2,sigma2,mubj,theta,beta,w,Ls,Lb,K,lambda = 0.1){
  #subject_id = unique(cellidx)
  alpha0 = alpha[1:K]
  alpha1 = alpha[(K+1):(2*K)]
  alpha2 = alpha[-(1:(2*K))]
  term1 = -Xs%*%alpha-w%*%mubj+w%*%sigmabj2/2-log(Ls)
  term2 = Xs%*%alpha+w%*%mubj+w%*%sigmabj2/2+log(Ls)
  
  Xb = model.matrix(~ytrt)
  alpha.matrix <- rbind(alpha0,alpha1)
  term.bulk = log(rowSums(WB*exp(Xb%*%alpha.matrix)))+log(Lb)+ifelse(!is.null(Ub),Ub%*%alpha2,0)
  
  
  tmp2 =  a2*(-(B+beta)*log(1+1/beta*exp(term.bulk))+ B*term.bulk)
  tmp1 = a1*(-(1-pi)*(S>0)*((S+theta)*log(1+theta*exp(term1))+theta*Xs%*%alpha)-
               (1-pi)*(S==0)*theta*log(1+1/theta*exp(term2)))
  Falpha <- sum(tmp1) + sum(tmp2) - a2*lambda*sum(WB%*%alpha1^2)
  return(-Falpha)
}
grad_alpha_weighted <-function(alpha,a1,a2,S,B,pi,Xs,WB,ytrt,ysc,Ub=NULL,sigmabj2,sigma2,mubj,theta,beta,w,Ls,Lb,K,lambda = 0.1){
  #subject_id = unique(cellidx)
  alpha0 = alpha[1:K]
  alpha1 = alpha[(K+1):(2*K)]
  alpha2 = alpha[-(1:(2*K))]
  Ws = Xs[,1:K]
  Wstrt = Xs[,(K+1):(2*K)]
  Us = Xs[,-(1:(2*K))]
  
  term1 = -Xs%*%alpha-w%*%mubj+w%*%sigmabj2/2-log(Ls)
  term2 = Xs%*%alpha+w%*%mubj+w%*%sigmabj2/2+log(Ls)
  
  
  tmp10 = as.vector((1-pi)*(S>0)*(S+theta)*theta*exp(term1))*Ws/c(1+theta*exp(term1))-as.vector((1-pi)*(S>0)*theta)*Ws-
    as.vector((1-pi)*(S==0)*exp(term2))*Ws/c(1+1/theta*exp(term2))
  
  tmp11 = as.vector((1-pi)*(S>0)*(S+theta)*theta*exp(term1))*Wstrt/c(1+theta*exp(term1))-as.vector((1-pi)*(S>0)*theta)*Wstrt-
    as.vector((1-pi)*(S==0)*exp(term2))*Wstrt/c(1+1/theta*exp(term2))
  
  Xb = model.matrix(~ytrt)
  alpha.matrix <- rbind(alpha0,alpha1)
  covariate_term = ifelse(!is.null(Ub),Ub%*%alpha2,0)
  term.bulk = log(rowSums(WB*exp(Xb%*%alpha.matrix)))+log(Lb)+ covariate_term
  term.bulk.gr1 = WB*exp(Xb%*%alpha.matrix)*Lb*exp(covariate_term)
  
  tmp20 = -(B+beta)*1/beta*term.bulk.gr1/c(1+1/beta*exp(term.bulk))+
    B*WB*exp(Xb%*%alpha.matrix)/rowSums(WB*exp(Xb%*%alpha.matrix))
  tmp21 =  -(B+beta)*1/beta*term.bulk.gr1*ytrt/c(1+1/beta*exp(term.bulk))+
    B*ytrt*WB*exp(Xb%*%alpha.matrix)/rowSums(WB*exp(Xb%*%alpha.matrix))
  
  
  if(is.null(Ub)){
    gradalpha <- c(a1*colSums(tmp10) + a2*colSums(tmp20),a1*colSums(tmp11)+a2*colSums(tmp21)-2*a2*lambda*colSums(WB)*alpha1) ####to be verified the dimension
    
  }else{
    tmp22 = -(B+beta)/beta*exp(term.bulk)*Ub/(1+1/beta*exp(term.bulk))+B*Ub
    tmp12 = as.vector((1-pi)*(S>0)*(S+theta)*theta*exp(term1))*Us/c(1+theta*exp(term1))-as.vector((1-pi)*(S>0)*theta)*Us-
      as.vector((1-pi)*(S==0)*exp(term2))*Us/c(1+1/theta*exp(term2))
    gradalpha <- c(a1*colSums(tmp10) + a2*colSums(tmp20),a1*colSums(tmp11)+a2*colSums(tmp21)-2*a2*lambda*colSums(WB)*alpha1, a1*colSums(tmp12)+a2*colSums(tmp22)) ####to be verified the dimension
    
  }
  
  return(-gradalpha)
}
# ELBO_total_weighted <- function(S,B,a1,a2,pi,Xs,Xsz,gamma,WB,ytrt,ysc,Ub=NULL,sigmabj2,sigma2,mubj,alpha,theta,beta,mum,sm2,w,Ls,Lb,K){
#   alpha0 = alpha[1:K]
#   alpha1 = alpha[(K+1):(2*K)]
#   alpha2 = alpha[-(1:(2*K))]
#   #####For zero inflation part#############
#   ##########################################
#   
#   Fspart = 0.5*sum((log(sigmabj2)-log(sigma2)-(mubj^2+sigmabj2)/sigma2)) -(log(sigma2)-mum)^2/2/sm2
#   
#   
#   Fzero = sum(pi*Xsz%*%gamma - log(1+exp(Xsz%*%gamma)))
#   
#   
#   Fpi = -sum(pi[pi>0]*log(pi[pi>0])) -sum((1-pi[pi<1])*log(1-pi[pi<1]))
#   
#   Xb = model.matrix(~ytrt)
#   alpha.matrix <- rbind(alpha0,alpha1)
#   covariate_term = ifelse(!is.null(Ub),Ub%*%alpha2,0)
#   term.bulk = log(rowSums(WB*exp(Xb%*%alpha.matrix)))+log(Lb)+ covariate_term
#   
#   term1 = -Xs%*%alpha-w%*%mubj+w%*%sigmabj2/2-log(Ls)
#   term2 = Xs%*%alpha+w%*%mubj+w%*%sigmabj2/2+log(Ls)
#   
#   tmp2 = lgamma(B+beta)-lgamma(beta)-lgamma(B+1)-(B+beta)*log(1+1/beta*exp(term.bulk))+B*term.bulk - B*log(beta)
#   
#   
#   
#   lgamma_diff <- lgamma(S+theta)-lgamma(theta)-lgamma(S+1)
#   
#   
#   logmu =  Xs%*%alpha+w%*%mubj+log(Ls)
#   common_term1 <- theta*exp(term1)
#   common_term2 <- 1/theta*exp(term2)
#   
#   tmp1 = (1-pi)*(S>0)*(lgamma_diff - (S+theta)*log(1+common_term1)+theta*log(theta)-theta*logmu)-(1-pi)*(S==0)*theta*log(1+common_term2)
#   
#   Ftotal <- a2*sum(tmp2) + a1*(sum(tmp1) + Fspart + Fzero + Fpi)
#   return(-Ftotal)
# }
ELBO_total_weighted <- function(S,B,a1,a2,pi,Xs,Xsz,gamma,WB,ytrt,ysc,Ub=NULL,sigmabj2,sigma2,mubj,alpha,theta,beta,w,Ls,Lb,K,lambda = 0.1){
  alpha0 = alpha[1:K]
  alpha1 = alpha[(K+1):(2*K)]
  alpha2 = alpha[-(1:(2*K))]
  #####For zero inflation part#############
  ##########################################
  
  Fspart = 0.5*sum((log(sigmabj2)-log(sigma2)-(mubj^2+sigmabj2)/sigma2)) 
  
  
  Fzero = sum(pi*Xsz%*%gamma - log(1+exp(Xsz%*%gamma)))
  
  
  Fpi = -sum(pi[pi>0]*log(pi[pi>0])) -sum((1-pi[pi<1])*log(1-pi[pi<1]))
  
  Xb = model.matrix(~ytrt)
  alpha.matrix <- rbind(alpha0,alpha1)
  covariate_term = ifelse(!is.null(Ub),Ub%*%alpha2,0)
  term.bulk = log(rowSums(WB*exp(Xb%*%alpha.matrix)))+log(Lb)+ covariate_term
  
  term1 = -Xs%*%alpha-w%*%mubj+w%*%sigmabj2/2-log(Ls)
  term2 = Xs%*%alpha+w%*%mubj+w%*%sigmabj2/2+log(Ls)
  
  tmp2 = lgamma(B+beta)-lgamma(beta)-lgamma(B+1)-
    (B+beta)*log(1+1/beta*exp(term.bulk))+B*term.bulk - B*log(beta) 
  
  
  
  lgamma_diff <- lgamma(S+theta)-lgamma(theta)-lgamma(S+1)
  
  
  logmu =  Xs%*%alpha+w%*%mubj+log(Ls)
  common_term1 <- theta*exp(term1)
  common_term2 <- 1/theta*exp(term2)
  
  tmp1 = (1-pi)*(S>0)*(lgamma_diff - (S+theta)*log(1+common_term1)+theta*log(theta)-theta*logmu)-(1-pi)*(S==0)*theta*log(1+common_term2)
  
  Ftotal <- a2*sum(tmp2) + a1*(sum(tmp1) + Fspart + Fzero + Fpi) - a2*lambda*sum(WB%*%alpha1^2)
  return(-Ftotal)
}

variational_weightedEM <- function(S,B,a1,a2,prior_alpha,param_initial = NULL,Xs,WB,Ub=NULL,ysc,
                           ytrt,cellidx,Xsz,Lb,mum,sm2,
                           Ls,K,eps = 1e-9,max.iter = 200, trace = 0,lambda = 0.1){
  ns = length(S)
  ord_id=table(cellidx)
  nsub = length(ord_id)
  w_l=list()
  for(i in 1:length(ord_id)){
    w_l[[i]]=matrix(1,ord_id[i],1)
  }
  w=as.matrix(Matrix::bdiag(w_l))
  
  if(!is.null(param_initial)){
    if("alpha" %in% names(param_initial)){
      alpha_initial = param_initial$alpha
    }else{
      alpha_initial = prior_alpha
    }
    
    if("beta" %in% names(param_initial)){
      beta_initial = param_initial$beta
    }else{
      beta_initial = 0.1
    }
    if("theta" %in% names(param_initial)){
      theta_initial = param_initial$theta
    }else{
      theta_initial = 0.1
    }
    if("gamma" %in% names(param_initial)){
      gamma_initial = param_initial$gamma
    }else{
      gamma_initial = coef(glm(S==0~ysc,family="binomial"))
    }
    if("sigma2" %in% names(param_initial)){
      sigma2_initial = param_initial$sigma2
    }else{
      sigma2_initial = 0.05
    }
  }else{
    
    gamma_initial = coef(glm(S==0~ysc,family="binomial"))
    alpha_initial = prior_alpha
    beta_initial = 0.1
    theta_initial = 0.1
    sigma2_initial = 0.05
  }
  
  ncov = ifelse(is.null(Ub),0,ncol(Ub))
  sigmabj_initial = sigmabj2 = rep(0.01,nsub)
  pi_initial = rep(mean(S==0),ns)
  mubj_initial = rnorm(nsub,0,0.1)
  #sigma2_initial = exp(rnorm(1,mum,sqrt(sm2)))
  
  gamma_current = gamma_initial
  alpha_current = alpha_initial
  beta_current = beta_initial
  theta_current = theta_initial
  mubj_current = mubj_initial
  sigmabj_current = sigmabj_initial
  sigma2_current = sigma2_initial
  pi_current = pi_initial
  
  
  current.loglik = -1e6
  loglik = c()
  iter = 1
  div = 1e6
  
  llgamma = llmubj = llsigmabj = lltheta = llbeta = llalpha = llsigma2 = -1e6
  

  while((div> eps*(abs(current.loglik)+eps)) && iter <= max.iter) {
    
    
    print(iter)
    # eps.gam = 1e-3
    # delta.gamma = gamma_current; p.iter = 1
    # while(!all(delta.gamma < eps.gam) & (p.iter < p.max.iter)){
    #   print(p.iter)
    q=try(nlminb(start =as.vector(gamma_current), ELBO_gamma,
                 gradient = grad_gamma,lower = rep(-Inf,ncol(Xsz)), upper = rep(Inf,ncol(Xsz)),
                 pi=pi_current,Xsz=Xsz,control = list(trace = 1, iter.max = 100)), silent = TRUE)
    if(!inherits(q, "try-error")) {
      if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
      if(llgamma > -q$objective) { if(trace) cat("Optimization gamma did not improve estimates of variational parameters.\n"); gamma_new <- gamma_current }
      else {
        if(trace) cat("Variational gamma parameters updated", "\n")
        gamma_new <- q$par
      }
    }else { gamma_new <- gamma_current}
    llgamma <- -q$objective
    #print(llgamma)
    
    pi_new = as.vector(ELBO_pi_update(alpha = alpha_current,S = S,
                                      mubj = mubj_current,sigmabj2 = sigmabj_current,
                                      gamma = gamma_new,theta = theta_current,w = w,
                                      Xsz=Xsz,Xs=Xs,Ls =Ls))
    if(any(pi_new[S>0]>0)) stop("see what happened")
    
    # pi_new = as.vector(ELBO_pi0(alpha = alpha_current,
    #                                   mubj = mubj_current,sigmabj2 = sigmabj_current,
    #                                   gamma = gamma_new,theta = theta_current,w = w,
    #                                   Xsz=Xsz,Xs=Xs,Ls =Ls))
    # pi_new = sapply(1:ns,function(i){
    #   expo_term <- 1/theta_current*Ls[i]*exp(Xs[i,]%*%alpha_current + mubj_current[cellidx[i]]+sigmabj_current[cellidx[i]]/2)
    #   denom_term <- exp(-Xsz[i,]%*%gamma_new - theta_current*log(1+expo_term))
    #   (S[i] == 0)*1/(1+denom_term)
    # })
    
    #   delta.gamma <- abs(gamma_current - gamma_new)
    #   p.iter <- p.iter+1
    #   pi_current = pi_new
    #   gamma_current = gamma_new
    # }
    
    # delta.sigmafij.required <- 1e-3; p.iter <- 1; p.max.iter <- 4; delta.sigma2<-sigma2_current ;
    # 
    # q0 = optim(par = mubj_current, fn = ELBO_mubj_total, lower = rep(-Inf,nsub), upper = rep(Inf,nsub),
    #        gr =grad_mubj_total,S=S,pi=pi_current, method = "BFGS",
    #        Xs=Xs,sigmabj2 =sigmabj_current,
    #        sigma2 = sigma2_current,alpha = alpha_current,
    #        theta = theta_current,cellidx = cellidx,
    #        Ls = libsizes )
    
    #while(!all(delta.sigma2 < delta.sigmafij.required) && p.iter < p.max.iter){
    # delta.sigmabj = 1e3; eps.bj = 1e-3; p.iter=1
    # while(!all(delta.sigmabj < eps.bj) & (p.iter < p.max.iter)){
    #   
    # }
    q=try(nlminb(start =as.vector(mubj_current), ELBO_mubj , gradient = grad_mubj,lower = rep(-Inf,nsub), upper = rep(Inf,nsub),
                 S=S,pi=pi_new,
                 Xs=Xs,sigmabj2 =sigmabj_current,
                 sigma2 = sigma2_current,alpha = alpha_current,
                 theta = theta_current,w=w,
                 Ls = Ls,control = list(trace = 0, iter.max = 100)), silent = TRUE)
    if(!inherits(q, "try-error")) {
      if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
      if( llmubj > -q$objective) { if(trace) cat("Optimization mubj did not improve estimates of variational parameters.\n"); mubj_new <- mubj_current}
      else {
        if(trace) cat("Variational mubj parameters updated", "\n")
        mubj_new <- q$par
      }
    }else{
      mubj_new<- mubj_current
    }
    llmubj <- -q$objective
    
    #print(llmubj)
    
    
    q=try(nlminb(start =as.vector(sigmabj_current), ELBO_sigmabj , gradient = grad_sigmabj,lower = rep(0,nsub), upper = rep(Inf,nsub),
                 S=S,pi=pi_new,
                 Xs=Xs,mubj =mubj_new,
                 sigma2 = sigma2_current,alpha = alpha_current,
                 theta = theta_current,w = w,
                 Ls = Ls,control = list(trace = 0, iter.max = 100)), silent = TRUE)
    if(!inherits(q, "try-error")) {
      if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
      if(llsigmabj > -q$objective) { if(trace) cat("Optimization sigbj did not improve estimates of variational parameters.\n"); sigmabj_new <- sigmabj_current}
      else {
        if(trace) cat("Variational sigbj parameters updated", "\n")
        sigmabj_new <- q$par
      }
    }else{
      sigmabj_new<- sigmabj_current
    }
    llsigmabj <- -q$objective
    
    #print(llsigmabj)
    
    # q=nlminb(start =as.vector(sigma2_current), ELBO_sigma2, gradient = grad_sigma2,lower = rep(0,nsub), upper = rep(Inf,nsub),
    #              mubj =mubj_new,
    #              sigmabj2 = sigmabj_new,mum = mum,
    #              sm2=sm2,control = list(trace = 0, iter.max = 100))
    # if(!inherits(q, "try-error")) {
    #   if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
    #   if(llsigma2 > -q$objective) { if(trace) cat("Optimization sigma2 did not improve estimates of variational parameters.\n"); sigma2_new <- sigma2_current}
    #   else {
    #     if(trace) cat("Variational sigma2 parameters updated", "\n")
    #     sigma2_new <- q$par
    #   }
    # }else{
    #   sigma2_new<- sigma2_current
    # }
    # llsigma2 <- -q$objective
    
    sigma2_new=sum(sigmabj_new+mubj_new^2)/nsub
    
  
    
    
    
    q=nlminb(start =as.vector(alpha_current), ELBO_alpha_weighted, 
             gradient = grad_alpha_weighted,
                 lower = rep(-Inf,2*K), upper = rep(Inf,2*K),
                 S=S,B=B,a1=a1,a2=a2,pi=pi_new,Xs=Xs,WB = WB, ytrt = ytrt,ysc = ysc,
                 Ub=NULL,mubj =mubj_new,sigmabj2 = sigmabj_new,sigma2 = sigma2_new,
                 theta = theta_current,beta = beta_current,w = w,
                 Ls = Ls,Lb = Lb,K=K,lambda = lambda,control = list(trace = 0, iter.max = 200))
    
    
    
    if(!inherits(q, "try-error")) {
      if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
      if(llalpha > -q$objective) { if(trace) cat("Optimization alpha did not improve estimates of variational parameters.\n"); alpha_new <- alpha_current }
      else {
        if(trace) cat("Variational parameters alpha updated", "\n")
        alpha_new <- q$par
      }
    }else { alpha_new <- alpha_current}
    llalpha <- -q$objective
    #print(llalpha)
    
    q=try(nlminb(start= theta_current, ELBO_theta, gradient = grad_theta,lower = 0, upper = Inf,
                 S=S,pi=pi_new,alpha = alpha_new,
                 Xs=Xs,mubj =mubj_new,sigmabj2 = sigmabj_new,sigma2 = sigma2_new,
                 w = w,Ls = Ls,K=K,control = list(trace = 0, iter.max = 100)), silent = TRUE)
    
    
    if(!inherits(q, "try-error")) {
      if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
      if( lltheta > -q$objective) { if(trace) cat("Optimization theta did not improve estimates of variational parameters.\n"); theta_new <- theta_current }
      else {
        if(trace) cat("Variational parameters theta updated", "\n")
        theta_new <- q$par
      }
    }else {theta_new <- theta_current}
    
    lltheta<- -q$objective
    
    
    q=try(nlminb(start =as.vector(beta_current), ELBO_beta, gradient = grad_beta,
                lower = 0, upper = Inf,
                 alpha = alpha_new,B=B,WB=WB,ytrt=ytrt,Ub=NULL,
                 Lb=Lb,K=K,control = list(trace = 0, iter.max = 100)), silent = TRUE)
    
    
    if(!inherits(q, "try-error")) {
      if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
      if(llbeta > -q$objective) { if(trace) cat("Optimization beta did not improve estimates of variational parameters.\n"); beta_new <- beta_current }
      else {
        if(trace) cat("Variational parameters beta updated", "\n")
        beta_new <- q$par
      }
    }else {beta_new <- beta_current}
    
    llbeta<- -q$objective
    
    
   
    objective =ELBO_total_weighted(S=S,B=B,a1=a1,a2=a2,pi=pi_new,
                          Xs=Xs,Xsz=Xsz,gamma=gamma_new,
                          WB = WB,ytrt = ytrt,ysc = ysc,alpha = alpha_new,
                          Ub=NULL,sigmabj2 = sigmabj_new,sigma2=sigma2_new,
                          mubj=mubj_new,theta=theta_new,beta = beta_new,
                          w = w,Ls = Ls,Lb = Lb,K = K,lambda = lambda)
    new.loglik <- objective
    div=abs(new.loglik-current.loglik)
    #err <- abs(new.loglik/current.loglik);
    
    current.loglik <- new.loglik
    
    iter <- iter + 1
    loglik=append(loglik,current.loglik)
    
    print(loglik)
    alpha_current = alpha_new
    beta_current = beta_new
    theta_current = theta_new
    mubj_current = mubj_new
    sigmabj_current = sigmabj_new
    sigma2_current = sigma2_new
    pi_current = pi_new
    gamma_current = gamma_new
    
  }
  return(list(alpha = alpha_new,beta = beta_new, 
              theta = theta_new, sigma2 = sigma2_new, 
              gamma = gamma_new, pi = pi_new, mubj = mubj_new,
              sigmabj2 = sigmabj_new,loglik = current.loglik))
}

variational_weightedEM_matrix <- function(Y,B0,a1,a2,prior_alpha,param_initial = NULL,Xs,WB,Ub=NULL,ysc,
                                   ytrt,cellidx,Xsz,Lb,
                                   Ls,K,eps = 1e-9,max.iter = 200, trace = 0){
  ns = ncol(Y)
  ng = nrow(Y)
  nb = ncol(B0)
  ord_id=table(cellidx)
  nsub = length(ord_id)
  w_l=list()
  for(i in 1:length(ord_id)){
    w_l[[i]]=matrix(1,ord_id[i],1)
  }
  w=as.matrix(Matrix::bdiag(w_l))
  
  if(!is.null(param_initial)){
    if("alpha" %in% names(param_initial)){
      alpha_initial = param_initial$alpha
    }else{
      alpha_initial = prior_alpha
    }
    
    if("beta" %in% names(param_initial)){
      beta_initial = param_initial$beta
    }else{
      beta_initial = 0.1
    }
    if("theta" %in% names(param_initial)){
      theta_initial = param_initial$theta
    }else{
      theta_initial = 0.1
    }
    if("gamma" %in% names(param_initial)){
      gamma_initial = param_initial$gamma
    }else{
      gamma_initial = t(apply(Y,1,function(x){coef(glm(x==0~ysc,family="binomial"))}))
    }
    if("sigma2" %in% names(param_initial)){
      sigma2_initial = param_initial$sigma2
    }else{
      sigma2_initial = 0.05
    }
  }else{
    
    gamma_initial = t(apply(Y,1,function(x){coef(glm(x==0~ysc,family="binomial"))}))
    alpha_initial = prior_alpha
    beta_initial = 0.1
    theta_initial = 0.1
    sigma2_initial = 0.05
  }
  
  ncov = ifelse(is.null(Ub),0,ncol(Ub))
  sigmabj_initial = sigmabj2 = rep(0.01,nsub)
  #pi_initial = matrix(rep(rowMeans(Y == 0),ns),nrow = ng, ncol = ns,byrow = FALSE)
  mubj_initial = rnorm(nsub,0,0.1)
  
  sc_predict = matrix(0,nrow=ng,ncol = ns)
  bulk_predict = matrix(0,nrow=ng,ncol = nb)
  gamma_all = gamma_initial
  alpha_all = prior_alpha
  beta_all = rep(beta_initial,ng)
  theta_all = rep(theta_initial,ng)
  mubj_all =  matrix(rep(mubj_initial,ng),nrow = ng, ncol = nsub,byrow=TRUE)
  sigmabj_all =matrix(rep(sigmabj_initial,ng),nrow = ng, ncol = nsub,byrow=TRUE)
  sigma2_all = rep(sigma2_initial,ng)
  pi_all = matrix(rep(rowMeans(Y == 0),ns),nrow = ng, ncol = ns,byrow = FALSE)
  current.loglik = -1e6
  loglik = c()
  iter = 1
  div = 1e6
  lla = llb = llth = lls = llm = llg = rep(-1e6,ng)
    while((div> eps*(abs(current.loglik)+eps)) && iter <= max.iter) {
      new.loglik = 0
      
      for(g in 1:nrow(B0)){
        gamma_current = gamma_all[g,]
        alpha_current = alpha_all[g,]
        beta_current = beta_all[g]
        theta_current = theta_all[g]
        mubj_current = mubj_all[g,]
        sigmabj_current = sigmabj_all[g,]
        sigma2_current = sigma2_all[g]
        pi_current = pi_all[g,]
        
        
        
        S = Y[g,]
        B = B0[g,]
        llgamma = llg[g]
        llmubj = llm[g]
        llsigmabj = lls[g]
          lltheta = llth[g]
            llbeta = llb[g]
              llalpha = lla[g]
                
        
      if(!is.null(param_initial)){
        if("alpha" %in% names(param_initial)){
          alpha_initial = param_initial$alpha[g,]
        }else{
          alpha_initial = prior_alpha[g,]
        }
      }else{
        alpha_initial = prior_alpha[g,]
      }
      
      #print(iter)
      # eps.gam = 1e-3
      # delta.gamma = gamma_current; p.iter = 1
      # while(!all(delta.gamma < eps.gam) & (p.iter < p.max.iter)){
      #   print(p.iter)
      q=try(nlminb(start =as.vector(gamma_current), ELBO_gamma,
                   gradient = grad_gamma,lower = rep(-Inf,ncol(Xsz)), upper = rep(Inf,ncol(Xsz)),
                   pi=pi_current,Xsz=Xsz,control = list(trace =0, iter.max = 100)), silent = TRUE)
      if(!inherits(q, "try-error")) {
        if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
        if(llgamma > -q$objective) { if(trace) cat("Optimization gamma did not improve estimates of variational parameters.\n"); gamma_new <- gamma_current }
        else {
          if(trace) cat("Variational gamma parameters updated", "\n")
          gamma_new <- q$par
        }
      }else { gamma_new <- gamma_current}
      llgamma <- -q$objective
      #print(llgamma)
      
      pi_new = as.vector(ELBO_pi_update(alpha = alpha_current,S = S,
                                        mubj = mubj_current,sigmabj2 = sigmabj_current,
                                        gamma = gamma_new,theta = theta_current,w = w,
                                        Xsz=Xsz,Xs=Xs,Ls =Ls))
      if(any(pi_new[S>0]>0)) stop("see what happened")
      
     
      q=try(nlminb(start =as.vector(mubj_current), ELBO_mubj , gradient = grad_mubj,lower = rep(-Inf,nsub), upper = rep(Inf,nsub),
                   S=S,pi=pi_new,
                   Xs=Xs,sigmabj2 =sigmabj_current,
                   sigma2 = sigma2_current,alpha = alpha_current,
                   theta = theta_current,w=w,
                   Ls = Ls,control = list(trace = 0, iter.max = 100)), silent = TRUE)
      if(!inherits(q, "try-error")) {
        if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
        if( llmubj > -q$objective) { if(trace) cat("Optimization mubj did not improve estimates of variational parameters.\n"); mubj_new <- mubj_current}
        else {
          if(trace) cat("Variational mubj parameters updated", "\n")
          mubj_new <- q$par
        }
      }else{
        mubj_new<- mubj_current
      }
      llmubj <- -q$objective
      
      #print(llmubj)
      
      
      q=try(nlminb(start =as.vector(sigmabj_current), ELBO_sigmabj , gradient = grad_sigmabj,lower = rep(0,nsub), upper = rep(Inf,nsub),
                   S=S,pi=pi_new,
                   Xs=Xs,mubj =mubj_new,
                   sigma2 = sigma2_current,alpha = alpha_current,
                   theta = theta_current,w = w,
                   Ls = Ls,control = list(trace = 0, iter.max = 100)), silent = TRUE)
      if(!inherits(q, "try-error")) {
        if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
        if(llsigmabj > -q$objective) { if(trace) cat("Optimization sigbj did not improve estimates of variational parameters.\n"); sigmabj_new <- sigmabj_current}
        else {
          if(trace) cat("Variational sigbj parameters updated", "\n")
          sigmabj_new <- q$par
        }
      }else{
        sigmabj_new<- sigmabj_current
      }
      llsigmabj <- -q$objective
      
      #print(llsigmabj)
      
      
      sigma2_new=sum(sigmabj_new+mubj_new^2)/nsub
      
      #   delta.sigma2 <- abs(sigma2_new-sigma2_current)
      #   sigma2_current <- sigma2_new
      #   p.iter <- p.iter+1
      # }    
      
      
      
      
      q=try(nlminb(start =as.vector(alpha_current), ELBO_alpha_weighted, gradient = grad_alpha_weighted,
                   lower = rep(-Inf,2*K), upper = rep(Inf,2*K),
                   S=S,B=B,a1=a1,a2=a2,pi=pi_new,Xs=Xs,WB = WB, ytrt = ytrt,ysc = ysc,
                   Ub=NULL,mubj =mubj_new,sigmabj2 = sigmabj_new,sigma2 = sigma2_new,
                   theta = theta_current,beta = beta_current,w = w,
                   Ls = Ls,Lb = Lb,K=K,control = list(trace = 0, iter.max = 200)), silent = TRUE)
      
      
      
      if(!inherits(q, "try-error")) {
        if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
        if(llalpha > -q$objective) { if(trace) cat("Optimization alpha did not improve estimates of variational parameters.\n"); alpha_new <- alpha_current }
        else {
          if(trace) cat("Variational parameters alpha updated", "\n")
          alpha_new <- q$par
        }
      }else { alpha_new <- alpha_current}
      llalpha <- -q$objective
      #print(llalpha)
      
      q=try(nlminb(start= theta_current, ELBO_theta , gradient = grad_theta,lower = 0, upper = Inf,
                   S=S,pi=pi_new,alpha = alpha_new,
                   Xs=Xs,mubj =mubj_new,sigmabj2 = sigmabj_new,sigma2 = sigma2_new,
                   w = w,Ls = Ls,K=K,control = list(trace = 0, iter.max = 100)), silent = TRUE)
      
      
      if(!inherits(q, "try-error")) {
        if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
        if( lltheta > -q$objective) { if(trace) cat("Optimization theta did not improve estimates of variational parameters.\n"); theta_new <- theta_current }
        else {
          if(trace) cat("Variational parameters theta updated", "\n")
          theta_new <- q$par
        }
      }else {theta_new <- theta_current}
      
      lltheta<- -q$objective
      
      
      q=try(nlminb(start =as.vector(beta_current), ELBO_beta, gradient = grad_beta,lower = 0, upper = Inf,
                   alpha = alpha_new,B=B,WB=WB,ytrt=ytrt,Ub=NULL,
                   Lb=Lb,K=K,control = list(trace = 0, iter.max = 100)), silent = TRUE)
      
      
      if(!inherits(q, "try-error")) {
        if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
        if(llbeta > -q$objective) { if(trace) cat("Optimization beta did not improve estimates of variational parameters.\n"); beta_new <- beta_current }
        else {
          if(trace) cat("Variational parameters beta updated", "\n")
          beta_new <- q$par
        }
      }else {beta_new <- beta_current}
      
      llbeta<- -q$objective
     
      alpha.matrix = rbind(alpha_new[1:(K)],alpha_new[(K+1):(2*K)])
      Xb = model.matrix(~ytrt)
      mus_predict <- exp(Xs%*%alpha_new + mubj_new)*(1-pi_new)
      mub_predict <- rowSums(WB*exp(Xb%*%alpha.matrix))
      
      sc_predict[g,] = mus_predict
      bulk_predict[g,] = mub_predict
      objective =ELBO_total_weighted(S=S,B=B,a1=a1,a2=a2,pi=pi_new,
                                     Xs=Xs,Xsz=Xsz,gamma=gamma_new,
                                     WB = WB,ytrt = ytrt,ysc = ysc,alpha = alpha_new,
                                     Ub=NULL,sigmabj2 = sigmabj_new,sigma2=sigma2_new,
                                     mubj=mubj_new,theta=theta_new,beta = beta_new,
                                     w = w,Ls = Ls,Lb = Lb,K = K)
      new.loglik <- new.loglik + objective
      
      
      gamma_all[g,] = gamma_new
      alpha_all[g,] = alpha_new
      beta_all[g] = beta_new
      theta_all[g] = theta_new
      mubj_all[g,] =  mubj_new
      sigmabj_all[g,] =sigmabj_new
      sigma2_all[g] = sigma2_new
      pi_all[g,] = pi_new
      
      llg[g]=llgamma 
      llm[g]=llmubj 
      lls[g]=llsigmabj 
      llth[g]=lltheta
      llb[g]=llbeta
      lla[g]=llalpha
      }
      
      #err <- abs(new.loglik/current.loglik);
      
      q=try(nlminb(start =as.vector(Ls), ELBO_Ls_total, gradient = grad_Ls_total,
                   lower = 0, upper = Inf,
                   alpha_all = alpha_all,theta_all=theta_all,pi_all=pi_all,
                   ytrt=ytrt,Ub=NULL,
                   Lb=Lb,K=K,control = list(trace = 0, iter.max = 100)), silent = TRUE)
      
      
      if(!inherits(q, "try-error")) {
        if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
        if(llbeta > -q$objective) { if(trace) cat("Optimization beta did not improve estimates of variational parameters.\n"); beta_new <- beta_current }
        else {
          if(trace) cat("Variational parameters beta updated", "\n")
          beta_new <- q$par
        }
      }else {beta_new <- beta_current}
      
       #<- function(Ls,alpha_all,theta_all,pi_all,sigmbj_all,sigma2_all,mubj_all,w,K)
      Ls = colMeans(Y)/colMeans(sc_predict)
      Lb = colMeans(B0)/colMeans(bulk_predict)
      
      div=abs(new.loglik-current.loglik)
      current.loglik <- new.loglik
      
      iter <- iter + 1
      loglik=append(loglik,current.loglik)
      
      print(loglik)
      
      
      
    }
  
  return(list(alpha = alpha_all,beta = beta_all, 
              theta = theta_all, sigma2 = sigma2_all, 
              gamma = gamma_all, pi = pi_all, mubj = mubj_all,
              sigmabj2 = sigmabj_all,loglik = current.loglik))
}
