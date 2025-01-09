###############All variational function and gradient function##############
###########################################################################

ELBO_gamma <- function(gamma,pi,Xsz){
  tmp = as.vector(pi)*(Xsz%*%gamma) - log(1+exp(Xsz%*%gamma)) 
  return(-sum(tmp))
}

ELBO_pi_update <- function(alpha,mubj,sigmabj2,gamma,S,theta,w,Xsz,Xs,Ls){
  term2 = Xs%*%alpha+w%*%mubj+w%*%sigmabj2/2+log(Ls)
  
  expo_term <- 1/theta*exp(term2)
  denom_term <- exp(-Xsz%*%gamma - theta*log(1+expo_term))
  (S == 0)*1/(1+denom_term)
} 

ELBO_pi0 <- function(alpha,mubj,sigmabj2,gamma,S,theta,w,Xsz,Xs,Ls){

  lgamma_diff <- lgamma(S+theta)-lgamma(theta)
  term1 = -Xs%*%alpha-w%*%mubj+w%*%sigmabj2/2-log(Ls)
  term2 = Xs%*%alpha+w%*%mubj+w%*%sigmabj2/2+log(Ls)
  
  logmu = Xs%*%alpha+w%*%mubj+log(Ls)
  common_term1 <- theta*exp(term1)
  common_term2 <- 1/theta*exp(term2)
  
  
  Ftheta = (S>0)*(lgamma_diff - (S+theta)*log(1+common_term1)+
                           theta*log(theta)-theta*logmu)-(S==0)*theta*log(1+common_term2)
  
  denom_term <- exp(-Xsz%*%gamma + Ftheta)
  1/(1+denom_term)
} 
grad_gamma <- function(gamma,pi,Xsz){
  tmp = as.vector(pi)*Xsz - as.vector(exp(Xsz%*%gamma))*Xsz/c(1+exp(Xsz%*%gamma)) 
  
  return(-colSums(tmp))
}
# # hessian_gamma<- function(gamma,pi,Xsz){
#   tmp = Reduce("+",lapply(1:ns,function(i){as.vector(exp(Xsz[i,]%*%gamma)/(1+exp(Xsz[i,]%*%gamma))^2)*tcrossprod(Xsz[i,])}))
#   return(-tmp)
# }
ELBO_mubj <- function(mubj,S,pi,Xs,sigmabj2,sigma2,alpha,theta,w,Ls){
  
  
  Fmubj = c()
  term1 = -Xs%*%alpha-w%*%mubj+w%*%sigmabj2/2-log(Ls)
  term2 = Xs%*%alpha+w%*%mubj+w%*%sigmabj2/2+log(Ls)
  
  tmp1 = (1-pi)*(S>0)*(-(S+theta)*log(1+theta*exp(term1))-theta*w%*%mubj)-
    (1-pi)*(S==0)*theta*log(1+1/theta*exp(term2))
  Fmubj <- -sum(mubj^2/2/sigma2) + sum(tmp1)
  return(-Fmubj)
}

grad_mubj <- function(mubj,S,pi,Xs,sigmabj2,sigma2,alpha,theta,w,Ls){
  
  gradmubj = c()
  
  term1 = -Xs%*%alpha-w%*%mubj+w%*%sigmabj2/2-log(Ls)
  term2 = Xs%*%alpha+w%*%mubj+w%*%sigmabj2/2+log(Ls)
  tmp1 = as.vector((1-pi)*(S>0)*((S+theta)*theta*exp(term1)/(1+theta*exp(term1))-theta))*w-
    as.vector((1-pi)*(S==0)*exp(term2)/(1+1/theta*exp(term2)))*w
  gradmubj <- -mubj/sigma2 + colSums(tmp1)
  
  
  return(-gradmubj)
}

ELBO_sigmabj <- function(sigmabj2,mubj,S,pi,Xs,sigma2,alpha,theta,w,Ls){
  sgmubj = c()
  term1 = -Xs%*%alpha-w%*%mubj+w%*%sigmabj2/2-log(Ls)
  term2 = Xs%*%alpha+w%*%mubj+w%*%sigmabj2/2+log(Ls)
  tmp1 = -(1-pi)*(S>0)*(S+theta)*log(1+theta*exp(term1))-
    (1-pi)*(S==0)*theta*log(1+1/theta*exp(term2))
  sgmubj <- 0.5*sum(log(sigmabj2)-sigmabj2/sigma2) + sum(tmp1)
  return(-sgmubj)
}

grad_sigmabj <- function(sigmabj2,j,mubj,S,pi,Xs,sigma2,alpha,theta,w,Ls){
  #subject_id = unique(cellidx)
  gradsigmabj = c()
  term1 = -Xs%*%alpha-w%*%mubj+w%*%sigmabj2/2-log(Ls)
  term2 = Xs%*%alpha+w%*%mubj+w%*%sigmabj2/2+log(Ls)
  
  tmp1 = as.vector((1-pi)*(S>0)*(S+theta)*theta*exp(term1)/(1+theta*exp(term1)))*w+
    as.vector((1-pi)*(S==0)*exp(term2)/(1+1/theta*exp(term2)))*w
  gradsigmabj <- 0.5*(1/sigmabj2-1/sigma2) - 0.5*colSums(tmp1)
  
  return(-gradsigmabj)
}

ELBO_sigma2 <- function(mubj,sigmabj2,sigma2,mum,sm2){
  tmp = -(log(sigma2)-mum)^2/2/sm2 - length(mubj)/2*log(sigma2) - sum(mubj^2+sigmabj2)/sigma2 
  return(-tmp)
}
grad_sigma2 <- function(mubj,sigmabj2,sigma2,mum,sm2){
  tmp = -1/sigma2*((log(sigma2)-mum)/sm2 - length(mubj)/2 + sum(mubj^2+sigmabj2)/sigma2)
  return(-tmp)
}
#############################################################
ELBO_alpha <- function(alpha,S,B,pi,Xs,WB,ytrt,ysc,Ub=NULL,sigmabj2,sigma2,mubj,theta,beta,w,Ls,Lb,K){
  #subject_id = unique(cellidx)
  alpha0 = alpha[1:K]
  alpha1 = alpha[(K+1):(2*K)]
  alpha2 = alpha[-(1:(2*K))]
  term1 = -Xs%*%alpha-w%*%mubj+w%*%sigmabj2/2-log(Ls)
  term2 = Xs%*%alpha+w%*%mubj+w%*%sigmabj2/2+log(Ls)
  
  Xb = model.matrix(~ytrt)
  alpha.matrix <- rbind(alpha0,alpha1)
  term.bulk = log(rowSums(WB*exp(Xb%*%alpha.matrix)))+log(Lb)+ifelse(!is.null(Ub),Ub%*%alpha2,0)
  
  
  tmp2 =  (-(B+beta)*log(1+1/beta*exp(term.bulk))+ B*term.bulk)
  tmp1 = (-(1-pi)*(S>0)*((S+theta)*log(1+theta*exp(term1))+theta*Xs%*%alpha)-
               (1-pi)*(S==0)*theta*log(1+1/theta*exp(term2)))
  Falpha <- sum(tmp1) + sum(tmp2)
  return(-Falpha)
}


ELBO_alpha_sc <- function(alpha,S,B,pi,Xs,WB,ytrt,ysc,Ub=NULL,sigmabj2,sigma2,mubj,theta,beta,w,Ls,Lb,K){
  #subject_id = unique(cellidx)
  alpha0 = alpha[1:K]
  alpha1 = alpha[(K+1):(2*K)]
  alpha2 = alpha[-(1:(2*K))]
  term1 = -Xs%*%alpha-w%*%mubj+w%*%sigmabj2/2-log(Ls)
  term2 = Xs%*%alpha+w%*%mubj+w%*%sigmabj2/2+log(Ls)
  
  
  tmp1 = -(1-pi)*(S>0)*((S+theta)*log(1+theta*exp(term1))+theta*Xs%*%alpha)-
    (1-pi)*(S==0)*theta*log(1+1/theta*exp(term2))
  Falpha <- sum(tmp1) 
  return(-Falpha)
}
ELBO_alpha_bulk <- function(alpha,S,B,pi,Xs,WB,ytrt,ysc,lambda = 0.1,Ub=NULL,sigmabj2,sigma2,mubj,theta,beta,w,Ls,Lb,K){
  #subject_id = unique(cellidx)
  alpha0 = alpha[1:K]
  alpha1 = alpha[(K+1):(2*K)]
  alpha2 = alpha[-(1:(2*K))]

  
  Xb = model.matrix(~ytrt)
  alpha.matrix <- rbind(alpha0,alpha1)
  term.bulk = log(rowSums(WB*exp(Xb%*%alpha.matrix)))+log(Lb)+ifelse(!is.null(Ub),Ub%*%alpha2,0)
  
  
  tmp2 =  (-(B+beta)*log(1+1/beta*exp(term.bulk))+ B*term.bulk) 
  Falpha <- sum(tmp2)- lambda*sum(WB%*%alpha1^2)
  return(-Falpha)
}
grad_alpha <- function(alpha,S,B,pi,Xs,WB,ytrt,ysc,Ub=NULL,sigmabj2,sigma2,mubj,theta,beta,w,Ls,Lb,K){
  #subject_id = unique(cellidx)
  alpha0 = alpha[1:K]
  alpha1 = alpha[(K+1):(2*K)]
  alpha2 = alpha[-(1:(2*K))]
  Ws = Xs[,1:K]
  Wstrt = Xs[,(K+1):(2*K)]
  Us = Xs[,-(1:(2*K))]
  

  
  
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
    gradalpha <- c(colSums(tmp10) + colSums(tmp20),colSums(tmp11)+colSums(tmp21)) ####to be verified the dimension
    
  }else{
    tmp22 = -(B+beta)/beta*exp(term.bulk)*Ub/(1+1/beta*exp(term.bulk))+B*Ub
    tmp12 = as.vector((1-pi)*(S>0)*(S+theta)*theta*exp(term1))*Us/c(1+theta*exp(term1))-as.vector((1-pi)*(S>0)*theta)*Us-
      as.vector((1-pi)*(S==0)*exp(term2))*Us/c(1+1/theta*exp(term2))
    gradalpha <- c(colSums(tmp10) + colSums(tmp20),colSums(tmp11)+colSums(tmp21), colSums(tmp12)+colSums(tmp22)) ####to be verified the dimension
    
  }
  
  return(-gradalpha)
}

grad_alpha_sc <- function(alpha,S,B,pi,Xs,WB,ytrt,ysc,Ub=NULL,sigmabj2,sigma2,mubj,theta,beta,w,Ls,Lb,K){
  #subject_id = unique(cellidx)
  alpha0 = alpha[1:K]
  alpha1 = alpha[(K+1):(2*K)]
  alpha2 = alpha[-(1:(2*K))]
  Ws = Xs[,1:(2*K)]
  Us = Xs[,-(1:(2*K))]
  
  term1 = -Xs%*%alpha-w%*%mubj+w%*%sigmabj2/2-log(Ls)
  term2 = Xs%*%alpha+w%*%mubj+w%*%sigmabj2/2+log(Ls)
  
  
  tmp10 = as.vector((1-pi)*(S>0)*(S+theta)*theta*exp(term1))*Ws/c(1+theta*exp(term1))-as.vector((1-pi)*(S>0)*theta)*Ws-
    as.vector((1-pi)*(S==0)*exp(term2))*Ws/c(1+1/theta*exp(term2))
  
  gradalpha <- c(colSums(tmp10))
  
  return(-gradalpha)
}

grad_alpha_bulk <- function(alpha,S,B,pi,Xs,WB,ytrt,ysc,lambda = 0.1,Ub=NULL,sigmabj2,sigma2,mubj,theta,beta,w,Ls,Lb,K){
  #subject_id = unique(cellidx)
  alpha0 = alpha[1:K]
  alpha1 = alpha[(K+1):(2*K)]
  alpha2 = alpha[-(1:(2*K))]
  Ws = Xs[,1:K]
  Wstrt = Xs[,(K+1):(2*K)]
  Us = Xs[,-(1:(2*K))]
  
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
    gradalpha <- c(colSums(tmp20),colSums(tmp21)-2*lambda*colSums(WB)*alpha1) ####to be verified the dimension
    
  }else{
    tmp22 = -(B+beta)/beta*exp(term.bulk)*Ub/(1+1/beta*exp(term.bulk))+B*Ub
    gradalpha <- c(colSums(tmp20),colSums(tmp21)-2*lambda*colSums(WB)*alpha1,colSums(tmp22)) ####to be verified the dimension
    
  }
  
  return(-gradalpha)
}

hessian_alpha_weighted <- function(alpha,a1,a2,S,B,pi,Xs,WB,ytrt,ysc,Ub=NULL,sigmabj2,sigma2,mubj,theta,beta,w,Ls,Lb,K,lambda = 0.1){
  
  
  
  #subject_id = unique(cellidx)
  alpha0 = alpha[1:K]
  alpha1 = alpha[(K+1):(2*K)]
  alpha2 = alpha[-(1:(2*K))]
  Vs = Xs[,1:(2*K)]
  Wstrt = Xs[,(K+1):(2*K)]
  Us = Xs[,-(1:(2*K))]
  alpha.matrix <- rbind(alpha0,alpha1)
  
  A0 = -Xs%*%alpha-w%*%mubj+w%*%sigmabj2/2-log(Ls)
  G0 = Xs%*%alpha+w%*%mubj+w%*%sigmabj2/2+log(Ls)
  
  A = theta*exp(A0)
  G = 1/theta*exp(G0)
  
  covariate_term = ifelse(!is.null(Ub),Ub%*%alpha2,0)
  T0 = 1/beta*exp(log(Lb)+ covariate_term)
  
  Xb = model.matrix(~ytrt)
  C = rowSums(WB*exp(Xb%*%alpha.matrix))
  D = cbind(WB,WB*ytrt)*cbind(exp(Xb%*%alpha.matrix),exp(Xb%*%alpha.matrix))
  
  E1 = (WB*exp(Xb%*%alpha.matrix))
  E2 = (WB*ytrt*exp(Xb%*%alpha.matrix))
  E3 = (WB*ytrt^2*exp(Xb%*%alpha.matrix))
  
  cfbulk1 = -(B+beta)/(1+T0*C)
  cfbulk2 = (B+beta)*T0^2/(1+T0*C)^2-B/C^2
  VV3 = Reduce("+",lapply(1:nrow(WB),function(x){
    E = rbind(cbind(diag(E1[x,]),diag(E2[x,])),cbind(diag(E2[x,]),diag(E3[x,])))
    T0[x]*E*cfbulk1[x]+B[x]*E/C[x] + cfbulk2[x]*D[x,]%*%t(D[x,])
  }))
  
  cf1 = as.vector((1-pi)*(S>0)*(S+theta)*A/(1+A)^2)
  
  cf2 = as.vector((1-pi)*(S==0)*theta*G/(1+G)^2)
  VV1 = Reduce("+",lapply(1:nrow(Vs),function(x){-cf1[x]*Vs[x,]%*%t(Vs[x,])}))
  VV2 = Reduce("+",lapply(1:nrow(Vs),function(x){-cf2[x]*Vs[x,]%*%t(Vs[x,])}))
  
  Vd = matrix(0,nrow = ncol(Vs),ncol = ncol(Vs))
  diag(Vd)[(K+1):(2*K)] = 2
  hessianalpha =a1*(VV1+VV2)+a2*VV3-lambda*Vd
  return(hessianalpha)
}

ELBO_theta <- function(theta,alpha,S,pi,Xs,sigmabj2,sigma2,mubj,w,Ls,K){
  
  lgamma_diff <- lgamma(S+theta)-lgamma(theta)
  term1 = -Xs%*%alpha-w%*%mubj+w%*%sigmabj2/2-log(Ls)
  term2 = Xs%*%alpha+w%*%mubj+w%*%sigmabj2/2+log(Ls)
  
  logmu = Xs%*%alpha+w%*%mubj+log(Ls)
  common_term1 <- theta*exp(term1)
  common_term2 <- 1/theta*exp(term2)
  
  
  Ftheta = (1-pi)*(S>0)*(lgamma_diff - (S+theta)*log(1+common_term1)+
                           theta*log(theta)-theta*logmu)-
    (1-pi)*(S==0)*theta*log(1+common_term2)
  Ftheta0 = sum(Ftheta)
  return(-Ftheta0)
}
ELBO_Ls <- function(S,Xs,alpha_all,theta_all,pi_all,sigmabj_all,sigma2_all,mubj_all,w,Ls,K){
  term1 = -Xs%*%t(alpha_all)-w%*%t(mubj_all)+w%*%t(sigmabj_all)/2-log(Ls)
  term2 = Xs%*%t(alpha_all)+w%*%t(mubj_all)+w%*%t(sigmabj_all)/2+log(Ls)
  
  logmu = Xs%*%t(alpha_all)+w%*%t(mubj_all)+log(Ls)
  common_term1 <- as.vector(theta_all)*exp(t(term1))
  common_term2 <- 1/theta_all*exp(t(term2))
  loglsmat = matrix(rep(log(Ls),ng),nrow=ng,byrow=TRUE)
  tmp = (1-pi_all)*(S>0)*(-(S+theta_all)*log(1+common_term1)-theta_all*loglsmat)-
    (1-pi_all)*(S==0)*theta_all*log(1+common_term2)
  return(-sum(tmp))
}

ELBO_lb <- function(S,Xs,alpha,theta,pi,sigmbj2,sigma2,mubj,w,Ls,K){
  
  lgamma_diff <- lgamma(S+theta)-lgamma(theta)
  term1 = -Xs%*%alpha-w%*%mubj+w%*%sigmabj2/2-log(Ls)
  term2 = Xs%*%alpha+w%*%mubj+w%*%sigmabj2/2+log(Ls)
  
  logmu = Xs%*%alpha+w%*%mubj+log(Ls)
  common_term1 <- theta*exp(term1)
  common_term2 <- 1/theta*exp(term2)
  
  
  Ftheta = (1-pi)*(S>0)*(-(S+theta)*log(1+common_term1)-theta*log(Ls))-
    (1-pi)*(S==0)*theta*log(1+common_term2)
  Ftheta0 = sum(Ftheta)
  return(-Ftheta0)
}

grad_Ls <-  function(S,Xs,alpha_all,theta_all,pi_all,sigmabj_all,sigma2_all,mubj_all,w,Ls,K){
  term1 = -Xs%*%t(alpha_all)-w%*%t(mubj_all)+w%*%t(sigmabj_all)/2-log(Ls)
  term2 = Xs%*%t(alpha_all)+w%*%t(mubj_all)+w%*%t(sigmabj_all)/2+log(Ls)
  
  logmu = Xs%*%t(alpha_all)+w%*%t(mubj_all)+log(Ls)
  common_term1 <- as.vector(theta_all)*exp(t(term1))
  common_term2 <- 1/theta_all*exp(t(term2))
  Lsmat = matrix(rep(Ls,ng),nrow=ng,byrow=TRUE)
  tmp = (1-pi_all)*(S>0)*((S+theta_all)*common_term1/(1+common_term1)/Lsmat-theta_all/Lsmat)-
    (1-pi_all)*(S==0)*exp(t(term2-log(Ls)))/(1+common_term2)
  return(-colSums(tmp))
}

grad_Ls_total <- function(Ls,alpha_all,theta_all,pi_all,sigmbj_all,sigma2_all,mubj_all,w,K){
  tmp = rowSums(sapply(1:nrow(alpha_all),function(g){
    theta = theta_all[g]
    alpha = alpha_all[g,]
    pi = pi_all[g,]
    sigmabj2 = sigmabj_all[g,]
    mubj = mubj_all[g,]
    grad_Ls(theta,alpha,S[g,],pi,Xs,sigmabj2,sigma2,mubj,w,Ls,K)
  }))
  return(tmp)
}

grad_theta <- function(theta,alpha,S,pi,Xs,sigmabj2,sigma2,mubj,w,Ls,K){
  dgamma_diff <- digamma(S+theta)-digamma(theta)
  term1 = -Xs%*%alpha-w%*%mubj+w%*%sigmabj2/2-log(Ls)
  term2 = Xs%*%alpha+w%*%mubj+w%*%sigmabj2/2+log(Ls)
  
  logmu = Xs%*%alpha+w%*%mubj+log(Ls)
  common_term1 <- theta*exp(term1)
  common_term2 <- 1/theta*exp(term2)
  
  Ftheta = (1-pi)*(S>0)*(digamma(S+theta)-digamma(theta)-
                           log(1+common_term1)-
                           (S+theta)*exp(term1)/(1+common_term1)+
                           1+log(theta)-logmu)+
    (1-pi)*(S==0)*(exp(term2)/(theta+exp(term2))-log(1+common_term2))
  Ftheta0 = sum(Ftheta)
  return(-Ftheta0)
}


ELBO_beta <- function(beta,alpha,S,B,WB,ytrt,Ub,Lb,K){
  alpha0 = alpha[1:K]
  alpha1 = alpha[(K+1):(2*K)]
  alpha2 = alpha[-(1:(2*K))]
  
  lgamma_diff <- lgamma(B+beta)-lgamma(beta)
  # Common exponential term
  Xb = model.matrix(~ytrt)
  alpha.matrix <- rbind(alpha0,alpha1)
  covariate_term = ifelse(!is.null(Ub),Ub%*%alpha2,0)
  term.bulk = log(rowSums(WB*exp(Xb%*%alpha.matrix)))+log(Lb)+ covariate_term
  
  log_term <- log(1 + exp(term.bulk) / beta)
  
  Fbeta0 = lgamma_diff - (B+beta)*log_term - B*log(beta)
  
  Fbeta = sum(Fbeta0)
  return(-Fbeta)
}


grad_beta <- function(beta,alpha,S,B,WB,ytrt,Ub,Lb,K){
  alpha0 = alpha[1:K]
  alpha1 = alpha[(K+1):(2*K)]
  alpha2 = alpha[-(1:(2*K))]
  
  digamma_diff <- digamma(B + beta) - digamma(beta)
  
  
  
  Xb = model.matrix(~ytrt)
  alpha.matrix <- rbind(alpha0,alpha1)
  covariate_term = ifelse(!is.null(Ub),Ub%*%alpha2,0)
  term.bulk = log(rowSums(WB*exp(Xb%*%alpha.matrix)))+log(Lb)+ covariate_term
  
  
  # Logarithmic term
  log_term <- log(1 + exp(term.bulk) / beta)
  
  # Gradient of log-term with respect to beta
  grad_log_term <- exp(term.bulk) / (beta^2 + beta * exp(term.bulk))
  
  # Full gradient expression
  tmp2 =  digamma_diff - log_term + (B + beta) * grad_log_term - B / beta
  
  Fbeta = sum(tmp2)
  return(-Fbeta)
}

ELBO_Lb <- function(beta,alpha,B,WB,ytrt,Ub,Lb,K){
  alpha0 = alpha[1:K]
  alpha1 = alpha[(K+1):(2*K)]
  alpha2 = alpha[-(1:(2*K))]

  
  
  Xb = model.matrix(~ytrt)
  alpha.matrix <- rbind(alpha0,alpha1)
  covariate_term = ifelse(!is.null(Ub),Ub%*%alpha2,0)
  term.bulk0 = log(rowSums(WB*exp(Xb%*%alpha.matrix)))+ covariate_term
  
  
  log_term <- log(1 + Lb*exp(term.bulk0) / beta)
  Fbeta0 =  - (B+beta)*log_term + B*log(Lb)
  return(-sum(Fbeta0))
}
ELBO_Lb_total <- function(beta_all,alpha_all,B,WB,ytrt,Ub,Lb,K){
  tmp = sum(sapply(1:nrow(alpha_all),function(g){
    beta = beta_all[g]
    alpha = alpha_all[g,]
    ELBO_Lb(beta,alpha,B[g,],WB,ytrt,Ub,Lb,K)
  }))
  return(tmp)
}
grad_Lb <-  function(beta,alpha,B,WB,ytrt,Ub,Lb,K){
  alpha0 = alpha[1:K]
  alpha1 = alpha[(K+1):(2*K)]
  alpha2 = alpha[-(1:(2*K))]
  
  
  Xb = model.matrix(~ytrt)
  alpha.matrix <- rbind(alpha0,alpha1)
  covariate_term = ifelse(!is.null(Ub),Ub%*%alpha2,0)
  term.bulk0 = log(rowSums(WB*exp(Xb%*%alpha.matrix)))+ covariate_term
  
  
  log_term <- log(1 + exp(term.bulk0)*Lb / beta)
  grad_log_term <- exp(term.bulk0) / (beta + Lb*exp(term.bulk0))
  
  Fb =  - (B+beta)*grad_log_term + B/Lb
  return(-Fb)
}
grad_Lb_total <- function(beta_all,alpha_all,B,WB,ytrt,Ub,Lb,K){
  tmp = rowSums(sapply(1:nrow(alpha_all),function(g){
    beta = beta_all[g]
    alpha = alpha_all[g,]
    grad_Lb(beta,alpha,B[g,],WB,ytrt,Ub,Lb,K)
  }))
  return(as.vector(tmp))
}

ELBO_total <- function(S,B,pi,Xs,Xsz,gamma,WB,ytrt,ysc,Ub=NULL,sigmabj2,sigma2,mubj,alpha,theta,beta,w,Ls,Lb,K){
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
  
  tmp2 = lgamma(B+beta)-lgamma(beta)-lgamma(B+1)-(B+beta)*log(1+1/beta*exp(term.bulk))+B*term.bulk - B*log(beta)
  
  
  
  lgamma_diff <- lgamma(S+theta)-lgamma(theta)-lgamma(S+1)
  
  
  logmu =  Xs%*%alpha+w%*%mubj+log(Ls)
  common_term1 <- theta*exp(term1)
  common_term2 <- 1/theta*exp(term2)
  
  tmp1 = (1-pi)*(S>0)*(lgamma_diff - (S+theta)*log(1+common_term1)+theta*log(theta)-theta*logmu)-(1-pi)*(S==0)*theta*log(1+common_term2)
  
  Ftotal <- sum(tmp1) + sum(tmp2) + Fspart + Fzero + Fpi
  return(-Ftotal)
}


ELBO_total_sc <- function(S,B,pi,Xs,Xsz,gamma,WB,ytrt,ysc,Ub=NULL,sigmabj2,sigma2,mubj,alpha,theta,beta,w,Ls,Lb,K){
  alpha0 = alpha[1:K]
  alpha1 = alpha[(K+1):(2*K)]
  alpha2 = alpha[-(1:(2*K))]
  #####For zero inflation part#############
  ##########################################
  
  Fspart = 0.5*sum((log(sigmabj2)-log(sigma2)-(mubj^2+sigmabj2)/sigma2))
  
  
  Fzero = sum(pi*Xsz%*%gamma - log(1+exp(Xsz%*%gamma)))
  
  
  Fpi = -sum(pi[pi>0]*log(pi[pi>0])) -sum((1-pi[pi<1])*log(1-pi[pi<1]))
  
  # Xb = model.matrix(~ytrt)
  # alpha.matrix <- rbind(alpha0,alpha1)
  # covariate_term = ifelse(!is.null(Ub),Ub%*%alpha2,0)
  # term.bulk = log(rowSums(WB*exp(Xb%*%alpha.matrix)))+log(Lb)+ covariate_term
  
  term1 = -Xs%*%alpha-w%*%mubj+w%*%sigmabj2/2-log(Ls)
  term2 = Xs%*%alpha+w%*%mubj+w%*%sigmabj2/2+log(Ls)
  
  # tmp2 = lgamma(B+beta)-lgamma(beta)-lgamma(B+1)-(B+beta)*log(1+1/beta*exp(term.bulk))+B*term.bulk - B*log(beta)
  
  
  
  lgamma_diff <- lgamma(S+theta)-lgamma(theta)-lgamma(S+1)
  
  
  logmu =  Xs%*%alpha+w%*%mubj+log(Ls)
  common_term1 <- theta*exp(term1)
  common_term2 <- 1/theta*exp(term2)
  
  tmp1 = (1-pi)*(S>0)*(lgamma_diff - (S+theta)*log(1+common_term1)+theta*log(theta)-theta*logmu)-(1-pi)*(S==0)*theta*log(1+common_term2)
  
  Ftotal <- sum(tmp1) + Fspart + Fzero + Fpi
  return(-Ftotal)
}

ELBO_total_bulk <- function(S,B,pi,Xs,Xsz,gamma,WB,ytrt,ysc,Ub=NULL,lambda = 0.1,sigmabj2,sigma2,mubj,alpha,theta,beta,w,Ls,Lb,K){
  alpha0 = alpha[1:K]
  alpha1 = alpha[(K+1):(2*K)]
  alpha2 = alpha[-(1:(2*K))]
  #####For zero inflation part#############
  ##########################################
  
  
  Xb = model.matrix(~ytrt)
  alpha.matrix <- rbind(alpha0,alpha1)
  covariate_term = ifelse(!is.null(Ub),Ub%*%alpha2,0)
  term.bulk = log(rowSums(WB*exp(Xb%*%alpha.matrix)))+log(Lb)+ covariate_term
  
  tmp2 = lgamma(B+beta)-lgamma(beta)-lgamma(B+1)-(B+beta)*log(1+1/beta*exp(term.bulk))+B*term.bulk - B*log(beta) 
  
  
  

  Ftotal <- sum(tmp2)- lambda*sum(WB%*%alpha1^2)
  return(-Ftotal)
}

ELBO_alpha_hessian <- function(alpha,result,S,B,Xs,Xsz,WB,ytrt,ysc,Ub=NULL,w,Ls,Lb,K){
  beta = result$beta
  theta = result$theta
  sigma2 = result$sigma2
  alpha = result$alpha
  gamma = result$gamma
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  
  #subject_id = unique(cellidx)
  alpha0 = alpha[1:K]
  alpha1 = alpha[(K+1):(2*K)]
  alpha2 = alpha[-(1:(2*K))]
  
  Fspart = 0.5*sum((log(sigmabj2)-log(sigma2)-(mubj^2+sigmabj2)/sigma2))
  
  
  Fzero = sum(pi*Xsz%*%gamma - log(1+exp(Xsz%*%gamma)))
  
  
  Fpi = -sum(pi[pi>0]*log(pi[pi>0])) -sum((1-pi[pi<1])*log(1-pi[pi<1]))
  
  Xb = model.matrix(~ytrt)
  alpha.matrix <- rbind(alpha0,alpha1)
  covariate_term = ifelse(!is.null(Ub),Ub%*%alpha2,0)
  term.bulk = log(rowSums(WB*exp(Xb%*%alpha.matrix)))+log(Lb)+ covariate_term
  
  term1 = -Xs%*%alpha-w%*%mubj+w%*%sigmabj2/2-log(Ls)
  term2 = Xs%*%alpha+w%*%mubj+w%*%sigmabj2/2+log(Ls)
  
  tmp2 = lgamma(B+beta)-lgamma(beta)-lgamma(B+1)-(B+beta)*log(1+1/beta*exp(term.bulk))+B*term.bulk - B*log(beta)
  
  
  
  lgamma_diff <- lgamma(S+theta)-lgamma(theta)-lgamma(S+1)
  
  
  logmu =  Xs%*%alpha+w%*%mubj+log(Ls)
  common_term1 <- theta*exp(term1)
  common_term2 <- 1/theta*exp(term2)
  
  tmp1 = (1-pi)*(S>0)*(lgamma_diff - (S+theta)*log(1+common_term1)+theta*log(theta)-theta*logmu)-(1-pi)*(S==0)*theta*log(1+common_term2)
  
  Ftotal <- sum(tmp1) + sum(tmp2) + Fspart + Fzero + Fpi
  return(-Ftotal)
}
ELBO_partj_sc <- function(param,psi,S_ind,B,Xs_ind,Xsz_ind,WB,ytrt,ysc_ind,Ub=NULL,
                       Ls_ind,Lb,K){
  
  beta = param[1]
  theta = param[2]
  sigma2 = param[3]
  gamma = param[4:(3+ncol(Xsz_ind))]
  alpha0 = param[(4+ncol(Xsz_ind)):(3+K+ncol(Xsz_ind))]
  alpha1 = param[(4+K+ncol(Xsz_ind)):(3+2*K+ncol(Xsz_ind))]
  if(!is.null(Ub)){
    alpha2 = param[(4+2*K+ncol(Xsz_ind)):(3+2*K+ncol(Xsz_ind)+ncol(Ub))]
    alpha = c(alpha0,alpha1,alpha2)
  }else{
    alpha = c(alpha0,alpha1)
  }
  
  mubj_ind = psi[1]
  sigmabj2_ind = psi[2]
  pi_ind = psi[3:(2+length(S_ind))]
  
  #####For zero inflation part#############
  ##########################################
  Fspart = 0.5*(log(sigmabj2_ind)-log(sigma2)-(mubj_ind^2+sigmabj2_ind)/sigma2)
  
 
  Fzero = pi_ind*Xsz_ind%*%gamma - log(1+exp(Xsz_ind%*%gamma))
  
  Fpi = -sum(pi_ind[pi_ind>0]*log(pi_ind[pi_ind>0])) -sum((1-pi_ind[pi_ind<1])*log(1-pi_ind[pi_ind<1]))
  
  # Xb = model.matrix(~ytrt)
  # alpha.matrix <- rbind(alpha0,alpha1)
  # covariate_term = ifelse(!is.null(Ub),Ub%*%alpha2,0)
  # term.bulk = log(rowSums(WB*exp(Xb%*%alpha.matrix)))+log(Lb)+ covariate_term
  # tmp2 = lgamma(B+beta)-lgamma(beta)-lgamma(B+1)-(B+beta)*log(1+1/beta*exp(term.bulk))+B*term.bulk- B*log(beta)
  # 
  
    lgamma_diff <- lgamma(S_ind+theta)-lgamma(theta)-lgamma(S_ind+1)
    term1 <- -Xs_ind%*%alpha-mubj_ind+sigmabj2_ind/2 -log(Ls_ind)
    term2 <- Xs_ind%*%alpha+mubj_ind+sigmabj2_ind/2 + log(Ls_ind)
    
    logmu = Xs_ind%*%alpha+mubj_ind+log(Ls_ind)
    common_term1 <- theta*exp(term1)
    common_term2 <- 1/theta*exp(term2)
    
    tmp1 = (1-pi_ind)*(S_ind>0)*(lgamma_diff - (S_ind+theta)*log(1+common_term1)+theta*log(theta)-theta*logmu)-(1-pi_ind)*(S_ind==0)*theta*log(1+common_term2)
    
  
  Ftotal <- sum(tmp1) + sum(Fspart) + sum(Fzero) + Fpi
  return(Ftotal)
}


ELBO_partj_bulk <- function(param,S_ind,B,Xs_ind,Xsz_ind,WB,ytrt,ysc_ind,Ub=NULL,
                          Ls_ind,Lb,K){
  
  beta = param[1]
  theta = param[2]
  sigma2 = param[3]
  gamma = param[4:(3+ncol(Xsz_ind))]
  alpha0 = param[(4+ncol(Xsz_ind)):(3+K+ncol(Xsz_ind))]
  alpha1 = param[(4+K+ncol(Xsz_ind)):(3+2*K+ncol(Xsz_ind))]
  if(!is.null(Ub)){
    alpha2 = param[(4+2*K+ncol(Xsz_ind)):(3+2*K+ncol(Xsz_ind)+ncol(Ub))]
    alpha = c(alpha0,alpha1,alpha2)
  }else{
    alpha = c(alpha0,alpha1)
  }
  
  # mubj_ind = psi[1]
  # sigmabj2_ind = psi[2]
  # pi_ind = psi[3:(2+length(S_ind))]
  
  #####For zero inflation part#############
  ##########################################
  # Fspart = 0.5*(log(sigmabj2_ind)-log(sigma2)-(mubj_ind^2+sigmabj2_ind)/sigma2)
  # 
  # 
  # Fzero = pi_ind*Xsz_ind%*%gamma - log(1+exp(Xsz_ind%*%gamma))
  # 
  # Fpi = -sum(pi_ind[pi_ind>0]*log(pi_ind[pi_ind>0])) -sum((1-pi_ind[pi_ind<1])*log(1-pi_ind[pi_ind<1]))
  # 
  Xb = model.matrix(~ytrt)
  alpha.matrix <- rbind(alpha0,alpha1)
  covariate_term = ifelse(!is.null(Ub),Ub%*%alpha2,0)
  term.bulk = log(rowSums(WB*exp(Xb%*%alpha.matrix)))+log(Lb)+ covariate_term
  tmp2 = lgamma(B+beta)-lgamma(beta)-lgamma(B+1)-(B+beta)*log(1+1/beta*exp(term.bulk))+B*term.bulk- B*log(beta)

  
  # lgamma_diff <- lgamma(S_ind+theta)-lgamma(theta)-lgamma(S_ind+1)
  # term1 <- -Xs_ind%*%alpha-mubj_ind+sigmabj2_ind/2 -log(Ls_ind)
  # term2 <- Xs_ind%*%alpha+mubj_ind+sigmabj2_ind/2 + log(Ls_ind)
  # 
  # logmu = Xs_ind%*%alpha+mubj_ind+log(Ls_ind)
  # common_term1 <- theta*exp(term1)
  # common_term2 <- 1/theta*exp(term2)
  # 
  # tmp1 = (1-pi_ind)*(S_ind>0)*(lgamma_diff - (S_ind+theta)*log(1+common_term1)+theta*log(theta)-theta*logmu)-(1-pi_ind)*(S_ind==0)*theta*log(1+common_term2)
  # 
  
  Ftotal <- sum(tmp2) #+ sum(Fspart) + sum(Fzero) + Fpi
  return(Ftotal)
}

grad_partj_bulk<- function(param,S_ind,B,Xs_ind,Xsz_ind,WB,ytrt,ysc_ind,Ub=NULL,
                           Ls_ind,Lb,K){
  beta = param[1]
  theta = param[2]
  sigma2 = param[3]
  gamma = param[4:(3+ncol(Xsz_ind))]
  alpha0 = param[(4+ncol(Xsz_ind)):(3+K+ncol(Xsz_ind))]
  alpha1 = param[(4+K+ncol(Xsz_ind)):(3+2*K+ncol(Xsz_ind))]
  if(!is.null(Ub)){
    alpha2 = param[(4+2*K+ncol(Xsz_ind)):(3+2*K+ncol(Xsz_ind)+ncol(Ub))]
    alpha = c(alpha0,alpha1,alpha2)
  }else{
    alpha = c(alpha0,alpha1)
  }
  digamma_diff <- digamma(B + beta) - digamma(beta)
  
  Xb = model.matrix(~ytrt)
  alpha.matrix <- rbind(alpha0,alpha1)
  covariate_term = ifelse(!is.null(Ub),Ub%*%alpha2,0)
  term.bulk = log(rowSums(WB*exp(Xb%*%alpha.matrix)))+log(Lb)+ covariate_term
  #tmp2 = lgamma(B+beta)-lgamma(beta)-lgamma(B+1)-(B+beta)*log(1+1/beta*exp(term.bulk))+B*term.bulk- B*log(beta)
  term.bulk.gr1 = WB*exp(Xb%*%alpha.matrix)*Lb*exp(covariate_term)
  
  
  tmp20 = -(B+beta)*1/beta*term.bulk.gr1/c(1+1/beta*exp(term.bulk))+
    B*WB*exp(Xb%*%alpha.matrix)/rowSums(WB*exp(Xb%*%alpha.matrix))
  tmp21 =  -(B+beta)*1/beta*term.bulk.gr1*ytrt/c(1+1/beta*exp(term.bulk))+
    B*ytrt*WB*exp(Xb%*%alpha.matrix)/rowSums(WB*exp(Xb%*%alpha.matrix))
  
  if(is.null(Ub)){
    gradalpha <- cbind(tmp20,tmp21) ####to be verified the dimension
    
  }else{
    tmp22 = -(B+beta)/beta*exp(term.bulk)*Ub/(1+1/beta*exp(term.bulk))+B*Ub
    gradalpha <- cbind(tmp20,tmp21,tmp22) ####to be verified the dimension
    
  }
  
  # Logarithmic term
  log_term <- log(1 + exp(term.bulk) / beta)
  
  # Gradient of log-term with respect to beta
  grad_log_term <- exp(term.bulk) / (beta^2 + beta * exp(term.bulk))
  
  # Full gradient expression
  tmp2 =  digamma_diff - log_term + (B + beta) * grad_log_term  - B / beta
  
  final = cbind(tmp2,0,0,0,0,gradalpha)
  return(final)
}
ELBO_combined <- function(parampsi,S_ind,B,Xs_ind,Xsz_ind,WB,ytrt,ysc_ind,Ub=NULL,
                          Ls_ind,Lb,K){
  psi = parampsi[1:(2+length(S_ind))]
  param = parampsi[-(1:(2+length(S_ind)))]
  
  ELBO_partj_sc(param,psi,S_ind,B,Xs_ind,Xsz_ind,WB,ytrt,ysc_ind,Ub,
             Ls_ind,Lb,K)
}






sandwich_cov <- function(result,nsub,a1,a2,S,B,Xs,Xsz0,WB,ytrt,ysc,
                         Ub=NULL,cellidx,Ls,Lb,K) {
  library(parallel)
  library(MASS)
  library(numDeriv)
  nb = length(B)
  beta = result$beta
  theta = result$theta
  sigma2 = result$sigma2
  alpha = result$alpha
  gamma = result$gamma
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  
  param <- c(beta,theta,sigma2,gamma,alpha)
  param_ln <- length(param)
  system.time(abs <- foreach(t = 1:nsub, .combine='rbind', .multicombine=TRUE,
                             .export = c("ELBO_partj_sc","ELBO_combined","ELBO_partj_bulk","grad_partj_bulk")) %dopar%{
    #nodes=nodes
    library(numDeriv)
    
    index = which(cellidx == t)
    mubj_ind = mubj[t]
    sigmabj2_ind = sigmabj2[t]
    pi_ind = pi[index]
    psi = c(mubj_ind,sigmabj2_ind,pi_ind)
    parampsi = c(psi,param)
    psi_ln <- length(psi)
    
    grad <- grad(ELBO_partj_sc,x = param,psi = psi,S_ind = S[index],
                 B = B,Xs_ind = Xs[index,],Xsz_ind=Xsz0[index,],
                 WB = WB,ytrt = ytrt,ysc_ind = ysc[index],Ub=NULL,
                 Ls_ind = Ls[index],Lb = Lb,K = K)
    hess <- hessian(ELBO_combined,x=parampsi,S_ind = S[index],
                    B = B,Xs_ind = Xs[index,],Xsz_ind=Xsz0[index,],
                    WB = WB,ytrt = ytrt,ysc_ind = ysc[index],Ub=NULL,
                    Ls_ind = Ls[index],Lb = Lb,K = K)
    
    Ai <- hess[-(1:psi_ln),-(1:psi_ln)] - hess[-(1:psi_ln), 1:psi_ln] %*% solve(hess[1:psi_ln, 1:psi_ln]) %*% hess[1:psi_ln, -(1:psi_ln)]
    Bi <- tcrossprod(grad) 
    c(c(Ai), c(Bi))
  })
  
  #abs <- matrix(unlist(abs), nrow=length(abs), byrow=TRUE)
  abhat <- colMeans(abs)
  ahat <- matrix(abhat[1:param_ln^2], nrow=param_ln)
  bhat <- matrix(abhat[-(1:param_ln^2)], nrow=param_ln)
  
  
  bulkhess <- hessian(ELBO_partj_bulk,x=param,S_ind = S,
                     B = B,Xs_ind = Xs,Xsz_ind=Xsz0,
                     WB = WB,ytrt = ytrt,ysc_ind = ysc,Ub=NULL,
                     Ls_ind = Ls,Lb = Lb,K = K)
  A_final <- (a1*ahat+a2*bulkhess)/(a2*nb+a1*nsub)
  gr_bulk <- grad_partj_bulk(param = param,S_ind = S,
                             B = B,Xs_ind = Xs,Xsz_ind=Xsz0,
                             WB = WB,ytrt = ytrt,ysc_ind = ysc,Ub=NULL,
                             Ls_ind = Ls,Lb = Lb,K = K)
  
  B4_final <- Reduce("+",lapply(1:nrow(gr_bulk),function(x){tcrossprod(gr_bulk[x,])}))
  
  B_final <- (a1^2*bhat + a2^2*B4_final)/(a1*nsub+a2*nb)
  
  A_inv = ginv(A_final)
  sand_cov= A_inv%*% B_final %*% A_inv /(a1*nsub+a2*nb)
  
  
  
  return(list(ahat=A_final, bhat=B_final, sand_cov = sand_cov, param = param))
}



DE_test <- function(sand_cov,param.est,C,b0=0){
  cov.inf <- C%*%sand_cov%*%t(C)
  tes_val <- t(matrix(C%*%param.est)-b0)%*%ginv(cov.inf)%*%matrix(C%*%param.est-b0)
  1-pchisq(as.numeric(tes_val),df=rankMatrix(C))
}


#empirical_cov <- function()

variational_EM <- function(S,B,prior_alpha,param_initial = NULL,Xs,WB,Ub=NULL,ysc,
                            ytrt,cellidx,Xsz,Lb,
                            Ls,K,eps = 1e-9,max.iter = 200, trace = 0){
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
  
  
  llgamma = llmubj = llsigmabj = lltheta = llbeta = llalpha = -1e6
  
  while((div> eps*(abs(current.loglik)+eps)) && iter <= max.iter) {
    
    
    print(iter)
    # eps.gam = 1e-3
    # delta.gamma = gamma_current; p.iter = 1
    # while(!all(delta.gamma < eps.gam) & (p.iter < p.max.iter)){
    #   print(p.iter)
    q=try(nlminb(start =as.vector(gamma_current), ELBO_gamma,
                 gradient = grad_gamma,hessian = hessian_gamma,lower = rep(-Inf,ncol(Xsz)), upper = rep(Inf,ncol(Xsz)),
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
    
    
    sigma2_new=sum(sigmabj_new+mubj_new^2)/nsub
    
    #   delta.sigma2 <- abs(sigma2_new-sigma2_current)
    #   sigma2_current <- sigma2_new
    #   p.iter <- p.iter+1
    # }    
    
    
    
    
    q=try(nlminb(start =as.vector(alpha_current), ELBO_alpha, gradient = grad_alpha,
                 lower = rep(-Inf,2*K), upper = rep(Inf,2*K),
                 S=S,B=B,pi=pi_new,Xs=Xs,WB = WB, ytrt = ytrt,ysc = ysc,
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
    print(llalpha)
    
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
    
    
    objective =ELBO_total(S=S,B=B,pi=pi_new,
                          Xs=Xs,Xsz=Xsz,gamma=gamma_new,
                          WB = WB,ytrt = ytrt,ysc = ysc,alpha = alpha_new,
                          Ub=NULL,sigmabj2 = sigmabj_new,sigma2=sigma2_new,
                          mubj=mubj_new,theta=theta_new,beta = beta_new,
                          w = w,Ls = Ls,Lb = Lb,K = K)
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
              sigmabj2 = sigmabj_new))
}


variational_EM_sc <- function(S,B,prior_alpha,param_initial = NULL,Xs,WB,Ub=NULL,ysc,
                              ytrt,cellidx,Xsz,Lb,
                              Ls,K,eps = 1e-9,max.iter = 200, trace = 0){
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
  
  
  llgamma = llmubj = llsigmabj = lltheta = llbeta = llalpha = -1e6
  
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
    
    
    sigma2_new=sum(sigmabj_new+mubj_new^2)/nsub
    
    #   delta.sigma2 <- abs(sigma2_new-sigma2_current)
    #   sigma2_current <- sigma2_new
    #   p.iter <- p.iter+1
    # }    
    
    
    
    
    q=try(nlminb(start =as.vector(alpha_current), ELBO_alpha_sc, gradient = grad_alpha_sc,
                 lower = rep(-Inf,2*K), upper = rep(Inf,2*K),
                 S=S,B=B,pi=pi_new,Xs=Xs,WB = WB, ytrt = ytrt,ysc = ysc,
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
    print(llalpha)
    
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
    
    
    beta_new = beta_current
    
    objective =ELBO_total_sc(S=S,B=B,pi=pi_new,
                             Xs=Xs,Xsz=Xsz,gamma=gamma_new,
                             WB = WB,ytrt = ytrt,ysc = ysc,alpha = alpha_new,
                             Ub=NULL,sigmabj2 = sigmabj_new,sigma2=sigma2_new,
                             mubj=mubj_new,theta=theta_new,beta = beta_new,
                             w = w,Ls = Ls,Lb = Lb,K = K)
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
              sigmabj2 = sigmabj_new))
}

variational_EM_bulk <- function(S,B,prior_alpha,param_initial = NULL,Xs,WB,Ub=NULL,ysc,
                                ytrt,cellidx,Xsz,Lb,
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
  
  
  llgamma = llmubj = llsigmabj = lltheta = llbeta = llalpha = -1e20
  
  while((div> eps*(abs(current.loglik)+eps)) && iter <= max.iter) {
    
    
    theta_new = theta_current 
    mubj_new = mubj_current 
    sigmabj_new = sigmabj_current 
    sigma2_new =  sigma2_current 
    pi_new=  pi_current 
    gamma_new  = gamma_current
    
    q=try(nlminb(start =as.vector(alpha_current), ELBO_alpha_bulk, gradient = grad_alpha_bulk,
                 lower = rep(-Inf,2*K), upper = rep(Inf,2*K),
                 S=S,B=B,pi=pi_new,Xs=Xs,WB = WB, ytrt = ytrt,ysc = ysc,
                 Ub=NULL,mubj =mubj_new,sigmabj2 = sigmabj_new,sigma2 = sigma2_new,
                 theta = theta_current,beta = beta_current,w = w,
                 Ls = Ls,Lb = Lb,K=K,lambda=lambda,control = list(trace = 1, iter.max = 200)), silent = TRUE)
    
    
    
    if(!inherits(q, "try-error")) {
      if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
      if(llalpha > -q$objective) { if(trace) cat("Optimization alpha did not improve estimates of variational parameters.\n"); alpha_new <- alpha_current }
      else {
        if(trace) cat("Variational parameters alpha updated", "\n")
        alpha_new <- q$par
      }
    }else { alpha_new <- alpha_current}
    llalpha <- -q$objective
    print(llalpha)
    
    q=try(nlminb(start =as.vector(beta_current), ELBO_beta, gradient = grad_beta,
                 hessian = hessian_beta0, lower = 0, upper = Inf,
                 alpha = alpha_new,B=B,WB=WB,ytrt=ytrt,Ub=NULL,
                 Lb=Lb,K=K,control = list(trace = 1, iter.max = 100)), silent = TRUE)
    
    
    if(!inherits(q, "try-error")) {
      if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
      if(llbeta > -q$objective) { if(trace) cat("Optimization beta did not improve estimates of variational parameters.\n"); beta_new <- beta_current }
      else {
        if(trace) cat("Variational parameters beta updated", "\n")
        beta_new <- q$par
      }
    }else {beta_new <- beta_current}
    print(q$objective)
    llbeta<- -q$objective
    
    
    objective =ELBO_total_bulk(S=S,B=B,pi=pi_new,
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
              sigmabj2 = sigmabj_new))
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
                   hessian = hessian_alpha_weighted,
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
      
      
      q=try(nlminb(start =as.vector(beta_current), ELBO_beta, gradient = grad_beta,
                   hessian = hessian_beta0,lower = 0, upper = Inf,
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
      
      # alpha.matrix = rbind(alpha_new[1:(K)],alpha_new[(K+1):(2*K)])
      # Xb = model.matrix(~ytrt)
      # mus_predict <- exp(Xs%*%alpha_new + mubj_new)*(1-pi_new)
      # mub_predict <- rowSums(WB*exp(Xb%*%alpha.matrix))
      # 
      # sc_predict[g,] = mus_predict
      # bulk_predict[g,] = mub_predict
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
    
    q=nlminb(start =as.vector(Ls), ELBO_Ls, gradient = grad_Ls,
                 lower = rep(0,ns), upper = rep(Inf,ns),S=Y,
                 alpha_all = alpha_all,theta_all=theta_all,pi_all=pi_all,
                 sigmabj_all = sigmabj_all,sigma2_all = sigma2_all,
                 mubj_all=mubj_all,w=w,K=K,Xs = Xs,
                 control = list(trace = 0, iter.max = 100))
    Ls = q$par
    q2=nlminb(start =as.vector(Lb), ELBO_Lb_total, gradient = grad_Lb_total,
                  lower = rep(0,nb), upper = rep(Inf,nb),B =B0,
                  alpha_all = alpha_all,beta_all=beta_all,
                  WB = WB,ytrt = ytrt,Ub = Ub,K = K,
                  control = list(trace = 0, iter.max = 100))
    Lb = q2$par
    # if(!inherits(q, "try-error")) {
    #   if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
    #   if(llbeta > -q$objective) { if(trace) cat("Optimization beta did not improve estimates of variational parameters.\n"); beta_new <- beta_current }
    #   else {
    #     if(trace) cat("Variational parameters beta updated", "\n")
    #     beta_new <- q$par
    #   }
    # }else {beta_new <- beta_current}
    
    #   Ls = colMeans(Y)/colMeans(sc_predict)
    # Lb = colMeans(B0)/colMeans(bulk_predict)
    
    div=abs(new.loglik-current.loglik)
    current.loglik <- new.loglik
    
    iter <- iter + 1
    loglik=append(loglik,current.loglik)
    
    print(loglik)
    
    
    
  }
  
  return(list(alpha = alpha_all,beta = beta_all, 
              theta = theta_all, sigma2 = sigma2_all, 
              gamma = gamma_all, pi = pi_all, mubj = mubj_all,
              sigmabj2 = sigmabj_all,loglik = current.loglik,Ls = Ls,Lb=Lb))
}
