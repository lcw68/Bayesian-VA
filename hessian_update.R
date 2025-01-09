########alpha related function##############
############################################

hessian_alpha <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz){
  
  
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
  Vs = Xs[,1:(2*K)]
  Wstrt = Xs[,(K+1):(2*K)]
  Us = Xs[,-(1:(2*K))]
  alpha.matrix <- rbind(alpha0,alpha1)
  
  A0 = -Xs%*%alpha-mubj+sigmabj2/2-log(Ls)
  G0 = Xs%*%alpha+mubj+sigmabj2/2+log(Ls)
  
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
  
  hessianalpha =a1*(VV1+VV2)+a2*VV3
  return(hessianalpha)
}

hessian_alpha_bulk <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz,lambda = 0.1){
  
  
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
  Vs = Xs[,1:(2*K)]
  Wstrt = Xs[,(K+1):(2*K)]
  Us = Xs[,-(1:(2*K))]
  alpha.matrix <- rbind(alpha0,alpha1)
  
  # A0 = -Xs%*%alpha-mubj+sigmabj2/2-log(Ls)
  # G0 = Xs%*%alpha+mubj+sigmabj2/2+log(Ls)
  # 
  # A = theta*exp(A0)
  # G = 1/theta*exp(G0)
  
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
  
  # cf1 = as.vector((1-pi)*(S>0)*(S+theta)*A/(1+A)^2)
  # 
  # cf2 = as.vector((1-pi)*(S==0)*theta*G/(1+G)^2)
  # VV1 = Reduce("+",lapply(1:nrow(Vs),function(x){-cf1[x]*Vs[x,]%*%t(Vs[x,])}))
  # VV2 = Reduce("+",lapply(1:nrow(Vs),function(x){-cf2[x]*Vs[x,]%*%t(Vs[x,])}))
  Vd = matrix(0,nrow = ncol(Vs),ncol = ncol(Vs))
  diag(Vd)[(K+1):(2*K)] = 2*colSums(WB)
  
  hessianalpha =VV3-lambda*Vd
  return(hessianalpha)
}

hessian_alpha_sc_ind <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz){
  
  
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
  Vs = Xs[,1:(2*K)]
  Wstrt = Xs[,(K+1):(2*K)]
  Us = Xs[,-(1:(2*K))]
  alpha.matrix <- rbind(alpha0,alpha1)
  
  A0 = -Xs%*%alpha-mubj+sigmabj2/2-log(Ls)
  G0 = Xs%*%alpha+mubj+sigmabj2/2+log(Ls)
  
  A = theta*exp(A0)
  G = 1/theta*exp(G0)
  
  # covariate_term = ifelse(!is.null(Ub),Ub%*%alpha2,0)
  # T0 = 1/beta*exp(log(Lb)+ covariate_term)
  # 
  # Xb = model.matrix(~ytrt)
  # C = rowSums(WB*exp(Xb%*%alpha.matrix))
  # D = cbind(WB,WB*ytrt)*cbind(exp(Xb%*%alpha.matrix),exp(Xb%*%alpha.matrix))
  # 
  # E1 = (WB*exp(Xb%*%alpha.matrix))
  # E2 = (WB*ytrt*exp(Xb%*%alpha.matrix))
  # E3 = (WB*ytrt^2*exp(Xb%*%alpha.matrix))
  # 
  # cfbulk1 = -(B+beta)/(1+T0*C)
  # cfbulk2 = (B+beta)*T0^2/(1+T0*C)^2-B/C^2
  # VV3 = Reduce("+",lapply(1:nrow(WB),function(x){
  #   E = rbind(cbind(diag(E1[x,]),diag(E2[x,])),cbind(diag(E2[x,]),diag(E3[x,])))
  #   T0[x]*E*cfbulk1[x]+B[x]*E/C[x] + cfbulk2[x]*D[x,]%*%t(D[x,])
  # }))
  
  cf1 = as.vector((1-pi)*(S>0)*(S+theta)*A/(1+A)^2)
  
  cf2 = as.vector((1-pi)*(S==0)*theta*G/(1+G)^2)
  VV1 = Reduce("+",lapply(1:nrow(Vs),function(x){-cf1[x]*Vs[x,]%*%t(Vs[x,])}))
  VV2 = Reduce("+",lapply(1:nrow(Vs),function(x){-cf2[x]*Vs[x,]%*%t(Vs[x,])}))
  
  hessianalpha =VV1+VV2
  return(hessianalpha)
}

hessian_alpha_theta_ind <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz)
{
  # ord_id=table(cellidx)
  # nsub = length(ord_id)
  # w_l=list()
  # for(i in 1:length(ord_id)){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  
  alpha =  result$alpha
  theta = result$theta
  sigma2 = result$sigma2
  gamma = result$gamma
  beta = result$beta
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  
  #subject_id = unique(cellidx)
  alpha0 = alpha[1:K]
  alpha1 = alpha[(K+1):(2*K)]
  alpha2 = alpha[-(1:(2*K))]
  Vs = Xs[,1:(2*K)]
  Wstrt = Xs[,(K+1):(2*K)]
  Us = Xs[,-(1:(2*K))]
  alpha.matrix <- rbind(alpha0,alpha1)
  
  A0 = -Xs%*%alpha-mubj+sigmabj2/2-log(Ls)
  G0 = Xs%*%alpha+mubj+sigmabj2/2+log(Ls)
  
  A = theta*exp(A0)
  G = 1/theta*exp(G0)
  
  
  cf1 = as.vector((1-pi)*(S>0)*(exp(A0)*S-1)/(1+A)^2)
  
  cf2 = as.vector((1-pi)*(S==0)*G^2/(1+G)^2)
  grad2 = cf1*Vs-cf2*Vs
  
  return(colSums(grad2))
}

hessian_alpha_beta <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz)
{
  # ord_id=table(cellidx)
  # nsub = length(ord_id)
  # w_l=list()
  # for(i in 1:length(ord_id)){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  
  alpha =  result$alpha
  theta = result$theta
  sigma2 = result$sigma2
  gamma = result$gamma
  beta = result$beta
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  
  #subject_id = unique(cellidx)
  alpha0 = alpha[1:K]
  alpha1 = alpha[(K+1):(2*K)]
  alpha2 = alpha[-(1:(2*K))]
  alpha.matrix <- rbind(alpha0,alpha1)
  
  covariate_term = ifelse(!is.null(Ub),Ub%*%alpha2,0)
  T0 = 1/beta*exp(log(Lb)+ covariate_term)
  
  Xb = model.matrix(~ytrt)
  C = rowSums(WB*exp(Xb%*%alpha.matrix))
  D = cbind(WB,WB*ytrt)*cbind(exp(Xb%*%alpha.matrix),exp(Xb%*%alpha.matrix))
  
  
  cb1 = -as.vector(T0/(1+T0*C))*D
  
  cb2 = as.vector((B+beta)/beta*T0/(1+T0*C)^2)*D
  grad2 = cb1+cb2
  
  return(colSums(grad2))
}
hessian_alpha_mubj_ind <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz)
{
  # ord_id=table(cellidx)
  # nsub = length(ord_id)
  # w_l=list()
  # for(i in 1:length(ord_id)){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  
  beta = result$beta
  sigma2 = result$sigma2
  gamma = result$gamma
  theta = result$theta
  alpha = result$alpha
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  
  #subject_id = unique(cellidx)
  alpha0 = alpha[1:K]
  alpha1 = alpha[(K+1):(2*K)]
  alpha2 = alpha[-(1:(2*K))]
  Vs = Xs[,1:(2*K)]
  Wstrt = Xs[,(K+1):(2*K)]
  Us = Xs[,-(1:(2*K))]
  alpha.matrix <- rbind(alpha0,alpha1)
  
  A0 = -Xs%*%alpha-mubj+sigmabj2/2-log(Ls)
  G0 = Xs%*%alpha+mubj+sigmabj2/2+log(Ls)
  
  A = theta*exp(A0)
  G = 1/theta*exp(G0)
  
  #grad3 = as.vector(-(1-pi)*(S>0)*((S+theta)*A/(1+A)^2))*Vs - as.vector((1-pi)*(S==0)*exp(G0)/(1+G)^2)*Vs
  
  #tmp = sapply(1:nsub,function(j){colSums(grad3[which(cellidx==j),])})
  tmp1 = -as.vector((1-pi)*(S>0)*(S+theta)*A/(1+A)^2)*Vs-as.vector((1-pi)*(S==0)*exp(G0)/(1+G)^2)*Vs
  tmp2 = colSums(tmp1)
  return(tmp2)
}
hessian_alpha_sigmabj2_ind <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz)
{
  # ord_id=table(cellidx)
  # nsub = length(ord_id)
  # w_l=list()
  # for(i in 1:length(ord_id)){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
 
  beta = result$beta
  sigma2 = result$sigma2
  gamma = result$gamma
  theta = result$theta
  alpha = result$alpha
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  
  #subject_id = unique(cellidx)
  alpha0 = alpha[1:K]
  alpha1 = alpha[(K+1):(2*K)]
  alpha2 = alpha[-(1:(2*K))]
  Vs = Xs[,1:(2*K)]
  Wstrt = Xs[,(K+1):(2*K)]
  Us = Xs[,-(1:(2*K))]
  alpha.matrix <- rbind(alpha0,alpha1)
  
  A0 = -Xs%*%alpha-mubj+sigmabj2/2-log(Ls)
  G0 = Xs%*%alpha+mubj+sigmabj2/2+log(Ls)
  
  A = theta*exp(A0)
  G = 1/theta*exp(G0)
  
  grad3 = 0.5*(as.vector((1-pi)*(S>0)*((S+theta)*A/(1+A)^2))*Vs - as.vector((1-pi)*(S==0)*exp(G0)/(1+G)^2)*Vs)
  
  tmp2 = colSums(grad3)
  return(tmp2)
}
hessian_alpha_pi_ind <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz)
{
  # ord_id=table(cellidx)
  # nsub = length(ord_id)
  # w_l=list()
  # for(i in 1:length(ord_id)){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  
  beta = result$beta
  sigma2 = result$sigma2
  gamma = result$gamma
  theta = result$theta
  alpha = result$alpha
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  
  #subject_id = unique(cellidx)
  alpha0 = alpha[1:K]
  alpha1 = alpha[(K+1):(2*K)]
  alpha2 = alpha[-(1:(2*K))]
  Vs = Xs[,1:(2*K)]
  Wstrt = Xs[,(K+1):(2*K)]
  Us = Xs[,-(1:(2*K))]
  alpha.matrix <- rbind(alpha0,alpha1)
  
  A0 = -Xs%*%alpha-mubj+sigmabj2/2-log(Ls)
  G0 = Xs%*%alpha+mubj+sigmabj2/2+log(Ls)
  
  A = theta*exp(A0)
  G = 1/theta*exp(G0)
  
  tmp = as.vector((S==0)*exp(G0)/(1+G))*Vs
  
  return(t(tmp))
}





################Theta function#########################
hessian_theta_ind <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz){
  
  
  beta = result$beta
  sigma2 = result$sigma2
  gamma = result$gamma
  theta = result$theta
  alpha = result$alpha
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  
  A0 = -Xs%*%alpha-mubj+sigmabj2/2-log(Ls)
  G0 = Xs%*%alpha+mubj+sigmabj2/2+log(Ls)
  
  A = theta*exp(A0)
  G = 1/theta*exp(G0)
  
  hestheta = (1-pi)*(S>0)*(trigamma(S+theta)-trigamma(theta)-
                             2*exp(A0)/(1+A)+
                             (S+theta)*exp(2*A0)/(1+A)^2+
                             1/theta)-
    (1-pi)*(S==0)*(exp(G0)/(theta+exp(G0))^2-G/(theta+exp(G0)))
  return(sum(hestheta))
  
}



hessian_gamma_ind <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz){
  # ord_id=table(cellidx)
  # nsub = length(ord_id)
  # w_l=list()
  # for(i in 1:length(ord_id)){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  ns = length(S)
  
  #beta = result$beta
  theta = result$theta
  sigma2 = result$sigma2
  alpha = result$alpha
  gamma = result$gamma
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  
  tmp = -Reduce("+",lapply(1:ns,function(i){as.vector(exp(Xsz[i,]%*%gamma)/(1+exp(Xsz[i,]%*%gamma))^2)*tcrossprod(Xsz[i,])}))
  return(tmp)
}




hessian_beta <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz){
  # ord_id=table(cellidx)
  # nsub = length(ord_id)
  # w_l=list()
  # for(i in 1:length(ord_id)){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  
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
  alpha.matrix <- rbind(alpha0,alpha1)
  
  
  covariate_term = ifelse(is.null(Ub),0,Ub%*%alpha2)
  
  Xb = model.matrix(~ytrt)
  C = exp(log(Lb)+ covariate_term)*rowSums(WB*exp(Xb%*%alpha.matrix))
  #T0 = 1/beta*C
  
  hesbeta = trigamma(B+beta)-trigamma(beta)+2*C/(beta^2+beta*C)-
    (B+beta)*C*(2*beta+C)/(beta^2+beta*C)^2+B/beta^2
  
  return(sum(hesbeta))
  
}

hessian_sigma2_ind <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz,mum,sm2){
  
  
  beta = result$beta
  theta = result$theta
  sigma2 = result$sigma2
  alpha = result$alpha
  gamma = result$gamma
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  

  
  tmp = 0.5*(1/sigma2^2-2*(mubj^2+sigmabj2)/sigma2^3)#-1/sm2*(1-log(sigma2)+mum)/sigma2^2
  return(sum(tmp))
}


hessian_mubj_ind <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz){

 
  beta = result$beta
  sigma2 = result$sigma2
  gamma = result$gamma
  theta = result$theta
  alpha = result$alpha
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  
  A0 = -Xs%*%alpha-mubj+sigmabj2/2-log(Ls)
  G0 = Xs%*%alpha+mubj+sigmabj2/2+log(Ls)
  
  A = theta*exp(A0)
  G = 1/theta*exp(G0)
  
  # tmp1 = as.vector((1-pi)*(S>0)*((S+theta)*A/(1+A)-theta))*w-
  #   as.vector((1-pi)*(S==0)*exp(G0)/(1+G))*w
  # gradmubj <- -mubj/sigma2 + colSums(tmp1)
  
  tmp1 = -as.vector((1-pi)*(S>0)*(S+theta)*A/(1+A)^2)-as.vector((1-pi)*(S==0)*exp(G0)/(1+G)^2)
  hesmubj = sum(tmp1) - 1/sigma2
  return(hesmubj)
  
}

hessian_sigmabj2_ind <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz){
  # ord_id=table(cellidx)
  # nsub = length(ord_id)
  # w_l=list()
  # for(i in 1:length(ord_id)){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  # sample_in = which(colSums(w)>0)
  # nsub =  length(sample_in)
  beta = result$beta
  sigma2 = result$sigma2
  gamma = result$gamma
  theta = result$theta
  alpha = result$alpha
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  
  A0 = -Xs%*%alpha-mubj+sigmabj2/2-log(Ls)
  G0 = Xs%*%alpha+mubj+sigmabj2/2+log(Ls)
  
  A = theta*exp(A0)
  G = 1/theta*exp(G0)
  
  # tmp1 = as.vector((1-pi)*(S>0)*((S+theta)*A/(1+A)-theta))*w-
  #   as.vector((1-pi)*(S==0)*exp(G0)/(1+G))*w
  # gradmubj <- -mubj/sigma2 + colSums(tmp1)
  
  tmp1 = -1/4*as.vector((1-pi)*(S>0)*(S+theta)*A/(1+A)^2)-1/4*as.vector((1-pi)*(S==0)*exp(G0)/(1+G)^2)
  hessgbj = sum(tmp1)-0.5/sigmabj2^2
  return(hessgbj)
  
}

hessian_pi_ind <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz){
 
  
  beta = result$beta
  sigma2 = result$sigma2
  gamma = result$gamma
  theta = result$theta
  alpha = result$alpha
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  
  tmp = -as.vector((S==0)*(1/pi+1/(1-pi)))
  tmp[is.na(tmp)] = 0
  tmp[tmp == -Inf | tmp == Inf] = 0
  # tmp1 = as.vector((1-pi)*(S>0)*((S+theta)*A/(1+A)-theta))*w-
  #   as.vector((1-pi)*(S==0)*exp(G0)/(1+G))*w
  # gradmubj <- -mubj/sigma2 + colSums(tmp1)
  
  return(diag(tmp))
  
}

hessian_mubj_sigmabj2_ind <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz){
  # ord_id=table(cellidx)
  # nsub = length(ord_id)
  # w_l=list()
  # for(i in 1:length(ord_id)){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  
  beta = result$beta
  sigma2 = result$sigma2
  gamma = result$gamma
  theta = result$theta
  alpha = result$alpha
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  
  A0 = -Xs%*%alpha-mubj+sigmabj2/2-log(Ls)
  G0 = Xs%*%alpha+mubj+sigmabj2/2+log(Ls)
  
  A = theta*exp(A0)
  G = 1/theta*exp(G0)
  
  # tmp1 = as.vector((1-pi)*(S>0)*((S+theta)*A/(1+A)-theta))*w-
  #   as.vector((1-pi)*(S==0)*exp(G0)/(1+G))*w
  # gradmubj <- -mubj/sigma2 + colSums(tmp1)
  
  tmp1 = 0.5*as.vector((1-pi)*(S>0)*(S+theta)*A/(1+A)^2)-0.5*as.vector((1-pi)*(S==0)*exp(G0)/(1+G)^2)
  hesmusgbj = sum(tmp1)
  return(hesmusgbj)
  
}

hessian_mubj_pi_ind <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz){
  # ord_id=table(cellidx)
  # nsub = length(ord_id)
  # w_l=list()
  # for(i in 1:length(ord_id)){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  # sample_in = which(colSums(w)>0)
  # nsub =  length(sample_in)
  
  beta = result$beta
  sigma2 = result$sigma2
  gamma = result$gamma
  theta = result$theta
  alpha = result$alpha
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  
  A0 = -Xs%*%alpha-mubj+sigmabj2/2-log(Ls)
  G0 = Xs%*%alpha+mubj+sigmabj2/2+log(Ls)
  
  A = theta*exp(A0)
  G = 1/theta*exp(G0)
  
  # tmp1 = as.vector((1-pi)*(S>0)*((S+theta)*A/(1+A)-theta))*w-
  #   as.vector((1-pi)*(S==0)*exp(G0)/(1+G))*w
  # gradmubj <- -mubj/sigma2 + colSums(tmp1)
  
  hesmupi = as.vector((S==0)*exp(G0)/(1+G)) #- as.vector((S>0)*((S+theta)*A/(1+A)+theta))
  return(hesmupi)
  
}

#############need attention for pi###########
#############################################

hessian_sigmabj2_pi_ind <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz){
  # ord_id=table(cellidx)
  # nsub = length(ord_id)
  # w_l=list()
  # for(i in 1:length(ord_id)){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
 
  
  beta = result$beta
  sigma2 = result$sigma2
  gamma = result$gamma
  theta = result$theta
  alpha = result$alpha
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  
  A0 = -Xs%*%alpha-mubj+sigmabj2/2-log(Ls)
  G0 = Xs%*%alpha+mubj+sigmabj2/2+log(Ls)
  
  A = theta*exp(A0)
  G = 1/theta*exp(G0)
  
  # tmp1 = as.vector((1-pi)*(S>0)*((S+theta)*A/(1+A)-theta))*w-
  #   as.vector((1-pi)*(S==0)*exp(G0)/(1+G))*w
  # gradmubj <- -mubj/sigma2 + colSums(tmp1)
  
  tmp1 = 0.5*as.vector((S==0)*exp(G0)/(1+G)) #- as.vector((S>0)*((S+theta)*A/(1+A)+theta))
  hesmupi = tmp1
  return(hesmupi)
  
}

hessian_theta_mubj_ind <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz){
  # ord_id=table(cellidx)
  # nsub = length(ord_id)
  # w_l=list()
  # for(i in 1:length(ord_id)){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
 
  beta = result$beta
  sigma2 = result$sigma2
  gamma = result$gamma
  theta = result$theta
  alpha = result$alpha
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  
  A0 = -Xs%*%alpha-mubj+sigmabj2/2-log(Ls)
  G0 = Xs%*%alpha+mubj+sigmabj2/2+log(Ls)
  
  A = theta*exp(A0)
  G = 1/theta*exp(G0)
  
  hesthetamu = (1-pi)*(S>0)*(exp(A0)*S-1)/(1+A)^2-
    (1-pi)*(S==0)*exp(2*G0)/(theta+exp(G0))^2
  tmp0 = sum(hesthetamu)
  return(tmp0)
  
}

hessian_theta_sigmabj2_ind <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz){
  # ord_id=table(cellidx)
  # nsub = length(ord_id)
  # w_l=list()
  # for(i in 1:length(ord_id)){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
 
  
  beta = result$beta
  sigma2 = result$sigma2
  gamma = result$gamma
  theta = result$theta
  alpha = result$alpha
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  
  A0 = -Xs%*%alpha-mubj+sigmabj2/2-log(Ls)
  G0 = Xs%*%alpha+mubj+sigmabj2/2+log(Ls)
  
  A = theta*exp(A0)
  G = 1/theta*exp(G0)
  
  hesthetas = -0.5*(1-pi)*(S>0)*(1+(exp(A0)*S-1)/(1+A)^2)-
    0.5*(1-pi)*(S==0)*exp(2*G0)/(theta+exp(G0))^2
  tmp0 = sum(hesthetas)
  return(tmp0)
  
}

hessian_theta_pi_ind<- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz){
  # ord_id=table(cellidx)
  # nsub = length(ord_id)
  # w_l=list()
  # for(i in 1:length(ord_id)){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
 
  beta = result$beta
  sigma2 = result$sigma2
  gamma = result$gamma
  theta = result$theta
  alpha = result$alpha
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  
  A0 = -Xs%*%alpha-mubj+sigmabj2/2-log(Ls)
  G0 = Xs%*%alpha+mubj+sigmabj2/2+log(Ls)
  
  A = theta*exp(A0)
  G = 1/theta*exp(G0)
  
  hesthetapi = as.vector((S==0)*(log(1+G)-G/(1+G)))
  
  return(hesthetapi)
  
}

hessian_sigma2_mubj_ind<- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz){
  # ord_id=table(cellidx)
  # nsub = length(ord_id)
  # w_l=list()
  # for(i in 1:length(ord_id)){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
 
  
  beta = result$beta
  sigma2 = result$sigma2
  gamma = result$gamma
  theta = result$theta
  alpha = result$alpha
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  
  
  tmp0 = mubj/sigma2^2
  tmp0
}

hessian_sigma2_sigmabj2_ind<- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz){
  # ord_id=table(cellidx)
  # nsub = length(ord_id)
  # w_l=list()
  # for(i in 1:length(ord_id)){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
 
  beta = result$beta
  sigma2 = result$sigma2
  gamma = result$gamma
  theta = result$theta
  alpha = result$alpha
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  
 
  tmp0 = 1/2/sigma2^2
  tmp0
}

hessian_gamma_pi_ind <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz){
  mat = Xsz
  mat[which(S>0),] = 0
  t(mat)
}
##############gradient##############
####################################
##########################################
#########################################

grad_beta_final <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz){
  # ord_id=table(cellidx)
  # nsub = length(ord_id)
  # w_l=list()
  # for(i in 1:length(ord_id)){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  
  theta = result$theta
  sigma2 = result$sigma2
  gamma = result$gamma
  beta = result$beta
  alpha = result$alpha
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  
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
  tmp2 =  digamma_diff - log_term + (B + beta) * grad_log_term  - B / beta
  
  Fbeta = sum(tmp2)
  return(Fbeta)
}

grad_theta_final <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz){
  
  # ord_id=table(cellidx)
  # nsub = length(ord_id)
  # w_l=list()
  # for(i in 1:length(ord_id)){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  
  theta = result$theta
  sigma2 = result$sigma2
  gamma = result$gamma
  beta = result$beta
  alpha = result$alpha
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  
  
  dgamma_diff <- digamma(S+theta)-digamma(theta)
  term1 = -Xs%*%alpha-mubj+sigmabj2/2-log(Ls)
  term2 = Xs%*%alpha+mubj+sigmabj2/2+log(Ls)
  
  logmu = Xs%*%alpha+mubj+log(Ls)
  common_term1 <- theta*exp(term1)
  common_term2 <- 1/theta*exp(term2)
  
  Ftheta = (1-pi)*(S>0)*(dgamma_diff-
                           log(1+common_term1)-
                           (S+theta)*exp(term1)/(1+common_term1)+
                           1+log(theta)-logmu)+(1-pi)*(S==0)*(exp(term2)/(theta+exp(term2))-log(1+common_term2))
  Ftheta0 = sum(Ftheta)
  return(Ftheta0)
}

grad_gamma_final <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz){
  # ord_id=table(cellidx)
  # nsub = length(ord_id)
  # w_l=list()
  # for(i in 1:length(ord_id)){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  ns = length(S)
  
  #beta = result$beta
  theta = result$theta
  sigma2 = result$sigma2
  alpha = result$alpha
  gamma = result$gamma
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  
  
  tmp = as.vector(pi)*Xsz - as.vector(exp(Xsz%*%gamma)/(1+exp(Xsz%*%gamma)))*Xsz 
  
  return(colSums(tmp))
}

grad_alpha_final <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz){
  # ord_id=table(cellidx)
  # nsub = length(ord_id)
  # w_l=list()
  # for(i in 1:length(ord_id)){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  
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
  Ws = Xs[,1:K]
  Wstrt = Xs[,(K+1):(2*K)]
  Us = Xs[,-(1:(2*K))]
  
  term1 = -Xs%*%alpha-mubj+sigmabj2/2-log(Ls)
  term2 = Xs%*%alpha+mubj+sigmabj2/2+log(Ls)
  
  
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
  
  return(gradalpha)
}

grad_alpha_bulk_final <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz,lambda = 0.1){
  # ord_id=table(cellidx)
  # nsub = length(ord_id)
  # w_l=list()
  # for(i in 1:length(ord_id)){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  
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
    gradalpha <- c(colSums(tmp20),colSums(tmp21)) ####to be verified the dimension
    
  }else{
    tmp22 = -(B+beta)/beta*exp(term.bulk)*Ub/(1+1/beta*exp(term.bulk))+B*Ub
    gradalpha <- c(colSums(tmp20),colSums(tmp21),colSums(tmp22)) ####to be verified the dimension
    
  }
  
  return(gradalpha)
}

grad_alpha_sc_final <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz){
  # ord_id=table(cellidx)
  # nsub = length(ord_id)
  # w_l=list()
  # for(i in 1:length(ord_id)){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  
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
  Vs = Xs[,1:(2*K)]
  #Wstrt = Xs[,(K+1):(2*K)]
  Us = Xs[,-(1:(2*K))]
  
  term1 = -Xs%*%alpha-mubj+sigmabj2/2-log(Ls)
  term2 = Xs%*%alpha+mubj+sigmabj2/2+log(Ls)
  
  
  tmp10 = as.vector((1-pi)*(S>0)*(S+theta)*theta*exp(term1))*Vs/c(1+theta*exp(term1))-as.vector((1-pi)*(S>0)*theta)*Vs-
    as.vector((1-pi)*(S==0)*exp(term2))*Vs/c(1+1/theta*exp(term2))
  
  
  gradalpha <- colSums(tmp10) ####to be verified the dimension
  
  
  
  return(gradalpha)
}


grad_sigma2_final <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz,mum,sm2){
  
  
  
  beta = result$beta
  theta = result$theta
  sigma2 = result$sigma2
  alpha = result$alpha
  gamma = result$gamma
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  
  
  tmp = 0.5*(-1/(sigma2)+(mubj^2+sigmabj2)/sigma2^2)#-(log(sigma2)-mum)/sm2/sigma2
  return(sum(tmp))
}

gradient_sc_total_ind <-function(result,S,B,Xs,WB,Ub,ysc,
                             ytrt,Lb,Ls,K,Xsz,mum,sm2){
  # nsub = length(unique(celltotal))
  # ord_id = sapply(1:nsub,function(x){sum(cellidx == x)})
  # w_l=list()
  # for(i in 1:nsub){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  
  x1 <- 0
  x2 <- grad_theta_final(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  x3 <- grad_sigma2_final(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz,mum,sm2)
  x4 <- grad_gamma_final(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  x5 <- grad_alpha_sc_final(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  return(c(x1,x2,x3,x4,x5))
}

gradient_bulk_total <-function(result = result0[[i]],
                               S = Y[i,],B= B0[i,],Xs = modmat,WB = W,Ub=NULL,ysc= trt,
                               ytrt = trt.bulk,
                               Lb = libsizeb,
                               Ls =libsizes,K=K,Xsz = Xsz,celltotal = cellidx0){
  # nsub = length(unique(cellidx0))
  # ord_id = sapply(1:nsub,function(x){sum(cellidx == x)})
  # w_l=list()
  # for(i in 1:nsub){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  # 
  x1 <- grad_beta_final(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  x2 <- 0
  x3 <- 0
  x4 <- rep(0,ncol(Xsz))
  x5 <- grad_alpha_bulk_final(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  return(c(x1,x2,x3,x4,x5))
}

########################################
########### Bulk separate ###############
########################################

grad_alpha_bulk_individual <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz,lambda = 0.1){
  # ord_id=table(cellidx)
  # nsub = length(ord_id)
  # w_l=list()
  # for(i in 1:length(ord_id)){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  
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
    gradalpha <- cbind(tmp20,tmp21) ####to be verified the dimension -2*lambda*colSums(WB)*alpha1
    
  }else{
    tmp22 = -(B+beta)/beta*exp(term.bulk)*Ub/(1+1/beta*exp(term.bulk))+B*Ub
    gradalpha <- cbind(tmp20,tmp21,tmp22) ####to be verified the dimension
    
  }
  
  return(gradalpha)
}

grad_beta_individual <- function(result,S,B,Xs,WB,Ub=NULL,ysc,ytrt,Lb,Ls,K,Xsz){
  # ord_id=table(cellidx)
  # nsub = length(ord_id)
  # w_l=list()
  # for(i in 1:length(ord_id)){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  
  theta = result$theta
  sigma2 = result$sigma2
  gamma = result$gamma
  beta = result$beta
  alpha = result$alpha
  
  mubj = result$mubj
  sigmabj2 = result$sigmabj2
  pi = result$pi
  
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
  tmp2 =  digamma_diff - log_term + (B + beta) * grad_log_term  - B / beta
  
  
  return(tmp2)
}

gradient_bulk_individual <- function(result = result0[[i]],
                                     S = Y[i,],B= B0[i,],Xs = modmat,WB = W,Ub=NULL,ysc= trt,
                                     ytrt = trt.bulk,
                                     Lb = libsizeb,
                                     Ls =libsizes,K=K,Xsz = Xsz,celltotal = cellidx0,lambda = 0.1){
  # nsub = length(unique(celltotal))
  # ord_id = sapply(1:nsub,function(x){sum(cellidx == x)})
  # w_l=list()
  # for(i in 1:nsub){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  
  nlen = length(B)
  
  x1 <- grad_beta_individual(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  x2 <- rep(0,nlen)
  x3 <- rep(0,nlen)
  x4 <- matrix(0,nrow = nlen, ncol = ncol(Xsz))
  x5 <- grad_alpha_bulk_individual(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz,lambda)
  return(cbind(x1,x2,x3,x4,x5))
}

gradient_total <- function(result = result0[[i]],
                           S = Y[i,],B= B0[i,],Xs = modmat,WB = W,Ub=NULL,ysc= trt,
                           ytrt = trt.bulk,cellidx = cellidx,
                           Lb = libsizeb,
                           Ls =libsizes,K=K,Xsz = Xsz,celltotal = cellidx0){
  nsub = length(unique(celltotal))
  ord_id = sapply(1:nsub,function(x){sum(cellidx == x)})
  w_l=list()
  for(i in 1:nsub){
    w_l[[i]]=matrix(1,ord_id[i],1)
  }
  w=as.matrix(Matrix::bdiag(w_l))
  
  x1 <- grad_beta_final(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  x2 <- grad_theta_final(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  x3 <- grad_sigma2_final(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  x4 <- grad_gamma_final(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  x5 <- grad_alpha_final(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  return(c(x1,x2,x3,x4,x5))
}

# gradient_total(result = result_part,
#                 S = Y[i,index],B= B0[i,],Xs = modmat[index,],WB = W,Ub=NULL,ysc= trt[index],
#                 ytrt = trt.bulk,cellidx = cellidx,
#                 Lb = Lb,
#                 Ls =Ls[index],K=K,Xsz = Xsz0[index,],celltotal = cellidx0)

# 
# h22 <- hessian(ELBO_partj,x =beta,param = param,psi=psi,S_ind = Y[i,index],
#                 B = B0[i,],Xs_ind = modmat[index,],Xsz_ind=Xsz0[index,],
#                 WB = WB,ytrt = trt.bulk,ysc_ind = trt[index],Ub=NULL,
#                 Ls_ind = Ls[index],Lb = Lb,K = K)
# 
# grad <- grad(ELBO_partj,x =param,psi = psi,S_ind = Y[i,index],
#              B = B0[i,],Xs_ind = modmat[index,],Xsz_ind=Xsz0[index,],
#              WB = WB,ytrt = trt.bulk,ysc_ind = trt[index],Ub=NULL,
#              Ls_ind = Ls[index],Lb = Lb,K = K)
# grad0 <- grad(ELBO_theta_adjusted,x = result_part$theta,result=result_part,
#               S = Y[i,index],
#               Xs = modmat[index,],w=w,
#              Ls=Ls,K = K)


ELBO_total_check <- function(param,psi,S,B,Xs,Xsz,WB,ytrt,ysc,Ub=NULL,cellidx,Ls,Lb,K,celltotal){
  
  nsub = length(unique(celltotal))
  ord_id = sapply(1:nsub,function(x){sum(cellidx == x)})
  w_l=list()
  for(i in 1:nsub){
    w_l[[i]]=matrix(1,ord_id[i],1)
  }
  w=as.matrix(Matrix::bdiag(w_l))
  
  
  beta = param[1]
  theta = param[2]
  sigma2 = param[3]
  gamma = param[4:(3+ncol(Xsz))]
  alpha0 = param[(4+ncol(Xsz)):(3+K+ncol(Xsz))]
  alpha1 = param[(4+K+ncol(Xsz)):(3+2*K+ncol(Xsz))]
  if(!is.null(Ub)){
    alpha2 = param[(4+2*K+ncol(Xsz)):(3+2*K+ncol(Xsz)+ncol(Ub))]
    alpha = c(alpha0,alpha1,alpha2)
  }else{
    alpha = c(alpha0,alpha1)
  }
  
  mubj = psi[1:nsub]
  sigmabj2 = psi[(nsub+1):(2*nsub)]
  pi = psi[(2*nsub+1):(2*nsub+length(S))]
  
  
  #####For zero inflation part#############
  ##########################################
  
  Fspart = 0.5*sum((log(sigmabj2)-log(sigma2)-(mubj^2+sigmabj2)/sigma2))
  
  
  Fzero = sum(pi*Xsz%*%gamma - log(1+exp(Xsz%*%gamma)))
  
  
  Fpi = -sum(pi[pi>0]*log(pi[pi>0])) -sum((1-pi[pi<1])*log(1-pi[pi<1]))
  
  Xb = model.matrix(~ytrt)
  alpha.matrix <- rbind(alpha0,alpha1)
  covariate_term = ifelse(!is.null(Ub),Ub%*%alpha2,0)
  term.bulk = log(rowSums(WB*exp(Xb%*%alpha.matrix)))+log(Lb)+ covariate_term
  
  term1 = -Xs%*%alpha-mubj+sigmabj2/2-log(Ls)
  term2 = Xs%*%alpha+mubj+sigmabj2/2+log(Ls)
  
  tmp2 = lgamma(B+beta)-lgamma(beta)-lgamma(B+1)-(B+beta)*log(1+1/beta*exp(term.bulk))+B*term.bulk - B*log(beta)
  
  
  
  lgamma_diff <- lgamma(S+theta)-lgamma(theta)-lgamma(S+1)
  
  
  logmu =  Xs%*%alpha+mubj+log(Ls)
  common_term1 <- theta*exp(term1)
  common_term2 <- 1/theta*exp(term2)
  
  tmp1 = (1-pi)*(S>0)*(lgamma_diff - (S+theta)*log(1+common_term1)+
                         theta*log(theta)-theta*logmu)-(1-pi)*(S==0)*theta*log(1+common_term2)
  
  Ftotal <- sum(tmp1) + sum(tmp2) + Fspart + Fzero + Fpi
  return(Ftotal)
}

ELBO_bulk_check <- function(param,psi,S,B,Xs,Xsz,WB,ytrt,ysc,Ub=NULL,cellidx,Ls,Lb,K,celltotal){
  
  nsub = length(unique(celltotal))
  ord_id = sapply(1:nsub,function(x){sum(cellidx == x)})
  w_l=list()
  for(i in 1:nsub){
    w_l[[i]]=matrix(1,ord_id[i],1)
  }
  w=as.matrix(Matrix::bdiag(w_l))
  
  
  beta = param[1]
  theta = param[2]
  sigma2 = param[3]
  gamma = param[4:(3+ncol(Xsz))]
  alpha0 = param[(4+ncol(Xsz)):(3+K+ncol(Xsz))]
  alpha1 = param[(4+K+ncol(Xsz)):(3+2*K+ncol(Xsz))]
  if(!is.null(Ub)){
    alpha2 = param[(4+2*K+ncol(Xsz)):(3+2*K+ncol(Xsz)+ncol(Ub))]
    alpha = c(alpha0,alpha1,alpha2)
  }else{
    alpha = c(alpha0,alpha1)
  }
  
  mubj = psi[1:nsub]
  sigmabj2 = psi[(nsub+1):(2*nsub)]
  pi = psi[(2*nsub+1):(2*nsub+length(S))]
  
  
  #####For zero inflation part#############
  ##########################################
  
  
  Xb = model.matrix(~ytrt)
  alpha.matrix <- rbind(alpha0,alpha1)
  covariate_term = ifelse(!is.null(Ub),Ub%*%alpha2,0)
  term.bulk = log(rowSums(WB*exp(Xb%*%alpha.matrix)))+log(Lb)+ covariate_term
  
  
  
  tmp2 = lgamma(B+beta)-lgamma(beta)-lgamma(B+1)-(B+beta)*log(1+1/beta*exp(term.bulk))+B*term.bulk - B*log(beta)
  
  
  
  Ftotal <- sum(tmp2)
  return(Ftotal)
}

ELBO_sc_check <- function(param,psi,pi,S,B,Xs,Xsz,WB,ytrt,ysc,Ub=NULL,cellidx,Ls,Lb,K,celltotal){
  
  # nsub = length(unique(celltotal))
  # ord_id = sapply(1:nsub,function(x){sum(cellidx == x)})
  # w_l=list()
  # for(i in 1:nsub){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  
  
  beta = param[1]
  theta = param[2]
  sigma2 = param[3]
  gamma = param[4:(3+ncol(Xsz))]
  alpha0 = param[(4+ncol(Xsz)):(3+K+ncol(Xsz))]
  alpha1 = param[(4+K+ncol(Xsz)):(3+2*K+ncol(Xsz))]
  if(!is.null(Ub)){
    alpha2 = param[(4+2*K+ncol(Xsz)):(3+2*K+ncol(Xsz)+ncol(Ub))]
    alpha = c(alpha0,alpha1,alpha2)
  }else{
    alpha = c(alpha0,alpha1)
  }
  
  ###va parameter for subject j###
  mubj = psi[1]
  sigmabj2 = psi[2]
  pi = psi[(3):(2+length(S))]
  
  
  #####For zero inflation part#############
  ##########################################
  
  # Fspart = 0.5*sum((log(sigmabj2[unique(cellidx)])-log(sigma2)-(mubj[unique(cellidx)]^2+sigmabj2[unique(cellidx)])/sigma2))
  
  Fspart = 0.5*(log(sigmabj2)-log(sigma2)-(mubj^2+sigmabj2)/sigma2)
  
  Fzero = sum(pi*Xsz%*%gamma - log(1+exp(Xsz%*%gamma)))
  
  
  Fpi = -sum(pi[pi>0]*log(pi[pi>0])) -sum((1-pi[pi<1])*log(1-pi[pi<1]))
  
  
  term1 = -Xs%*%alpha-mubj+sigmabj2/2-log(Ls)
  term2 = Xs%*%alpha+mubj+sigmabj2/2+log(Ls)
  
  
  
  
  lgamma_diff <- lgamma(S+theta)-lgamma(theta)-lgamma(S+1)
  
  
  logmu =  Xs%*%alpha+mubj+log(Ls)
  common_term1 <- theta*exp(term1)
  common_term2 <- 1/theta*exp(term2)
  
  tmp1 = (1-pi)*(S>0)*(lgamma_diff - (S+theta)*log(1+common_term1)+theta*log(theta)-theta*logmu)-(1-pi)*(S==0)*theta*log(1+common_term2)
  
  Ftotal <- sum(tmp1) + Fspart + Fzero + Fpi
  return(Ftotal)
}

Hess_modelparam <- function(result = result0[[i]],
                            S = Y[i,],B= B0[i,],Xs = modmat,WB = W,Ub=NULL,ysc= trt,
                            ytrt = trt.bulk,cellidx = cellidx,
                            Lb = libsizeb,
                            Ls =libsizes,K=K,Xsz = Xsz,celltotal = cellidx0){
  
  nsub = length(unique(celltotal))
  # ord_id = sapply(1:nsub,function(x){sum(cellidx == x)})
  # w_l=list()
  # for(i in 1:nsub){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  ncov = ifelse(is.null(Ub),0,ncol(Ub))
  Hessp = matrix(0,nrow = 2*K+ncov+ncol(Xsz)+3, ncol = 2*K+ncov+ncol(Xsz)+3)
  #2nd derivative of beta, theta,sigma2,gamma,alpha
  Hessp[1,1] = hessian_beta(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  Hessp[2,2] = hessian_theta(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  Hessp[3,3] = hessian_sigma2(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  Hessp[4:(3+ncol(Xsz)),4:(3+ncol(Xsz))] = hessian_gamma(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  Hessp[(4+ncol(Xsz)):nrow(Hessp),(4+ncol(Xsz)):nrow(Hessp)] = 
    hessian_alpha(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  
  ###beta * alpha
  Hessp[1,(1+3+ncol(Xsz)):ncol(Hessp)] = Hessp[(1+3+ncol(Xsz)):ncol(Hessp),1] = 
    hessian_alpha_beta(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  ###theta* alpha
  Hessp[2,(1+3+ncol(Xsz)):ncol(Hessp)] = Hessp[(1+3+ncol(Xsz)):ncol(Hessp),2]=
    hessian_alpha_theta(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  return(Hessp)
  ###theta* alpha
}

Hess_modelparam_bulk <- function(result = result0[[i]],
                                 S = Y[i,],B= B0[i,],Xs = modmat,WB = W,Ub=NULL,ysc= trt,
                                 ytrt = trt.bulk,
                                 Lb = libsizeb,
                                 Ls =libsizes,K=K,Xsz = Xsz,celltotal = cellidx0,lambda = 0.1){
  
  nsub = length(unique(celltotal))
  # ord_id = sapply(1:nsub,function(x){sum(cellidx == x)})
  # w_l=list()
  # for(i in 1:nsub){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  ncov = ifelse(is.null(Ub),0,ncol(Ub))
  Hessp = matrix(0,nrow = 2*K+ncov+ncol(Xsz)+3, ncol = 2*K+ncov+ncol(Xsz)+3)
  #2nd derivative of beta, theta,sigma2,gamma,alpha
  Hessp[1,1] = hessian_beta(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  Hessp[(4+ncol(Xsz)):nrow(Hessp),(4+ncol(Xsz)):nrow(Hessp)] = 
    hessian_alpha_bulk(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz,lambda)
  
  ###beta * alpha
  Hessp[1,(1+3+ncol(Xsz)):ncol(Hessp)] = Hessp[(1+3+ncol(Xsz)):ncol(Hessp),1] = 
    hessian_alpha_beta(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  
  return(Hessp)
  ###theta* alpha
}

Hess_modelparam_sc_ind <- function(result,S,B,Xs,WB,Ub,ysc,
                               ytrt,Lb,Ls,K,Xsz,celltotal,mum,sm2){
  
  nsub = length(unique(celltotal))
  # ord_id = sapply(1:nsub,function(x){sum(cellidx == x)})
  # w_l=list()
  # for(i in 1:nsub){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  ncov = ifelse(is.null(Ub),0,ncol(Ub))
  Hessp = matrix(0,nrow = 2*K+ncov+ncol(Xsz)+3, ncol = 2*K+ncov+ncol(Xsz)+3)
  #2nd derivative of beta, theta,sigma2,gamma,alpha
  Hessp[2,2] = hessian_theta_ind(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  Hessp[3,3] = hessian_sigma2_ind(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz,mum,sm2)
  Hessp[4:(3+ncol(Xsz)),4:(3+ncol(Xsz))] = hessian_gamma_ind(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  Hessp[(4+ncol(Xsz)):nrow(Hessp),(4+ncol(Xsz)):nrow(Hessp)] = 
    hessian_alpha_sc_ind(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  
  ###theta* alpha
  Hessp[2,(1+3+ncol(Xsz)):ncol(Hessp)] = Hessp[(1+3+ncol(Xsz)):ncol(Hessp),2]=
    hessian_alpha_theta_ind(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  return(Hessp)
  ###theta* alpha
}

Hess_variationalparam_ind <- function(result = result0[[i]],
                                  S = Y[i,],B= B0[i,],Xs = modmat,WB = W,Ub=NULL,ysc= trt,
                                  ytrt = trt.bulk,
                                  Lb = libsizeb,
                                  Ls =libsizes,K=K,Xsz = Xsz,celltotal = cellidx0){
  ##########mubj, sigmabj2, pi###############
  dimrow = 2+length(S)
  Hessv <- matrix(0,nrow = dimrow, ncol = dimrow)
  muidx = 1
  sgidx = 2
  piidx = 3:(dimrow)
  nsub0 = length(unique(celltotal))
  # ord_id = sapply(1:nsub0,function(x){sum(cellidx == x)})
  # w_l=list()
  # for(i in 1:nsub0){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  
  Hessv[muidx,muidx] <- 
    hessian_mubj_ind(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  Hessv[sgidx,sgidx]<-
    hessian_sigmabj2_ind(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  Hessv[piidx,piidx]<-
    hessian_pi_ind(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  
  ############mubj vs sigmabj################
  ##########################################
  
  Hessv[muidx,sgidx] = hessian_mubj_sigmabj2_ind(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  Hessv[sgidx,muidx]= t(Hessv[muidx,sgidx])
  ######Mubj vs pi ###################
  ####################################
  ####################################
  
  Hessv[piidx,muidx] =hessian_mubj_pi_ind(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  Hessv[muidx,piidx]= t(Hessv[piidx,muidx])
  ######sigmabj vs pi ###################
  ####################################
  ####################################     
  Hessv[piidx,sgidx] =hessian_sigmabj2_pi_ind(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  Hessv[sgidx,piidx]= t(Hessv[piidx,sgidx])
  return(Hessv)
}

Hess_varmodparam_ind<- function(result = result0[[i]],
                            S = Y[i,],B= B0[i,],Xs = modmat,WB = W,Ub=NULL,ysc= trt,
                            ytrt = trt.bulk,
                            Lb = libsizeb,
                            Ls =libsizes,K=K,Xsz = Xsz,celltotal = cellidx0){
  ######beta, theta,sigma2,gamma,alpha ###
  #####mubj, sigma2bj, pi ###############
  # nsub = length(unique(celltotal))
  # ord_id = sapply(1:nsub,function(x){sum(cellidx == x)})
  # w_l=list()
  # for(i in 1:nsub){
  #   w_l[[i]]=matrix(1,ord_id[i],1)
  # }
  # w=as.matrix(Matrix::bdiag(w_l))
  ncov = ifelse(is.null(Ub),0,ncol(Ub))
  param_ln = 2*K+ncov+ncol(Xsz)+3
  psi_ln = 2+length(S)
  muidx = 1
  sgidx = 2
  piidx = 3:(psi_ln)
  
  Hessmatcor = matrix(0,nrow = param_ln, ncol = psi_ln)
  
  ############theta vs mubj ######################
  Hessmatcor[2,muidx] = hessian_theta_mubj_ind(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  ############theta vs sigmabj ######################
  Hessmatcor[2,sgidx] = hessian_theta_sigmabj2_ind(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  ############theta vs pi ######################
  Hessmatcor[2,piidx] = hessian_theta_pi_ind(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  
  #########row 3#################
  ############sigma2 vs mubj ######################
  Hessmatcor[3,muidx] = hessian_sigma2_mubj_ind(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  ############sigma2 vs sigma2bj ######################
  Hessmatcor[3,sgidx] = hessian_sigma2_sigmabj2_ind(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  
  
  ############gamma vs pi ######################
  Hessmatcor[(4:(3+ncol(Xsz))),piidx] = hessian_gamma_pi_ind(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  
  
  
  ########alpha vs mubj###############
  Hessmatcor[(1+3+ncol(Xsz)):param_ln,muidx] = 
    hessian_alpha_mubj_ind(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  
  ########alpha vs sigmabj###############
  Hessmatcor[(1+3+ncol(Xsz)):param_ln,sgidx] = 
    hessian_alpha_sigmabj2_ind(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  
  ###########alpha vs pi ##################
  Hessmatcor[(1+3+ncol(Xsz)):param_ln,piidx] = 
    hessian_alpha_pi_ind(result,S,B,Xs,WB,Ub,ysc,ytrt,Lb,Ls,K,Xsz)
  
  return(Hessmatcor)
}

VEM_cov<- function(result,a1=1,a2=1,
                   S,B,Xs,WB,Ub=NULL,ysc,
                   ytrt,cellidx,Lb,mum,sm2,
                   Ls,K,Xsz,celltotal,lambda){
  nsub = length(unique(celltotal))
  nb =length(B)
  
  Lmat <- lapply(1:nsub,function(t){
    index = which(cellidx == t)
    # k=i-900
    result_part = result
    result_part$mubj = result_part$mubj[t]
    result_part$sigmabj2 = result_part$sigmabj2[t]
    result_part$pi = result_part$pi[index]
    Hesspt = Hess_modelparam_sc_ind(result_part,S[index],B,Xs[index,],WB,Ub,ysc[index],
                                    ytrt,Lb,
                                    Ls[index],K,Xsz[index,],celltotal,mum,sm2)
    Hessv = Hess_variationalparam_ind(result_part,S[index],B,Xs[index,],WB,Ub,ysc[index],
                                      ytrt,Lb,
                                      Ls[index],K,Xsz[index,],celltotal)
    Hesscov = Hess_varmodparam_ind(result_part,S[index],B,Xs[index,],WB,Ub,ysc[index],
                                   ytrt,Lb,
                                   Ls[index],K,Xsz[index,],celltotal)
    if(sum(S[index]>0)==0)
    {Hessv0 = Hessv
    Hesscov0 <- Hesscov}else{
      Hessv0 <- Hessv[-(which(colSums(Hessv>0)==0)),-(which(colSums(Hessv>0)==0))]
      Hesscov0 <- Hesscov[,-which(colSums(Hessv>0)==0)]
    }
    
    #print(diag(Hessv0))
    #A_cov = Hesspt - Hesscov%*%ginv(Hessv)%*%t(Hesscov)
    A_cov0 = Hesspt - Hesscov0%*%solve(Hessv0)%*%t(Hesscov0)
    
    vgr = gradient_sc_total_ind(result_part,S[index],B,Xs[index,],WB,Ub,ysc[index],
                                ytrt,Lb,
                                Ls[index],K,Xsz[index,],mum,sm2)
    
    
    B_cov = tcrossprod(vgr)
    list(A_cov = A_cov0, B_cov = B_cov, gr_sc = vgr)
  })
  
  #sc_gr_mean <- 1/nsub*Reduce("+",lapply(Lmat,function(x){x$gr_sc}))
  B1_final = Reduce("+",lapply(Lmat,function(x){x$B_cov}))
  
  A1_final <- Reduce("+",lapply(Lmat,function(x){x$A_cov}))
  
  A2_final <- Hess_modelparam_bulk(result,S,B,Xs,WB,Ub,ysc,
                                   ytrt,Lb,
                                   Ls,K,Xsz,celltotal,lambda)
  
  A_final <- (a1*A1_final+a2*A2_final)/(a1*nsub+a2*nb)
  
  #covt = -solve(A_final)
  #B1_final <- Reduce("+",lapply(Lmat,function(x){x$B_cov}))
  
  gr_bulk=gradient_bulk_individual(result,S,B,Xs,WB,Ub,ysc,
                                   ytrt,Lb,
                                   Ls,K,Xsz,celltotal,lambda)
  # bu_gr = colMeans(gr_bulk)
  
  #B2_final <- sc_gr%*%t(bu_gr)
  #B3_final <- bu_gr%*%t(sc_gr)
  
  B2_final <- Reduce("+",lapply(1:nrow(gr_bulk),function(x){tcrossprod(gr_bulk[x,])}))
  B_final <- (a1^2*B1_final + a2^2*B2_final)/(a1*nsub+a2*nb)
  
  #c(1,6:11),c(1,6:11)
  A_f0 <- A_final
  B_f0 <- B_final
  covt = ginv(A_f0)%*%B_f0%*%ginv(A_f0)/(a2*nb+a1*nsub)
  library(Matrix)
  pvalue <- sapply((K+1):(2*K),function(z){
    C = matrix(rep(0,2*K),nrow=1)
    C[,z]=1
    DE_test(covt[6:(5+2*K),6:(5+2*K)],result$alpha,C)
  })
  
  # pvalue <- sapply(1:2,function(z){
  #   C = matrix(rep(0,2),nrow=1)
  #   C[,z]=1
  #   DE_test(covt[4:5,4:5],result$gamma,C)
  # })
  
  var3 = diag(covt)
  
  return(list(var=var3,pvalue =pvalue))
  
}


VEM_cov_bulk <- function(result,a1=1,a2=1,
                   S,B,Xs,WB,Ub=NULL,ysc,
                   ytrt,cellidx,Lb,mum,sm2,
                   Ls,K,Xsz,celltotal,lambda){
  nsub = length(unique(celltotal))
  nb =length(B)
 
  A2_final <- Hess_modelparam_bulk(result,S,B,Xs,WB,Ub,ysc,
                                   ytrt,Lb,
                                   Ls,K,Xsz,celltotal,lambda)
  
  A_final <- A2_final/nb
  
  #covt = -solve(A_final)
  #B1_final <- Reduce("+",lapply(Lmat,function(x){x$B_cov}))
  
  gr_bulk=gradient_bulk_individual(result,S,B,Xs,WB,Ub,ysc,
                                   ytrt,Lb,
                                   Ls,K,Xsz,celltotal,lambda)
  bu_gr = colSums(gr_bulk)
  
  #B2_final <- sc_gr%*%t(bu_gr)
  #B3_final <- bu_gr%*%t(sc_gr)
  
  B2_final <- Reduce("+",lapply(1:nrow(gr_bulk),function(x){tcrossprod(gr_bulk[x,])}))
  B_final <- B2_final/nb
  
  #c(1,6:11),c(1,6:11)
  A_f0 <- A_final[c(1,6:11),c(1,6:11)]
  B_f0 <- B_final[c(1,6:11),c(1,6:11)]
  covt = ginv(A_f0)%*%B_f0%*%ginv(A_f0)/nb
  #covt = -t(bu_gr)%*%solve(A_f0)%*%bu_gr
  library(Matrix)
 
  
  pvalue <- sapply(4:6,function(z){
    C = matrix(rep(0,6),nrow=1)
    C[,z]=1
    DE_test(covt[2:7,2:7],result$alpha,C)
  })
  
  var3 = diag(covt)
  
  return(list(var=var3,pvalue =pvalue))
  
}
