args = commandArgs(trailingOnly=TRUE)  # passed from script
vb <- as.numeric(args[1])
vc <- as.numeric(args[2])
theta.bulk <- as.numeric(args[3]) 
theta.sc <- as.numeric(args[4])
sg <- as.numeric(args[5])
typeoffset <- as.numeric(args[6])

setwd("D:/UNC/research/Xiaojing/Project1")
source("hessian_update.R")
source("weighted_vem.R")
source("variationalfunction_1013.R")
vb=1
vc=0.6
sg = 0
theta.sc = 0.4
theta.bulk = 1
typeoffset=2

# nb = 200
# ng = 1000
# K=3


ng=1000
K=3
nb=300
prop = c(30,30,30)
prob0 = prop/sum(prop)
W=LaplacesDemon::rdirichlet(nb,prop)
trt.bulk = rep(c(0,1),each=nb/2)
nsub = 80
ncell = 100
cell=ncell
ns = nsub*ncell
trt.sub = c(rep(0,nsub/2),rep(1,nsub/2))
trt = rep(trt.sub,each = ncell)

singlecmat = t(rmultinom(ns,1,prob = prob0))
singlec = apply(singlecmat,1,function(x){which(x==1)})
modmat = model.matrix(~-1+factor(singlec)+factor(singlec):trt)
#modmat = cbind(singlecmat,singlecmat * trt)

marker1 = 1:floor(ng/10)
marker2 = (floor(ng/10)+1):(floor(ng*2/10))
marker3 = ((floor(ng*2/10))+1):(floor(ng*3/10))

dd=3.5
d1=3
d2=4
d3=5
#-5, -7
alpha0 = matrix(rnorm(ng*K,dd,0.5),nrow = ng,ncol = K)
# alpha0[,1] = d1#marker1
# alpha0[,2] = d2#marker2
# alpha0[,3] = d3#marker3


alpha1 = matrix(0,nrow = ng,ncol = K)
deind1 = ((floor(ng*3/10))+1):(floor(ng*4/10))
deind2 = ((floor(ng*4/10))+1):(floor(ng*5/10))
deind3 = ((floor(ng*5/10))+1):(floor(ng*6/10))

ediff1 = vc
ediff2 = 0
ediff3 = vc
alpha1[deind1,1] = (1)^(deind1)*ediff1
alpha1[deind2,2] = (1)^(deind2)*ediff2
alpha1[deind3,3] = (1)^(deind3)*ediff3


alpha.true = cbind(alpha0,alpha1)

logmusc = alpha.true %*% t(modmat)

cellidx0 = rep(1:nsub,each = ncell)

libsizes.true = rnorm(ns,2,0.1)#rnorm(ns,5000,500)
libsizes.mat <- matrix(rep(libsizes.true,ng),nrow = ng, ncol = ns, byrow = TRUE)
sigmag = sg
b = matrix(0,nrow = ng,ncol = nsub, byrow=FALSE)
for(j in 1:ng){
  b[j,] = rnorm(nsub,0,sigmag)
}

lambda.sc.part2 = t(sapply(1:ng,function(j){Reduce("+",lapply(1:nsub,function(i){b[j,i]*(cellidx0 == i)}))}))
lambda.sc = exp(logmusc)*exp(lambda.sc.part2)*libsizes.mat

#theta.sc = 0.5
library(MASS)
Z0 = matrix(rnegbin(ng*ns,mu = lambda.sc, theta = theta.sc),nrow = ng,byrow=FALSE)

Y = Z0

Xsz = model.matrix(~trt)
gamma_true = c(1.5,0) ###let's model zero proportion

for(g in 1:ng)
{
  dropout = rbinom(ns,1,prob = exp(Xsz%*%gamma_true)/(1+exp(Xsz%*%gamma_true)))
  Y[g,] = Z0[g,]*(1-dropout)
}
Y = ceiling(Y)

libsizeb.true = rnorm(nb,100,5)
K_ir = apply(Y,1,function(x){exp(mean(log(x)))})
mu_bulk = matrix(0,nrow = ng,ncol = nb)
for(g in 1:ng){
  coefbulk = model.matrix(~trt.bulk) %*% rbind(alpha0[g,],alpha1[g,])
  for(i in 1:nb)
  {
    mu_bulk[g,i] = c(W[i,]%*%exp(coefbulk[i,]+log(libsizeb.true[i])))
  }
}

B0 = matrix(rnegbin(ng*nb,mu=mu_bulk,theta = theta.bulk),nrow = ng,ncol = nb,byrow = FALSE)



colnames(W) = paste0("c",1:K)
colnames(B0) = paste0("sample",1:nb)
colnames(Y) = paste0("cell",1:ns)

expmean = exp(apply(Y,1,function(x){mean(log(x[x>0]))}))
apply(Y,2,function(x){tmp = x/expmean
   median(tmp[tmp>0])})
#typeoffset=2
if(typeoffset == 1){
  
  library(DESeq2)
  ddsb = DESeqDataSetFromMatrix(countData = cbind(B0,Y+1), colData = data.frame(colnames(B0)), design = ~1)
  ddss = DESeqDataSetFromMatrix(countData = (Y+1), colData n n= data.frame(colnames(Y)), design = ~1)
  ddsb = estimateSizeFactors(ddsb)
  ddss = estimateSizeFactors(ddss)
  
  
  libsizeb = sizeFactors(ddsb)
  libsizes = sizeFactors(ddss)
}
 pseudo_expect = apply(Y,2,function(x){sapply(1:length(x),function(z){
   mean(x[-z])
 })})
 pseudo_expectB = apply(B0,2,function(x){sapply(1:length(x),function(z){
   mean(x[-z])
 })})
pt = exp(Xsz%*%gamma_true)/(1+exp(Xsz%*%gamma_true))
if(typeoffset==2){
  
  libsizeb = exp(colMeans(log(pseudo_expectB)))
  libsizes = exp(apply(pseudo_expect,2,function(x){mean(log(x))})-log(1-colMeans(Y==0)))
  libsizes = colMeans(Y)/(1-pt[1,1])
}

Xsz0 = model.matrix(~trt)

expmean = exp(apply(Y,1,function(x){mean(log(x[x>0]))}))
expmean_b = exp(apply(B0,1,function(x){mean(log(x[x>0]))}))


# 
libsizeb = libsizeb.true
libsizes = libsizes.true

B0 = floor(edgeR::cpm(B0))
Y = floor(edgeR::cpm(Y))
K=3
library(foreach)
library(doParallel)
ncores = detectCores()
cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl,cores = 10)
#mod = pscl::zeroinfl(Y[g,]~-1+factor(singlec)+factor(singlec):trt|trt,offset = libsizes,dist = "negbin")
prior_int =  sapply(1:K,function(x){
  tmp = log(1+rowMeans(Y[,singlec == x & trt == 0]))-log(rowMeans(1-(Y[,singlec==x & trt == 0]==0)))-mean(log(libsizes[singlec==x & trt ==0]))
 # tmp = rowMeans(log(1+Z0[,singlec == x])-log(libsizes[singlec==x]))
  
  tmp
  
}) ## 500 cells

prior_int1 =  sapply(1:K,function(x){
  tmp = rowMeans(log(1+pseudo_expect[,singlec == x]))-log(rowMeans(1-(Y[,singlec==x]==0)))-mean(log(libsizes[singlec==x]))
  # tmp = rowMeans(log(1+Z0[,singlec == x])-log(libsizes[singlec==x]))
  
  tmp
  
}) ## 500 cells

#prior_int= matrix(0,nrow=ng,ncol = K)
prior_trt = matrix(0,nrow=ng,ncol = K)

prior_mean= cbind(prior_int,prior_trt)


# S=Y[g,];B=B0[g,];a1=1;a2=1;param_initial = list(sigma2 = 0.7,beta = 3,theta = 1);
# prior_alpha = prior_mean[g,];Xs = modmat;WB = W;Ub=NULL;ysc= trt;
# ytrt = trt.bulk;cellidx = cellidx0;
# Xsz = model.matrix(~trt);Lb = libsizeb;eps = 1e-5;
# Ls =libsizes;K=K;max.iter=160

# system.time(result <- variational_weightedEM_matrix(S=Y[1:100,],B=B0[1:100,],a1=1,a2=1,param_initial = list(sigma2 = 0.7,beta = 3,theta = 1),
#                                 prior_alpha = prior_mean[1:100,],Xs = modmat,WB = W,Ub=NULL,ysc= trt,
#                                 ytrt = trt.bulk,cellidx = cellidx0,
#                                 Xsz = model.matrix(~trt),Lb = libsizeb,eps = 1e-5,
#                                 Ls =libsizes,K=K,max.iter=160))
# result04 = list()
# alpha_all = result$alpha
# beta_all = result$beta
# gamma_all = result$gamma
# theta_all = result$theta
# mubj_all = result$mubj
# sigmabj_all = result$sigmabj2
# pi_all = result$pi
# sigma2_all = result$sigma2
# for(i in 1:100){
#  
#   result04[[i]] = list(alpha = alpha_all[i,],beta = beta_all[i],
#                        gamma = gamma_all[i,],theta = theta_all[i],mubj = mubj_all[i,],
#                        sigmabj2 = sigmabj_all[i,],pi = pi_all[i,],sigma2 = sigma2_all[i])
# }

system.time(result07<- foreach (g =1:1000, .combine='c', .multicombine=TRUE) %dopar%{
  result = NULL
  result = try(variational_weightedEM(S=Y[g,],B=B0[g,],a1=1,a2=1,param_initial = list(sigma2 = 0.03,beta = 3,theta = 0.3),
                            prior_alpha = prior_mean[g,],Xs = modmat,WB = W,Ub=NULL,ysc= trt,
                            ytrt = trt.bulk,cellidx = cellidx0,
                            Xsz = model.matrix(~trt),Lb = libsizeb,eps = 1e-6,#mum = -0.5,sm2 = 0.5,
                            Ls =libsizes,K=K,max.iter=160,trace=1,lambda = 1e-4),silent = TRUE)
  list(result)})
result08 = c(result08,result09)
result10 = c(result10,result09)
apply(tmp2,2,var)
res09 <- NULL
resstd1 <- NULL
for(g in 1:1000){
  if(inherits(result07[[g]], "try-error")){next;}else{
    der =  VEM_cov(result = result07[[g]],
                        S=Y[g,],B=B0[g,],a1=1,a2=1,
                        Xs = modmat,WB = W,Ub=NULL,ysc= trt,
                        ytrt = trt.bulk,cellidx = cellidx0,mum=-0.5,sm2=0.5,
                        Xsz = model.matrix(~trt),Lb = libsizeb,Ls =libsizes,K=K,celltotal = cellidx0,
                        lambda = 1e-4)
    print(der$pvalue)
    resstd1 <- rbind(resstd1,der$var)
    res09 <- rbind(res09,der$pvalue)
  }
  
}

################ singlecell or bulk independent################
#################################################################
#################################################################
system.time(result04_sc <- foreach (g = 1:ng, .combine='c', .multicombine=TRUE) %dopar%{
  result = variational_weightedEM(S=Y[g,],B=B0[g,],a1 = 1,a2 = 0, param_initial = list(sigma2 = sg+0.1, beta = 3, theta = 1),
                                  prior_alpha = prior_mean[g,],Xs = modmat,WB = W,Ub=NULL,ysc= trt,
                                  ytrt = trt.bulk,cellidx = cellidx0,
                                  Xsz = model.matrix(~trt),Lb = libsizeb,eps = 1e-8,
                                  Ls =libsizes,K=K,max.iter=160,lambda = lambda0/nb)
  list(result)})




resf3 <- NULL
resstd0 <- NULL
for(g in 1:ng){
  der =  VEM_cov(result = result04_sc[[g]],
                 S=Y[g,],B=B0[g,],a1=1,a2=0,
                 Xs = modmat,WB = W,Ub=NULL,ysc= trt,
                 ytrt = trt.bulk,cellidx = cellidx0,
                 Xsz = model.matrix(~trt),Lb = libsizeb,Ls =libsizes,K=K,celltotal = cellidx0,lambda = lambda0/nb,mum=-0.5,sm2=0.5)
  print(der$pvalue)
  resstd0 <- rbind(resstd0,der$var)
  resf3 <- rbind(resf3,der$pvalue)
}

system.time(result04_bulk <- foreach (g = 1:ng, .combine='c', .multicombine=TRUE) %dopar%{
  result = variational_weightedEM(S=Y[g,],B=B0[g,],a1 = 0,a2 = 1, param_initial = list(sigma2 = sg+0.1, beta = 3, theta = 1),
                                  prior_alpha = prior_mean[g,],Xs = modmat,WB = W,Ub=NULL,ysc= trt,
                                  ytrt = trt.bulk,cellidx = cellidx0,
                                  Xsz = model.matrix(~trt),Lb = libsizeb,eps = 1e-8,
                                  Ls =libsizes,K=K,max.iter=160,lambda = lambda0/nb)
  list(result)})

resfbu <- NULL
resstd0 <- NULL
for(g in 1:ng){
  der =  VEM_cov(result = result04_bulk[[g]],
                 S=Y[g,],B=B0[g,],a1=0,a2=1,
                 Xs = modmat,WB = W,Ub=NULL,ysc= trt,
                 ytrt = trt.bulk,cellidx = cellidx0,
                 Xsz = model.matrix(~trt),Lb = libsizeb,Ls =libsizes,K=K,celltotal = cellidx0,lambda = lambda0/nb,mum=-0.5,sm2=0.5)
  print(der$pvalue)
  resstd0 <- rbind(resstd0,der$var)
  resfbu <- rbind(resfbu,der$pvalue)
}

#########################################
#############################################
#######################################################




system.time(result01 <- foreach (g =101:200, .combine='c', .multicombine=TRUE) %dopar%{
  result = variational_EM_sc(S=Y[g,],B=B0[g,],
                                  prior_alpha = prior_mean[g,],param_initial = list(beta = 5),Xs = modmat,WB = W,Ub=NULL,ysc= trt,
                                  ytrt = trt.bulk,cellidx = cellidx0,
                                  Xsz = model.matrix(~trt),Lb = libsizeb,eps = 1e-9,
                                  Ls =libsizes,K=K,max.iter=200)
  list(result)})

#########comparison with bMIND, TCA, ENIGMA###############
###########################################################

library(MIND)

prior_int =  sapply(1:K,function(x){
  tmp = Y[,singlec == x ]/libsizes[singlec==x]
  apply(log(tmp+1),1,function(x){mean(x[x>0])})
}) ## 500 cells

prior_case0 =  sapply(1:K,function(x){rowMeans(Y[,singlec == x & trt == 1])}) ## 500 cells
prior_control0 =  sapply(1:K,function(x){rowMeans(Y[,singlec == x & trt == 0])})
prior0 =  sapply(1:K,function(x){rowMeans(Y[,singlec == x])}) ## 500 cells

colnames(prior_control0) = colnames(prior_case0) = colnames(prior0) = paste0("c",1:K)
rownames(prior_control0) = rownames(prior_case0) = rownames(prior0) = 1:ng

colnames(B0) = paste0("sample",1:nb)
rownames(B0) = rownames(Y) = 1:ng

colnames(W) = paste0("c",1:K)
rownames(W) = colnames(B0)
colnames(Y) = 1:ns
library(TCA)
tca.mdl <- tca(X = as.matrix(log2(1+edgeR::cpm(B0))),
               W = W,
               C1 = NULL,
               C2 = NULL,
               parallel = TRUE,
               num_cores=10,               
               max_iter = 2)
tca.exp <- tensor(X = as.matrix(log2(1+edgeR::cpm(B0))), tca.mdl)
names(tca.exp) <- colnames(profile)
ptca = lapply(tca.exp,function(x){
  apply(x,1, function(y){summary(lm(y~trt.bulk))$coefficients[2,4]})
})

frac=  est_frac(sig = prior0[c(marker1,marker2,marker3),], bulk = B)
cov_calculate <- function(idx)
{
  cdx = cellidx[idx]
  sing = singlec[idx]
  Y0 = Y[,idx]
  cts = array(NA, dim = c(ng, length(idx), K))
  # rownames(cts) = rownames(sc)
  # colnames(cts) = sample
  # dimnames(cts)[[3]] = cell_type
  for (j in 1:nsub) {
    for (k in 1:K) {
      id = which(cdx == j & sing == 
                   k)
      if (length(id) > 0) 
        cts[, j, k] = rowMeans(Y0[, id, drop = F])
    }
  }
  cov = array(NA, dim = c(nrow(Y0), K, K))
  rownames(cov) = rownames(Y0)
  colnames(cov) = dimnames(cov)[[3]] = 1:K
  for (i in 1:nrow(cov)) {
    cov[i, , ] = cov(cts[i, , ], use = "pairwise")
  }
  return(cov)
}


covariance_co = cov_calculate(which(trt==0))
covariance_ca = cov_calculate(which(trt==1))
prior_case0 =  sapply(1:K,function(x){rowMeans(Y[,singlec == x & trt == 1])}) ## 500 cells
prior_control0 =  sapply(1:K,function(x){rowMeans(Y[,singlec == x & trt == 0])})
prior0 =  sapply(1:K,function(x){rowMeans(Y[,singlec == x]/(1-colMeans(Y[,singlec==x]==0)))}) ## 500 cells

prior0 = sapply(1:K,function(x){rowMeans(log(1+Y[,singlec == x]/(1-colMeans(Y[,singlec==x]==0))))})
colnames(prior_control0) = colnames(prior_case0) = colnames(prior0) = paste0("c",1:K)
rownames(prior_control0) = rownames(prior_case0) = rownames(prior0) = 1:ng

colnames(B0) = paste0("sample",1:nb)
rownames(B0) = rownames(Y) = 1:ng

colnames(W) = paste0("c",1:K)
rownames(W) = colnames(B0)
colnames(Y) = 1:ns
system.time(mod30 <- bmind_de(bulk = B0, frac = W, 
                              
                              max_samp = 1e+4, np = TRUE,y =trt.bulk, ncore = 10))

#mind0 = bMIND(bulk = B0, frac = W, profile = prior0, ncore = 10)
#pmind = MIND::test(mind0$A,trt.bulk)

system.time(mod30 <- bmind_de(bulk = log2(1+B0), frac = W, 
                            
                              max_samp = 1e+4, np = TRUE,y =trt.bulk, ncore = 10))


system.time(mod30 <- bmind_de(bulk = B0, frac = W, 
                              profile_co =  prior_control0, profile_ca = prior_case0,
                              covariance_co = covariance_co, covariance_ca = covariance_ca,
                              max_samp = 1e+4, np = F,y =trt.bulk, ncore = 10))

system.time(mod31 <- bmind_de(bulk = B, frac = frac, 
                              profile_co =  prior_control0, profile_ca = prior_case0,
                              covariance_co = covariance_co, covariance_ca = covariance_ca,
                              max_samp = 1e+4, np = F,y =trt.bulk, ncore = 10))

mind0 = bMIND(bulk = B0, frac = W, profile = prior0, ncore = 10)
pmind = MIND::test(mind0$A,trt.bulk)
library(MAST)
library(data.table)
mast.de <- function(dat,jm)
{
  sca <- FromMatrix(
    dat,
    cData = data.frame(t2d = trt[which(singlec == jm)],sample = cellidx0[which(singlec == jm)]),
    fData = data.frame(genename = rownames(dat)),
    class = "SingleCellAssay",check_sanity = FALSE
  )
  
  cond<-factor(colData(sca)$t2d)
  #cond<-relevel(cond,"Unstim")
  colData(sca)$condition<-cond
  zlmCond <- zlm(~condition+(1|sample), method = "glmer", ebayes = FALSE,sca)
  FCTHRESHOLD = log2(1.5)
  summaryCond <- summary(zlmCond, doLRT='condition1') 
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[contrast=='condition1' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast=='condition1' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  fcHurdle = fcHurdle[match(rownames(sim1.sc),primerid),]
  fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
  #setorder(fcHurdleSig, fdr)
  qval <- fcHurdle$fdr
  return(fcHurdle)
}


sim1.sc = Y[,singlec == 1]
sim2.sc = Y[,singlec == 2]
sim3.sc = Y[,singlec == 3]


fc1 <- mast.de(sim1.sc,1)
fc2 <- mast.de(sim2.sc,2)
fc3 <- mast.de(sim3.sc,3)


pval1 <- fc1$`Pr(>Chisq)`
pval2 <- fc2$`Pr(>Chisq)`
pval3 <- fc3$`Pr(>Chisq)`                              

library(ENIGMA)
egm = create_ENIGMA(bulk = B0, ref = prior0, ref_type = "aggre")

egm@result_cell_proportion = W
egm = ENIGMA_trace_norm(egm,solver = "proximalpoint", pos = TRUE,alpha = 0.1,preprocess= "none" )
DEG = FindCSE_DEG(egm,trt.bulk,FDR_control = FALSE)
glist = 1:ng

dc1 = DEG$c1[match(glist,rownames(DEG$c1)),]
dc2 = DEG$c2[match(glist,rownames(DEG$c2)),]
dc3 = DEG$c3[match(glist,rownames(DEG$c3)),]

####Type-1 error ########
m1 = mean(dc1$pvalue[-deind1]<.05)
m2 = mean(dc2$pvalue[-deind2]<.05)
m3 = mean(dc3$pvalue[common_null]<.05)

