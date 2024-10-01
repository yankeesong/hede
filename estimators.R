require(pacman)
options(rgl.useNULL=TRUE)
p_load(rgl,dplyr,glmnet,matrixStats,MASS,pracma,Rcpp,GenOrd,cccp,SHT,CovTools,PDSCE,RPtests)

###### our method
HEDE <- function(y,X,lambda.L.list=10^seq(2,-1,-0.01),
                       lambda.R.list=10^seq(3,-1,-0.02),
                       thresholds=c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1),
                       constraint=FALSE){
  
  # Remark: this function is designed to be optimized for high dimensions. 
  # i.e. the calculation for ridge DoF is the asymptotic formula.
  # Need to be modified if want more accurate solution in low dimension.
  
  # Author: Yanke Song (ysong@g.harvard.edu)
  # Runs ensembled debiased estimators for signal strength, noise and heritability
  #   estimation.
  # Model: y = X %*% beta + e, rows of X iid~ N(0,I),  e iid~ N(0,sigma^2)
  # Note: if constraint=TRUE, estimates are automatically clipped at 0.
  
  # Inputs:
  #   y: response
  #   X: covariates, each column is already normalized
  #   lambda.L.list: range of lambda for Lasso to search for
  #   lambda.R.list: range of lambda for Ridge to search for
  #   thresholds: a list of DoF-related thresholds for our method
  
  # Note: the input lambda is defined as in the loss function 1/2 RSS / n + lambda * penalty
  # where penalty is 1/2 ||_2^2 for ridge and ||_1 / sqrt(n) for lasso.
  # Need to rescale lambda to feed in glmnet
  p = ncol(X)
  n = nrow(X)
  kappa = p / n
  n_lambda_L = length(lambda.L.list)
  n_lambda_R = length(lambda.R.list)
  sd_y =  sqrt(sum(y^2)/n)
  
  
  em_s2_hat = rep(NA,n_lambda_L * n_lambda_R)
  helper_list = rep(NA,n_lambda_L * n_lambda_R)
  tau_hatsqC = rep(NA,n_lambda_L * n_lambda_R)
  
  ret = numeric(length(thresholds))
  
  
  # Compute Lasso
  fit_L=glmnet(X,y,family="gaussian",alpha=1,lambda=sd_y*lambda.L.list/sqrt(n),intercept=FALSE,standardize=FALSE,thresh=1e-10);
  beta_hatL_list=as.matrix(fit_L$beta)
  residualsL_list = matrix(replicate(n_lambda_L,y),nrow=n) - X%*%beta_hatL_list
  scoreL_list = t(X)%*%residualsL_list
  helperL_list = numeric(n_lambda_L)
  for (j in 1:n_lambda_L)
  {
    l0norm=nnzero(beta_hatL_list[,j])
    if (l0norm==n){l0norm = n-1}
    helperL_list[j] = 1 - l0norm/n
  }
  em_beta_debiasedL_list=beta_hatL_list+t(t(scoreL_list)/helperL_list/n)
  
  # Compute Ridge
  fit_R=glmnet(X,y,family="gaussian",alpha=0,lambda=sd_y*lambda.R.list,intercept=FALSE,standardize=FALSE,thresh=1e-10);
  beta_hatR_list=as.matrix(fit_R$beta)
  residualsR_list = matrix(replicate(n_lambda_R,y),nrow=n) - X%*%beta_hatR_list
  scoreR_list = t(X)%*%residualsR_list
  helperR_list = numeric(n_lambda_R)
  for (k in 1:n_lambda_R)
  {
    x = (-(lambda.R.list[k]+1-kappa) + sqrt((lambda.R.list[k]+1-kappa)^2+4*kappa*lambda.R.list[k]))/2/kappa
    helperR_list[k] = 1 - kappa * (1-x)
  }
  em_beta_debiasedR_list=beta_hatR_list+t(t(scoreR_list)/helperR_list/n)
  
  # Consider all linear combinations
  for (j in 1:n_lambda_L)
  {
    for (k in 1:n_lambda_R)
    {
      ind = (j-1)*n_lambda_R+k
      
      # Lasso
      residualsL= residualsL_list[,j]
      helper_L = helperL_list[j]
      tau_hatsqL=sum(residualsL^2)/(n^2*helper_L^2)
      
      # Ridge
      residualsR = residualsR_list[,k]
      helper_R = helperR_list[k]
      tau_hatsqR=sum(residualsR^2)/(n^2*helper_R^2)
      
      helper_list[ind] = min(helper_L,helper_R)
      
      # Then combination
      tauhatRL=sum(residualsL*residualsR)/(n^2*helper_L*helper_R)
      a_1 = (tau_hatsqR-tauhatRL)/(tau_hatsqL+tau_hatsqR-2*tauhatRL)
      a_1 = pmin(pmax(a_1, 0), 1)
      
      # Update value lists
      em_est_comb=a_1*em_beta_debiasedL_list[,j]+(1-a_1)*em_beta_debiasedR_list[,k]
      tau_hatsqC[ind] = a_1^2*tau_hatsqL+(1-a_1)^2*tau_hatsqR+2*a_1*(1-a_1)*tauhatRL
      em_s2_hat[ind] = sum(em_est_comb^2)-tau_hatsqC[ind]*p
      if (constraint==TRUE) {em_s2_hat[ind] = min(var(y),max(0,em_s2_hat[ind]))}
    }
  }
  
  
  # Filter by helper value
  for (i in 1:length(thresholds))
  {
    pval_top_index = which(helper_list > thresholds[i])
    if (length(pval_top_index) > 0)
    {
      indexC_p_h=pval_top_index[which.min(tau_hatsqC[pval_top_index])]
      s2 = em_s2_hat[indexC_p_h]
      if (var(y)==s2)
      {
        ret[i] = 1
      } else
      {
        ret[i] = s2 / var(y)
      }
    } else
    {
      ret[i] = NA
    }
  }
  return (ret)
}

#### Other fixed effect based methods
## Method of Moments in Dicker 2014.
Dicker <- function(y,X){
  n = nrow(X)
  p = ncol(X)
  norm1 = sum(y^2)
  norm2 = sum((t(X) %*% y)^2)
  tau2 = -p/n/(n+1) * norm1 + 1/n/(n+1) * norm2
  return (n*tau2/norm1)
}

## Method in Verzelen and Gassiat 18
VG18 <- function(y,X,lambda_scale=13){
  # We want to solve the scaled Lasso: ||y-X beta||_2 + 13* sqrt(log p) ||beta||_1
  # However the sqrt_lasso is solving: ||y-X beta||_2 / sqrt(n) + lambda ||beta||_1
  n = nrow(X)
  p = ncol(X)
  lambda = lambda_scale * sqrt(log(p) / n)
  beta_hat = sqrt_lasso(X,y,lam0=lambda,intercept=FALSE)
  return (1 - sum((y-X%*%beta_hat)^2) / sum(y^2))
}

## AMP in Bayati and Montanari 13
BM13 <- function(y,X,lambda=1){
  n = nrow(X)
  p = ncol(X)
  fit=glmnet(X,y,family="gaussian",alpha=1,lambda=lambda/sqrt(n),intercept=FALSE,standardize=FALSE);
  beta_hat=as.matrix(fit$beta)
  residual = y - X%*%beta_hat
  smallscore = t(X)%*%residual

  l0norm=nnzero(beta_hat)
  if (l0norm==n){l0norm = n-1}
  helper = 1 - l0norm/n
  tauhat2 = sum(residual^2)/(n^2*helper^2)
  sigmahat2 = tauhat2 * (n + p - 2*l0norm) - sum(smallscore^2) / (n-l0norm)^2
  return (1 - sigmahat2 / var(y))
}


## Estimating Equation method in Chen 22
EstEqn <- function(y,X){
  lambda = 0.1
  n = nrow(X)
  p = ncol(X)
  y = y / (sqrt((n-1)/n) * sd(y))
  M = X %*% t(X) / p
  for (t in 1:6)
  {
    helper = solve(diag(n)+lambda*M,M-diag(n))
    nume = helper %*% solve(diag(n)+lambda*M,y %*% t(y)-diag(n))
    deno = helper %*% helper
    r2hat = sum(diag(nume)) / sum(diag(deno))
    lambda = r2hat / (1-r2hat)
  }
  return (r2hat)
}

## CHIVE in Guo 2020
CHIVE <- function(y,X,lambda_scale=1){
  # We want to solve the scaled Lasso: ||y-X beta||_2 + scale * sqrt(2.01 * log p) ||beta||_1
  # However the sqrt_lasso is solving: ||y-X beta||_2 / sqrt(n) + lambda ||beta||_1
  n = nrow(X)
  p = ncol(X)
  lambda = lambda_scale * sqrt(2.01 * log(p) / n)
  beta_hat = sqrt_lasso(X,y,lam0=lambda,intercept=FALSE,standardize=FALSE)
  residual = y-X%*%beta_hat
  s2 = sum(beta_hat^2) + 2/n * t(beta_hat) %*% t(X) %*% residual
  return (s2 / var(y))
}

## EigenPrism in Janson 2015
EigenPrism <- function(y,X,invsqrtSig=NULL,alpha=0.05,target='beta2',zero.ind=c(),diagnostics=T,constraint=FALSE){
  # Author: Lucas Janson (statweb.stanford.edu/~ljanson)
  # Runs EigenPrism procedure for estimating and generating confidence
  #  intervals for variance components in high-dimensional linear model:
  #       y = X%*%beta + e,   rows of X iid~ N(0,Sig),   e iid~ N(0,sigma^2)
  #  Requires cccp package for solving second order cone optimization.
  #  Note confidence interval endpoints may lie outside parameter domain, so it may be appropriate
  #   to clip them after the fact.
  # 
  # Inputs:
  #  y: response vector of length n (will automatically be centered)
  #  X: n by p design matrix; columns will automatically be centered and scaled to variance 1;
  #      should not contain intercept column, since both y and X will be centered
  #  invsqrtSig: if columns of X not independent, p by p positive definite matrix which is the square-root
  #               of the inverse of Sig, where Sig is the *correlation* matrix of the X (default is identity)
  #  alpha: significance level for confidence interval (default = 0.05)
  #  target: target of estimation/inference
  #		  'beta2' (default) is the squared 2-norm of the coefficient vector: sum(beta^2)
  #           'sigma2' is the noise variance sigma^2
  #           'heritability' is the fraction of variance of y explained by X%*%beta: t(beta)%*%Sig%*%beta/var(y)
  #  zero.ind: vector of which indices of the weight vector w to constrain to zero (default is none)
  #  diagnostics: boolean (default = T) for whether to generate diagnostic plots for the V_i, lambda_i, and w_i
  #  
  # Outputs:
  #  estimate: unbiased estimate of the target (for heritability, only approximately unbiased)
  #  CI: 100*(1-alpha)% confidence interval for target
  
  # Get dimensionality of problem
  n = nrow(X)
  p = ncol(X)
  
  # Transform y and X to proper form
  #y = y-mean(y)
  X = scale(X,T,T)*n/(n-1)
  if(!is.null(invsqrtSig)) X = X%*%invsqrtSig
  
  # Take singular value decomposition and rescale singular values
  svd = svd(X)
  lambda = svd$d^2/p
  
  # Defined cone-constrained linear problem to optimize weights; [v; w] is vector of optimization variables
  q = c(1,rep(0,n)) #coefficient vector in objective function
  A = rbind(c(0,rep(1,n)),c(0,lambda)) #matrix for linear constraints
  b = c(0,1) #vector for linear constraints
  if(target=='sigma2') b = c(1,0) #switch constraints if target is sigma^2
  # Constrain some weights to be zero if desired
  if(!is.null(zero.ind)){
    A = rbind(A,cbind(rep(0,length(zero.ind)),diag(rep(1,n))[zero.ind,]))
    b = c(b,rep(0,length(zero.ind)))
  }
  # Define second-order cone constraints
  soc1 = socc(diag(c(1/4,rep(1,n))),c(-1/2,rep(0,n)),c(1/4,rep(0,n)),1/2)
  soc2 = socc(diag(c(1/4,lambda)),c(-1/2,rep(0,n)),c(1/4,rep(0,n)),1/2)
  prob = dlp(as.vector(q),A,as.vector(b),list(soc1,soc2))
  
  
  # Solve optimization problem and extract variables
  opt = cps(prob,ctrl(trace=F))
  v = getx(opt)[1]
  w = getx(opt)[-1]
  
  # Compute estimate and y's variance
  est = sum(w*(t(svd$u)%*%y)^2)
  #yvar = sum(y^2)/n
  yvar = var(y)
  if (constraint==TRUE)
  {
    est = max(0,min(yvar,est))
  }
  # Compute confidence interval
  CI = est + yvar*sqrt(v)*qnorm(1-alpha/2)*c(-1,1)
  if(target=='heritability'){
    est = est/yvar
    CI = CI/yvar
  }
  
  # Generate list with results
  result=list()
  result$estimate = est
  result$CI = CI
  
  # Generate diagnostic plots
  if(diagnostics){
    par(mfrow=c(1,3))
    
    # Check that eigenvectors are approximately Gaussian
    nV = floor(log10(n))
    srtV = svd$v[,10^(0:nV)]
    labs = c()
    for(i in 1:(nV+1)){
      srtV[,i] = sort(srtV[,i])
      ind = 10^(i-1)
      labs = c(labs,bquote(V[.(ind)]))
    }
    matplot(qnorm((1:p)/(p+1)),srtV,type="l",lwd=2,
            ylab="Quantiles of Eigenvectors",xlab="Gaussian Quantiles",
            main=expression(paste("Check Gaussianity of Eigenvectors ",V[i])))
    legend("topleft",as.expression(labs),col=1:(nV+1),lty=1:(nV+1),lwd=2)
    
    # Check that there are no outliers in the eigenvalues
    hist(lambda,main=expression(paste("Histogram of Normalized Eigenvalues ",lambda[i])),
         xlab=expression(lambda[i]))
    
    # Check that the weights are not dominated by just a few values
    srtw = sort(abs(w),T)
    plot(1:n,cumsum(srtw)/sum(srtw),type="l",lwd=2,
         main=expression(paste("Fraction of Total Weight in Largest k ",w[i])),
         xlab="k",ylab="Fraction of Total Weight")
  }
  
  return(result)
}