# parameters: T, sample size; p1,p2, dimensions; k1,k2, factor numbers;
# tau: change point location, change: 0, 1 or 2, indicating null/alternative 1/ alternative 2
# pp: proportion (magnitude) of change, taking value from (0,1)
# a: control size of off-diagonal covariance entries
# cc: AR coefficients
gen.data<-function(Sample_T,p1,p2,k1,k2,tau=0.5,change=0,pp=0.3,a=0,cc=0){
  G1=matrix(a/p1,p1,p1)+diag(rep(1-a/p1,p1))#covariance matrices for errors
  G2=matrix(a/p2,p2,p2)+diag(rep(1-a/p2,p2))

  R=matrix(runif(p1*k1,min=-1,max=1),p1,k1)# loading matrices
  C=matrix(runif(p2*k2,min=-1,max=1),p2,k2)

  Y=array(0,c(Sample_T,p1,p2))

  if(change==0){## null hypothesis
    X=Y# common component, array
    E=Y# error, array
    F=matrix(rnorm(k1*k2),k1,k2)# factor matrix, rolling
    Er=rmatrixnorm(matrix(0,p1,p2),G1,G2)# error matrix, rolling
    for(t in 1:Sample_T){
      F=cc*F+sqrt(1-cc^2)*matrix(rnorm(k1*k2),k1,k2)# update factor
      X[t,,]=R%*%F%*%t(C)# update common component
      Er=cc*Er+sqrt(1-cc^2)*rmatrixnorm(matrix(0,p1,p2),G1,G2)# update error
      E[t,,]=Er
    }
    Y=X+E# data array
  } else if(change==1){## alternative 1
    # R.new=matrix(runif(p1*k1,min=-1,max=1),p1,k1)# loading matrix after change point
    rpn=ceiling(p1*pp)
    R.new=R
    R.new[sample(1:p1,rpn),]=runif(rpn,min=-1,max=1)
    X=Y# common component
    E=Y# error, array
    F=matrix(rnorm(k1*k2),k1,k2)#factor matrix, rolling
    Er=rmatrixnorm(matrix(0,p1,p2),G1,G2)# error matrix, rolling
    T1=ceiling(Sample_T*tau)# change point location
    for(t in 1:T1){# data sample before change point
      F=cc*F+sqrt(1-cc^2)*matrix(rnorm(k1*k2),k1,k2)# update F
      X[t,,]=R%*%F%*%t(C)# update X
      Er=cc*Er+sqrt(1-cc^2)*rmatrixnorm(matrix(0,p1,p2),G1,G2)# update E
      E[t,,]=Er
    }
    for(t in (T1+1):Sample_T){# data after change point
      F=cc*F+sqrt(1-cc^2)*matrix(rnorm(k1*k2),k1,k2)# update F
      X[t,,]=R.new%*%F%*%t(C)# update X using new loading
      Er=cc*Er+sqrt(1-cc^2)*rmatrixnorm(matrix(0,p1,p2),G1,G2)# update E
      E[t,,]=Er
    }
    Y=X+E
  } else {## alternative 2
    rpn=ceiling(p1*pp)
    R.add=matrix(0,p1,1)# loading correspondong to additional factor
    R.add[sample(1:p1,rpn),]=runif(rpn,min=-1,max=1)
    X=Y# common component
    E=Y# error
    F=matrix(rnorm(k1*k2),k1,k2)
    Er=rmatrixnorm(matrix(0,p1,p2),G1,G2)
    T1=ceiling(Sample_T*tau)# change point location
    for(t in 1:T1){#data before change
      F=cc*F+sqrt(1-cc^2)*matrix(rnorm(k1*k2),k1,k2)
      X[t,,]=R%*%F%*%t(C)
      Er=cc*Er+sqrt(1-cc^2)*rmatrixnorm(matrix(0,p1,p2),G1,G2)
      E[t,,]=Er
    }
    for(t in (T1+1):Sample_T){# data after change
      F=cc*F+sqrt(1-cc^2)*matrix(rnorm(k1*k2),k1,k2)
      X[t,,]=R%*%F%*%t(C)+R.add%*%t(rnorm(k2))%*%t(C)## add a new factor
      Er=cc*Er+sqrt(1-cc^2)*rmatrixnorm(matrix(0,p1,p2),G1,G2)
      E[t,,]=Er
    }
    Y=X+E
  }
  return(Y)
}

## generate psi_tau in a rolling scheme, without projection
# output: three series, corresponding to original eigenvalues, rescaled eigenvalues, and psi.tau
# parameters: Y, data array; r, which factor we want to monitor,
# m, training sample size; delta, shrinking; p, order of transformation function
gen.psi.tau.flat<-function(Y,r,m,delta,p){
  # this outputs a (Tm*3) matrix,
  # the first column corresponde to original sample eigenvalue lambda_{k+1} (rolling),
  # the second column corresponds to rescaled p1^{-delta}*lambda_{k+1}/trace,
  # the third for column is for psi_tau
  T=dim(Y)[1]
  p1=dim(Y)[2]
  p2=dim(Y)[3]
  Tm=T-m

  # a function to calculate sample covariance (unnormalized), x is p by n sample matrix
  f<-function(x){
    return(x%*%t(x))
  }

  #calculate psi.tau
  if (Tm>0){# ensure Tm is positive
    ## calculate M in rolling scheme
    # initial M
    temp=apply(Y[2:(m+1),,],1,f)
    # the apply function calculate faster than for-loop, it outputs a (p1^2)*m matrix, each
    # column corresponds to the vectorized sample covariance at time t (unnormalized)

    Mr.temp=matrix(apply(temp,1,sum),p1,p1)
    # this sums up all the columns of Y and output a p1*p1 matrix

    temp2=Mr.temp/(m*p2)# rescale, the initial M

    eigval=svds(temp2,r,r)$d # calculate the original leading r eigenvalues
    hat.lambda=eigval[r] # record the original eigenvalue

    lambda=p1^(1-delta)*eigval[r]/sum(diag(temp2))# rescaled eigenvalues
    res.lambda=lambda# record the rescaled eigenvalue

    psi.tau=lambda^p# psi_tau

    # rolling M
    for (i in 2:Tm){
      Mr.temp=Mr.temp-Y[i,,]%*%t(Y[i,,])+Y[m+i,,]%*%t(Y[m+i,,])# rolling M
      temp2=Mr.temp/(m*p2)# rescale
      eigval=svds(temp2,r,r)$d # update eigenvalues
      hat.lambda=c(hat.lambda,eigval[r]) # record

      lambda=p1^(1-delta)*eigval[r]/sum(diag(temp2))# update rescaled eigenvalue
      res.lambda=c(res.lambda,lambda)# record

      psi.tau=c(psi.tau,lambda^p)# update psi_tau
    }
    return(cbind(hat.lambda,res.lambda,psi.tau))
  } else {
    stop("error: too large m")
  }
}

## calculate psi_tau with projection
# parameters: Y, data array; r, which factor we want to monitor corresponding to (k+1);
# m, training sample size; delta, shrinking; p, transformation order
# kmax, upper-bound for the number of factors
gen.psi.tau.proj<-function(Y,r,m,delta,p,kmax){
  # this outputs a (Tm*3) matrix,
  # the first column corresponde to original sample eigenvalue lambda_{k+1} (rolling),
  # the second column corresponds to rescaled p1^{-delta}*lambda_{k+1}/trace,
  # the third for column is for psi_tau
  T=dim(Y)[1]
  p1=dim(Y)[2]
  p2=dim(Y)[3]
  Tm=T-m

  # a function to calculate sample covariance (unnormalized), x is p by n sample matrix
  f<-function(x){
    return(x%*%t(x))
  }

  #calculate psi.tau
  if (Tm>0){# ensure Tm is positive

    ## calculate M in rolling scheme
    # estimate C
    Y.tran=aperm(Y,c(1,3,2))# transpose Y[t,,]
    temp.tran=apply(Y.tran[2:(m+1),,],1,f)# calculate sample covariance for each t
    Mr.temp.tran=matrix(apply(temp.tran,1,sum),p2,p2)# sum up all t
    temp2.tran=Mr.temp.tran/(m*p1)# rescale
    C.hat=svds(temp2.tran,kmax,kmax)$u## estimate C, assuming kmax factor numbers

    # projected sample covariance function
    f.temp<-function(x){
      return(x%*%C.hat%*%t(C.hat)%*%t(x))
    }# similar to f, but need to update C each time

    # initial M
    temp=apply(Y[2:(m+1),,],1,f.temp)# sample covariance at each t
    Mr.temp=matrix(apply(temp,1,sum),p1,p1)# sum up all t
    temp2=Mr.temp/(m*p2)# rescale, the projected sample covariance matrix

    eigval=svds(temp2,r,r)$d # calculate the original leading r eigenvalues
    hat.lambda=eigval[r] # record the original eigenvale

    # lambda=p1^(-delta)*eigval[r]/((sum(diag(temp2))-sum(eigval))/min(m,p1))# rescaled eigenvalues
    lambda=p1^(-delta)*eigval[r]/((sum(diag(temp2)))/p1)# rescaled eigenvalues

    res.lambda=lambda# record the rescaled eigenvalue

    psi.tau=lambda^p

    # rolling M
    for (i in 2:Tm){
      # re-estimate C
      Mr.temp.tran=Mr.temp.tran-t(Y[i,,])%*%Y[i,,]+t(Y[m+i,,])%*%Y[m+i,,]# update
      temp2.tran=Mr.temp.tran/(m*p1)# rescale
      C.hat=svds(temp2.tran,kmax,kmax)$u## update C

      # projected sample covariance function
      f.temp<-function(x){
        return(x%*%C.hat%*%t(C.hat)%*%t(x))
      }# update
      temp=apply(Y[(i+1):(i+m),,],1,f.temp)# update
      Mr.temp=matrix(apply(temp,1,sum),p1,p1)# sum up
      temp2=Mr.temp/(m*p2)# rescale, update projected sample covarianc ematrix

      eigval=svds(temp2,r,r)$d # update eigenvalues
      hat.lambda=c(hat.lambda,eigval[r]) # record

      # lambda=p1^(-delta)*eigval[r]/((sum(diag(temp2))-sum(eigval))/min(m,p1))# rescaled eigenvalues
      lambda=p1^(-delta)*eigval[r]/((sum(diag(temp2)))/p1)# rescaled eigenvalues
      res.lambda=c(res.lambda,lambda)# record

      psi.tau=c(psi.tau,lambda^p)# update psi.tau

      # print(i)
    }
    return(cbind(hat.lambda,res.lambda,psi.tau))
  } else {
    stop("error: too large m")
  }
}

## test and locate change point, given psi.tau
# psi: the psi.tau series, length of Tm
# method: "ps"-partial sum, or "wc"-worst case
# eta: parameter for "ps"
# cv: a vector of critical values, corresonding to a series of significance level
test.matrix.psi<-function(m,psi,method="ps",eta=0.25,cv=2.3860){
  Tm=length(psi)
  y.tau=psi+rnorm(Tm)
  S.tau=cumsum(y.tau)#calculate S

  if(method=="ps"){
    re=cv
    loc=cv
    if(eta<0.5){
      Test=max(abs(S.tau)/((1:Tm)^(eta)))*Tm^(eta-0.5)
      for(i in 1:length(cv)){
        re[i]=(Test>cv[i])
        loc[i]=which((abs(S.tau)/((1:Tm)^(eta)))*Tm^(eta-0.5)>cv[i])[1]+m
      }
    }else if (eta==0.5){
      Test=max(abs(S.tau)/((1:Tm)^(0.5)))
      for(i in 1:length(cv)){
        re[i]=(Test>cv[i])
        loc[i]=which((abs(S.tau)/((1:Tm)^(0.5)))>cv[i])[1]+m
      }
    }else{
      rm=floor(log(Tm))
      Test=max(abs(S.tau[rm:Tm])/((rm:Tm)^(eta)))*rm^(eta-0.5)
      for(i in 1:length(cv)){
        re[i]=(Test>cv[i])
        loc[i]=which((abs(S.tau[rm:Tm])/((rm:Tm)^(eta)))*rm^(eta-0.5)>cv[i])[1]+m+rm
      }
    }
  }else{
    re=cv
    loc=cv
    for(i in 1:length(cv)){
      re[i]=(max(y.tau)>cv[i])
      loc[i]=which(y.tau>cv[i])[1]+m
    }
  }
  return(list(test=re,loc=loc))
}

## calculate critical values for eta=0.5 or "wc"
getcv<-function(Tm,alpha=0.05,method="ps"){
  if(method=="ps"){
    dm=2*log(log(Tm))+0.5*log(log(log(Tm)))-0.5*log(pi)
    cm=sqrt(2*log(log(Tm)))
    cv=(-log(-log(1-alpha))+dm)/cm
  }else{
    bm=sqrt(2*log(Tm))-log(4*pi*log(Tm))/(2*sqrt(2*log(Tm)))
    am=1/sqrt(2*log(Tm))
    cv=bm-log(-log(1-alpha))*am
  }
  return(cv)
}

## test and locate change point, given observation Y
# Y: T by pa by p2 array
#with projection
test.matrix.proj<-function(Y,r,m,epsilon,p,kmax,decrease=0,method="ps",eta=0.25,cv=2.3860){
  # Y: a T by p1 by p2 array
  # decrease: if True, we are testing whether factor number decreases
  T=dim(Y)[1]
  p1=dim(Y)[2]
  p2=dim(Y)[3]
  Tm=T-m
  if(Tm>0){
    temp=log(p1)/log(m*p2)
    delta=epsilon*(temp<=0.5)+(epsilon+1-1/(2*temp))*(temp>0.5)

    psi.tau=gen.psi.tau.proj(Y,r,m,delta,p,kmax)
    if(decrease==0){
      psi.tau=psi.tau[,3]
    }else{
      psi.tau=1/psi.tau[,3]
    }

    result=test.matrix.psi(m,psi.tau,method,eta,cv)
    return(result)
  } else {
    stop("error: too large m")
  }
}

#without projection
test.matrix.flat<-function(Y,r,m,epsilon,p,decrease=0,method="ps",eta=0.25,cv=2.3860){
  # Y: a T by p1 by p2 array
  # decrease: if True, we are testing whether factor number decreases
  T=dim(Y)[1]
  p1=dim(Y)[2]
  p2=dim(Y)[3]
  Tm=T-m
  if(Tm>0){
    temp=log(p1)/log(m*p2)
    delta=epsilon*(temp<=0.5)+(epsilon+1-1/(2*temp))*(temp>0.5)

    psi.tau=gen.psi.tau.flat(Y,r,m,delta,p)
    if(decrease==0){
      psi.tau=psi.tau[,3]
    }else{
      psi.tau=1/psi.tau[,3]
    }

    result=test.matrix.psi(m,psi.tau,method,eta,cv)
    return(result)

  } else {
    stop("error: too large m")
  }
}
