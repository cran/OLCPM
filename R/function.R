
#' generate data
#'
#' This function generates matrix-valued time series under a two-way factor
#' structure with/without a change point.
#'
#' See the paper He et al. (2021).
#'
#' @param Sample_T positive integer indicating the length of series.
#' @param p1 positive integer indicating the row dimension.
#' @param p2 positive integer indicating the column dimension.
#' @param k1 positive integer indicating the number of row factors.
#' @param k2 positive integer indicating the number of column factors.
#' @param tau a real number in \eqn{(0,1)}, indicating the location of change
#' point, i.e., (\eqn{\tau T}).
#' @param change the type of change, taking 0 for no change point, taking 1 for
#' the case that the loading matrix \bold{R} changes, taking other values for the case
#' that a new row factor occurs.
#' @param pp a number in \eqn{(0,1]}, indicating the magnitude of the break.
#' When \eqn{change=1}, \emph{pp} is the proportion of entries in \bold{R} that changes;
#' when \emph{change} is not equal to 0 or 1, \emph{pp} is the proportion of non-zero entries in the new factor
#' loading.
#' @param a a number in \eqn{[0,min(p_1,p_2))}, indicating the
#' cross-sectional correlations of the idiosyncratic errors.
#' @param cc a number in \eqn{[0,1)}, indicating the AR(1) coefficient of the
#' factor and error processes.
#' @return  a \eqn{T\times p1 \times p2} array.
#' @author Yong He, Xinbing Kong, Lorenzo Trapani, Long Yu
#' @references He Y, Kong X, Trapani L, & Yu L(2021). Online change-point
#' detection for matrix-valued time series with latent two-way factor
#' structure. \emph{arXiv preprint}, arXiv:2112.13479.
#'
#' @export
#'
#' @examples
#' # set parameters
#' k1=3
#' k2=3
#' epsilon=0.05
#' Sample_T=50
#' p1=40
#' p2=20
#'
#'
#' # generate data
#' Y=gen.data(Sample_T,p1,p2,k1,k2,tau=0.5,change=1,pp=0.3)
#' print("the dimension of Y is:")
#' print(dim(Y))
#'
#' @importFrom stats rnorm runif
#' @import RSpectra
#' @import LaplacesDemon

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


#' remove outliers
#'
#' This function removes outliers in the data, which are far from the sample
#' median.
#'
#' An outlier is detected if it deviates from the sample median more than
#' \emph{rg} times interquantile range.
#'
#' @param x a numeric data vector.
#' @param rg a positive number indicating how the outliers are defined; see
#' details.
#' @return  a list containing the following: \item{x}{a new
#' data vector with outliers replaced with \eqn{NA}. The original
#' missing values in the data are preserved.} \item{id}{the locations of the
#' outliers in the data vector.}
#' @author Yong He, Xinbing Kong, Lorenzo Trapani, Long Yu
#' @export
#' @examples
#'
#' a=c(1:5,NA,10000)
#' outlier.remove(a,3)


outlier.remove<-function(x,rg=3){
  upper=quantile(x,0.75,na.rm=T)
  lower=quantile(x,0.25,na.rm=T)
  me=median(x,na.rm = T)
  idx=abs(x-me)>rg*(upper-lower)
  x[idx]=NA
  return(list(x=x,id=which(idx)))
}


#' impute missing entries by linear interpolation
#'
#' This function imputes missing entries in a numeric vector by linear
#' interpolation.
#'
#' A missing entry will be imputed by linear interpolation with the two nearest
#' values before and after it in the vector. When all the values before (after)
#' it are missing, use the two nearest values after (before) it, instead.
#'
#' @param x a numeric data vector, where \eqn{NA} indicates missing
#' entries. The vector should contain at least two non-missing values.
#' @return  a new numeric vector, with the same size of the
#' original vector, while all the missing entries have been imputed.
#' @author Yong He, Xinbing Kong, Lorenzo Trapani, Long Yu
#' @examples
#'
#' a=c(NA,NA,3,NA,NA,6,NA,NA)
#' b=c(1,2,3,4.5,5,NA,6.5,7,NA)
#'
#' impute.linear(a)
#' impute.linear(b)
#'
#' @export

impute.linear<-function(x){
  missing_id=which(is.na(x)==1)
  n=length(x)
  nm=length(missing_id)
  if(nm<n-1){
    y=x
    for(i in missing_id){
      y1=x[1:(i-1)]
      y2=x[(i+1):n]
      id1=which(is.na(y1)==0)
      id2=which(is.na(y2)==0)
      if(length(id1)>0){
        if (length(id2)>0){
          d1=tail(id1,1)
          d2=head(id2,1)+i
          y[i]=((i-d1)*x[d2]+(d2-i)*x[d1])/(d2-d1)
        } else {
          d1=tail(id1,2)[1]
          d2=tail(id1,2)[2]
          y[i]=x[d2]+(i-d2)/(d2-d1)*(x[d2]-x[d1])
        }
      } else {
        d1=head(id2,2)[1]+i
        d2=head(id2,2)[2]+i
        y[i]=x[d1]-(d1-i)/(d2-d1)*(x[d2]-x[d1])
      }
    }
    return(y)
  } else {
    stop("error: too less valid entries")
  }
}

#' test whether the k-th moment exists
#'
#' This function tests the existence of k-th moment by randomized method in
#' Trapani (2016).
#'
#' The procedure is adapted from Trapani (2016) with \eqn{\psi=2}, where
#' \eqn{\psi} is a tuning parameter to scale the sample moments defined
#' in Section 3.1 of Trapani (2016). For simplicity, we only test the 4th, 6th,
#' ... 2c-th moments.
#'
#' @param x a numeric vector of data samples.
#' @param k a number no smaller than 4, indicating that the procedure will test
#' the existence of the k-th moment when k is even. Otherwise, the procedure
#' will test the existence of the \eqn{k'}-th moment, with
#' \eqn{k'=round(k/2,0)\times 2}.
#' @param R the number of standard Gaussian variables generated in the
#' randomized test; see more details in Trapani (2016).
#' @return a scalar in \eqn{[0,1]}, indicating the p-value
#' of the test. The null hypothese is that the k-th moment doesn't exist.
#' Therefore, a small p-value indicates the existense of the k-th moment.
#' @author Yong He, Xinbing Kong, Lorenzo Trapani, Long Yu
#' @references Trapani, L. (2016). Testing for (in) finite moments.
#' \emph{Journal of Econometrics}, 191(1), 57-68.
#' @examples
#'
#' x=rt(10000,5)
#' moment.test(x,4)
#'
#' x=rt(10000,4)
#' moment.test(x,4)
#'
#' @export
moment.test<-function(x,k=16,R=400){
  k=round(k/2,0)*2
  if(k>0){
    if(sum(is.na(x))>0){
      x=x[-(which(is.na(x)))]
    }
    mu1=mean(abs(x)^k)
    mu2=mu1/((mean(abs(x)^2))^(k/2))
    temp=2*(1:(k/2))-1
    mu3=mu2/(prod(temp))
    mu4=sqrt(exp(mu3))

    xi=rnorm(R)
    z1=((mu4*xi)<1)
    z2=((mu4*xi)<(-1))
    t1=sum(z1-0.5)*2/sqrt(R)
    t2=sum(z2-0.5)*2/sqrt(R)

    test=mean(t1^2,t2^2)

    return(1-pchisq(test,1))
  }else{
    stop("error: too small k")
  }
}

#' determine the moment (largest) of the data samples
#'
#' This function reports the largest moment that exists for a collection of
#' data samples.
#'
#' The procedure will sequentially test the existence of the \eqn{4th, 6th,
#' 8th, ... k.max-th} moment, using the function \eqn{moment.test} in the same
#' package. As soon as the procedure finds that the \eqn{k-th} moment does not
#' exist, it stops and reports at most \eqn{(k-1)-th} moment.
#'
#' @param x a numeric vector of data samples.
#' @param k.max a number indicating the upper bound, i.e., at most k.max-th
#' moment exists.
#' @param alpha a number in \eqn{(0,1)}, indicating the significance level of
#' the test.
#' @param R the number of standard Gaussian variables generated in the
#' randomized test; see also \code{\link{moment.test}}.
#' @return  an integer, indicating the largest moment that
#' exists for the data samples.
#' @author Yong He, Xinbing Kong, Lorenzo Trapani, Long Yu
#' @examples
#'
#' x=rt(10000,5)
#' moment.determine(x,10)
#'
#' x=rt(10000,4)
#' moment.determine(x,10)
#'
#'
#' @export
#'
moment.determine<-function(x,k.max=8,alpha=0.05,R=400){
  k.opt=c(0,(2:round(k.max/2,0))*2)
  if(max(k.opt)>0){
    for(i in 2:length(k.opt)){
      temp.p=moment.test(x,k.opt[i],R)
      if(temp.p>alpha){
        return(k.opt[i-1])
        break
      }
    }
    return(k.opt[i])
  }else{
    stop("error: too small k.max")
  }
}


#' calculate eigenvalue series by ``flat'' method
#'
#' This function calculates the rolling eigenvalue series for the monitoring
#' process, based on the ``flat'' version of sample covanriance matrix.
#'
#' The rolling eigenvalue series will start at the stage \eqn{m+1}, with length
#' \eqn{T-m}.
#'
#' @param Y the observed \eqn{T\times p1\times p2} array. \eqn{T} is the
#' sample size, \eqn{p1} and \eqn{p2} are the row and column dimensions,
#' respectively.
#' @param k a positive integer determining which eigenvalue to monitor.
#' \eqn{k=1} for the largest eigenvalue.
#' @param m a positive integer (\eqn{>1}) indicating the bandwidth of the
#' rolling windom.
#' @param delta a number in \eqn{(0,1)} indicating the rescaling parameter for
#' the eigenvalue. The default approach to calcualte delta is in the paper He
#' et al. (2021).
#' @param r a positive integer indicating the order of the transformation
#' function \eqn{g(x)=|x|^r}. Motivated by the paper, \eqn{r} should be chosen
#' according to the moments of the data; see more details in He et al. (2021).
#' @return  a \eqn{(T-m)\times 3} matrix, whose three columns are the original,
#'  rescaled, and transformed eigenvalue series, respectively.
#' @author Yong He, Xinbing Kong, Lorenzo Trapani, Long Yu
#' @references He Y, Kong X, Trapani L, & Yu L(2021). Online change-point
#' detection for matrix-valued time series with latent two-way factor
#' structure. \emph{arXiv preprint}, arXiv:2112.13479.
#' @examples
#'
#' ## generate data
#' k1=3
#' k2=3
#' epsilon=0.05
#' Sample_T=50
#' p1=40
#' p2=20
#' kmax=8
#' r=8
#' m=p2
#'
#' # generate data
#' Y=gen.data(Sample_T,p1,p2,k1,k2,tau=0.5,change=1,pp=0.3)
#'
#' # calculate delta
#' temp=log(p1)/log(m*p2)
#' delta=epsilon*(temp<=0.5)+(epsilon+1-1/(2*temp))*(temp>0.5)
#'
#' # calculate psi.tau
#' psi2=gen.psi.tau.flat(Y,k1+1,m,delta,r)
#' print(psi2)
#'
#' @export
#'
gen.psi.tau.flat<-function(Y,k,m,delta,r){
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

    eigval=svds(temp2,k,k)$d # calculate the original leading r eigenvalues
    hat.lambda=eigval[k] # record the original eigenvalue

    lambda=p1^(1-delta)*eigval[k]/sum(diag(temp2))# rescaled eigenvalues
    res.lambda=lambda# record the rescaled eigenvalue

    psi.tau=lambda^r# psi_tau

    # rolling M
    for (i in 2:Tm){
      Mr.temp=Mr.temp-Y[i,,]%*%t(Y[i,,])+Y[m+i,,]%*%t(Y[m+i,,])# rolling M
      temp2=Mr.temp/(m*p2)# rescale
      eigval=svds(temp2,k,k)$d # update eigenvalues
      hat.lambda=c(hat.lambda,eigval[k]) # record

      lambda=p1^(1-delta)*eigval[k]/sum(diag(temp2))# update rescaled eigenvalue
      res.lambda=c(res.lambda,lambda)# record

      psi.tau=c(psi.tau,lambda^r)# update psi_tau
    }
    return(cbind(hat.lambda,res.lambda,psi.tau))
  } else {
    stop("error: too large m")
  }
}

#' calculate eigenvalue series by projected method
#'
#' This function calculates the rolling eigenvalue series for the monitoring
#' process, based on the projected version of sample covanriance matrix.
#'
#' The rolling eigenvalue series will start at the stage \eqn{m+1}, with length
#' \eqn{T-m}.
#'
#' @param Y the observed \eqn{T\times p1\times p2} array. \eqn{T} is the
#' sample size, \eqn{p1} and \eqn{p2} are the row and column dimensions,
#' respectively.
#' @param k a positive integer determining which eigenvalue to monitor.
#' \eqn{k=1} for the largest eigenvalue.
#' @param m a positive integer (\eqn{>1}) indicating the bandwidth of the
#' rolling windom.
#' @param delta a number in \eqn{(0,1)} indicating the rescaling parameter for
#' the eigenvalue. The default approach to calcualte delta is in the paper He
#' et al. (2021).
#' @param r a positive integer indicating the order of the transformation
#' function \eqn{g(x)=|x|^r}. Motivated by the paper, \eqn{r} should be chosen
#' according to the moments of the data; see more details in He et al. (2021).
#' @param kmax a positive integer indicating the column number of the
#' projection matrix, should be larger than 0 but smaller than \eqn{p2}.
#' @return  a \eqn{(T-m)\times 3} matrix, whose three columns are the original,
#'  rescaled, and transformed eigenvalue series, respectively.
#' @author Yong He, Xinbing Kong, Lorenzo Trapani, Long Yu
#' @references He Y, Kong X, Trapani L, & Yu L(2021). Online change-point
#' detection for matrix-valued time series with latent two-way factor
#' structure. \emph{arXiv preprint}, arXiv:2112.13479.
#' @examples
#'
#' ## generate data
#' k1=3
#' k2=3
#' epsilon=0.05
#' Sample_T=50
#' p1=40
#' p2=20
#' kmax=8
#' r=8
#' m=p2
#'
#' # generate data
#' Y=gen.data(Sample_T,p1,p2,k1,k2,tau=0.5,change=1,pp=0.3)
#'
#' # calculate delta
#' temp=log(p1)/log(m*p2)
#' delta=epsilon*(temp<=0.5)+(epsilon+1-1/(2*temp))*(temp>0.5)
#'
#' # calculate psi.tau
#' psi2=gen.psi.tau.proj(Y,k1+1,m,delta,r,kmax)
#' print(psi2)
#'
#' @export
#'
gen.psi.tau.proj<-function(Y,k,m,delta,r,kmax){
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

    eigval=svds(temp2,k,k)$d # calculate the original leading r eigenvalues
    hat.lambda=eigval[k] # record the original eigenvale

    # lambda=p1^(-delta)*eigval[r]/((sum(diag(temp2))-sum(eigval))/min(m,p1))# rescaled eigenvalues
    lambda=p1^(-delta)*eigval[k]/((sum(diag(temp2)))/p1)# rescaled eigenvalues

    res.lambda=lambda# record the rescaled eigenvalue

    psi.tau=lambda^r

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

      eigval=svds(temp2,k,k)$d # update eigenvalues
      hat.lambda=c(hat.lambda,eigval[k]) # record

      lambda=p1^(-delta)*eigval[k]/((sum(diag(temp2)))/p1)# rescaled eigenvalues
      res.lambda=c(res.lambda,lambda)# record

      psi.tau=c(psi.tau,lambda^r)# update psi.tau

      # print(i)
    }
    return(cbind(hat.lambda,res.lambda,psi.tau))
  } else {
    stop("error: too large m")
  }
}

#' calculate critical values
#'
#' This function calculates critical values for the partial-sum and worst-case
#' statistics.
#'
#' For the partial-sum statistic with \eqn{\eta=0.5} or the worst-case
#' statistic, the critical value is simply \eqn{-log(-log(1-alpha))}. For the
#' partial-sum statistic with \eqn{\eta} not equal to 0.5, the critical
#' value of the scaled Wiener process is approximated by simulated data or from
#' our preserved table \emph{cv.table}, covering \eqn{\eta} in
#' \eqn{[0.01,0.49]} with step size equal to 0.01 and \eqn{\alpha} in
#' \eqn{[0.001,0.500]} with step size equal to 0.001. See more details for the
#' test statistics in He et al. (2021).
#'
#' @param alpha a number in \eqn{(0,1)}, indicating the significance level of
#' the test.
#' @param method ``ps'' for the partial-sum staistic, others for the worst-case
#' statistic.
#' @param eta a number in \eqn{[0,1]}, a scaling parameter required for "ps"
#' method; see more details in He et al. (2021).
#' @param simul logical value, woking only for "ps" method with
#' \eqn{\eta} not equal to 0.5.  When \emph{simul} is true, the
#' function will return approximated critical values based on 50000
#' replications of simulated Wiener process on a grid of 10000 points in
#' \eqn{[0,1]}. Otherwise, the function first checks for the nearest pair of
#' \eqn{(\eta,\alpha)} in the preserved \emph{cv.table}, and
#' then returns the corresponding critical value.
#' @return  a real number.
#' @author Yong He, Xinbing Kong, Lorenzo Trapani, Long Yu
#' @references He Y, Kong X, Trapani L, & Yu L(2021). Online change-point
#' detection for matrix-valued time series with latent two-way factor
#' structure. \emph{arXiv preprint}, arXiv:2112.13479.
#' @examples
#' getcv(0.05,method="ps",eta=0.25)
#' getcv(0.05,method="ps",eta=0.25,simul=1)
#' getcv(0.10,method="wc")
#' @export
#'
getcv<-function(alpha=0.05,method="ps",eta=0.5,simul=0){
  cv=-log(-log(1-alpha))
  if(method=="ps"){
    if(eta<0.5){
      if(simul==0){
        eta.id=which.min(abs(eta-seq(0,0.49,0.01)))
        alpha.id=c()
        for(alpha.temp in alpha){
          alpha.id=c(alpha.id,which.min(abs(alpha-seq(0.001,0.500,0.001))))
        }
        cv=OLCPM::cv.table[eta.id,alpha.id]
      }else{
        simul.data=c()
        for(j in 1:20){
          temp.cv=apply(matrix(rnorm(2500*10000),2500,10000)/sqrt(10000),1,cumsum)
          simul.data=c(simul.data,apply(abs(temp.cv)/(((1:10000)/10000)^eta),2,max))
          remove(temp.cv)
        }
        cv=quantile(simul.data,1-alpha)
      }
    }else if(eta>0.5){
      eta1=1-eta
      if(simul==0){
        eta.id=which.min(abs(eta1-seq(0,0.49,0.01)))
        alpha.id=c()
        for(alpha.temp in alpha){
          alpha.id=c(alpha.id,which.min(abs(alpha-seq(0.001,0.500,0.001))))
        }
        cv=OLCPM::cv.table[eta.id,alpha.id]
      }else{
        simul.data=c()
        for(j in 1:20){
          temp.cv=apply(matrix(rnorm(2500*10000),2500,10000)/sqrt(10000),1,cumsum)
          simul.data=c(simul.data,apply(abs(temp.cv)/(((1:10000)/10000)^eta1),2,max))
          remove(temp.cv)
        }
        cv=quantile(simul.data,1-alpha)
      }
    }
  }
  return(cv)
}

#' test single change point for matrix-valued online data given rolling
#' eigenvalue series
#'
#' This function tests single change point for matrix-valued online time
#' series, under a two-way factor structure, given the transformed eigenvalue
#' series.
#'
#' See He et al. (2021).
#'
#' @param m a positive integer (\eqn{>1}) indicating the bandwidth of the
#' rolling windom.
#' @param psi the transformed eigenvalue series, produced by gen.psi.tau.flat
#' or gen.psi.tau.proj, with length \eqn{T-m}.
#' @param method indicating the test statistic, ``ps'' for the partial-sum
#' method, while others for the worst-case method.
#' @param eta a number between \eqn{[0,1)}, indicating the parameter \eqn{\eta}
#' used in the partial-sum statistic.
#' @param cv critical value, related to the significance level and test
#' statistic. The default cv is from Horvath et al. (2004), and only works for
#' \eqn{\eta=0.25} or \eqn{\eta=0.75}. For other cases, generate the critical
#' value first by function \code{\link{getcv}}. Note that for the partial-sum
#' statistic with \eqn{\eta} not equal to 0.5, the critical values are
#' approximated by simulated data, thus can be slightly different from those in
#' Horvath et al. (2004).
#' @return  a list containing:
#' \item{test}{a logical value. 1 indicating the existence of change point, 0
#' indicating no change point.}
#' \item{loc}{an integer larger than m, indicating
#' the location of change point; or \eqn{NA} when no change point is reported.}
#' @author Yong He, Xinbing Kong, Lorenzo Trapani, Long Yu
#' @references Horvath L, Huskova M, Kokoszka P, et al (2004). Monitoring
#' changes in linear models. \emph{Journal of statistical Planning and
#' Inference}, 126(1): 225-251.
#'
#' He Y, Kong X, Trapani L, & Yu L(2021). Online change-point detection for
#' matrix-valued time series with latent two-way factor structure. \emph{arXiv
#' preprint}, arXiv:2112.13479.
#' @examples
#'
#' k1=3
#' k2=3
#' epsilon=0.05
#' Sample_T=50
#' p1=40
#' p2=20
#' kmax=8
#' r=8
#' m=p2
#'
#' # generate data
#' Y=gen.data(Sample_T,p1,p2,k1,k2,tau=0.5,change=1,pp=0.5)
#'
#' # calculate delta
#' temp=log(p1)/log(m*p2)
#' delta=epsilon*(temp<=0.5)+(epsilon+1-1/(2*temp))*(temp>0.5)
#'
#' # calculate psi.tau
#' psi1=gen.psi.tau.proj(Y,k1+1,m,delta,r,kmax)
#' psi2=gen.psi.tau.flat(Y,k1+1,m,delta,r)
#'
#' # calculate cv for "ps" with eta=0.45 and "wc"
#' cv1=getcv(0.05,method="ps",eta=0.45)
#' cv2=getcv(0.05,method="wc")
#'
#'
#' # test with psi1
#' test.once.psi(m,psi1[,3],method="ps",eta=0.45,cv1)
#'
#' test.once.psi(m,psi1[,3],method="wc",eta=0.5,cv2)
#'
#'
#' # test with psi2
#' test.once.psi(m,psi2[,3],method="ps",eta=0.45,cv1)
#'
#' test.once.psi(m,psi2[,3],method="wc",eta=0.5,cv2)
#'
#'
#' @export
#'
test.once.psi<-function(m=20,psi,method="ps",eta=0.25,cv=2.3860){
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
        if(re[i]==1){
          loc[i]=which((abs(S.tau)/((1:Tm)^(eta)))*Tm^(eta-0.5)>cv[i])[1]+m
        }else{
          loc[i]=NA
        }
      }
    }else if (eta==0.5){
      Test=max(abs(S.tau)/((1:Tm)^(0.5)))
      dm=2*log(log(Tm))+0.5*log(log(log(Tm)))-0.5*log(pi)
      cm=sqrt(2*log(log(Tm)))
      Test=Test*cm-dm
      for(i in 1:length(cv)){
        re[i]=(Test>cv[i])
        if(re[i]==1){
          loc[i]=which((abs(S.tau)/((1:Tm)^(0.5))*cm-dm)>cv[i])[1]+m
        }else{
          loc[i]=NA
        }
      }
    }else{
      rm=floor(log(Tm))
      if(rm>Tm){
        re=rep(FALSE,length(cv))
        loc=rep(NA,length(cv))
      }else{
        Test=max(abs(S.tau[rm:Tm])/((rm:Tm)^(eta)))*rm^(eta-0.5)
        for(i in 1:length(cv)){
          re[i]=(Test>cv[i])
          if(re[i]==1){
            loc[i]=which((abs(S.tau[rm:Tm])/((rm:Tm)^(eta)))*rm^(eta-0.5)>cv[i])[1]+m+rm
          }else{
            loc[i]=NA
          }
        }
      }
    }
  }else{
    re=cv
    loc=cv
    bm=sqrt(2*log(Tm))-log(4*pi*log(Tm))/(2*sqrt(2*log(Tm)))
    am=bm/(1+bm^2)
    for(i in 1:length(cv)){
      re[i]=((max(y.tau)-bm)/am>cv[i])
      if(re[i]==1){
        loc[i]=which((y.tau-bm)/am>cv[i])[1]+m
      }else{
        loc[i]=NA
      }
    }
  }
  return(list(test=re,loc=loc))
}

#' robust test of single change point for matrix-valued online data given
#' rolling eigenvalue series
#'
#' Based on \code{\link{test.once.psi}}, this function
#' repeats the randomized procedure multiple times and reports the majority
#' vote, thus more robust.
#'
#' See He et al. (2021).
#'
#' @param m a positive integer (\eqn{>1}) indicating the bandwidth of the
#' rolling windom.
#' @param psi the transformed eigenvalue series, produced by gen.psi.tau.flat
#' or gen.psi.tau.proj, with length \eqn{T-m}.
#' @param method indicating the test statistic, ``ps'' for the partial-sum
#' method, while others for the worst-case method.
#' @param eta a number between \eqn{[0,1)}, indicating the parameter \eqn{\eta}
#' used in the partial-sum statistic.
#' @param cv critical value; see also \code{\link{test.once.psi}}.
#' @param S an integer indicating the number of replications.
#' @param pr an number in \eqn{(0,1])}. The procedure reports a change point
#' only when the proportion of positive votes is over \emph{pr} in the
#' \emph{S} replications.
#' @return a list containing:
#' \item{test}{a logical value. 1 indicating the existence of change point, 0
#' indicating no change point.} \item{loc}{an integer larger than m, indicating
#' the median location of the change point among the positive votes in the
#' \emph{S} replications; or \eqn{NA} when no change point is reported.}
#' @author Yong He, Xinbing Kong, Lorenzo Trapani, Long Yu
#' @references  He Y, Kong X, Trapani L, & Yu L(2021). Online change-point detection for
#' matrix-valued time series with latent two-way factor structure. \emph{arXiv
#' preprint}, arXiv:2112.13479.
#' @examples
#'
#' k1=3
#' k2=3
#' epsilon=0.05
#' Sample_T=50
#' p1=40
#' p2=20
#' kmax=8
#' r=8
#' m=p2
#'
#' # generate data
#' Y=gen.data(Sample_T,p1,p2,k1,k2,tau=0.5,change=1,pp=0.5)
#'
#' # calculate delta
#' temp=log(p1)/log(m*p2)
#' delta=epsilon*(temp<=0.5)+(epsilon+1-1/(2*temp))*(temp>0.5)
#'
#' # calculate psi.tau
#' psi1=gen.psi.tau.proj(Y,k1+1,m,delta,r,kmax)
#' psi2=gen.psi.tau.flat(Y,k1+1,m,delta,r)
#'
#' # calculate cv for "ps" with eta=0.45 and "wc"
#' cv1=getcv(0.05,method="ps",eta=0.45)
#' cv2=getcv(0.05,method="wc")
#'
#'
#' # test with psi1
#' test.once.psi.robust(m,psi1[,3],method="ps",eta=0.45,cv1,S=100,pr=0.75)
#'
#' test.once.psi.robust(m,psi1[,3],method="wc",eta=0.5,cv2,S=100,pr=0.75)
#'
#'
#' # test with psi2
#' test.once.psi.robust(m,psi2[,3],method="ps",eta=0.45,cv1,S=100,pr=0.75)
#'
#' test.once.psi.robust(m,psi2[,3],method="wc",eta=0.5,cv2,S=100,pr=0.75)
#'
#'
#' @export
#'
test.once.psi.robust<-function(m=20,psi,method="ps",eta=0.25,cv=2.3860,S=100,pr=0.75){
  result=matrix(0,S,2)
  for(i in 1:S){
    temp=test.once.psi(m,psi,method,eta,cv)
    result[i,]=c(temp$test,temp$loc)
  }
  if(mean(result[,1])>pr){
    return(list(test=1,loc=median(result[,2],na.rm = T)))
  }else{
    return(list(test=0,loc=NA))
  }
}

#' test single change point for matrix-valued online time series-projected
#' version
#'
#' This function tests single change point for matrix-valued online time
#' series, under a two-way factor structure, using projected sample covariance
#' matrix.
#'
#' See He et al. (2021).
#'
#' @param Y data, a \eqn{T\times p1\times p2} array.
#' @param k a positive integer indicating which eigenvalue to monitor.
#' \eqn{k=1} for the largest eigenvalue.
#' @param m a positive integer (\eqn{>1}) indicating the bandwidth of the
#' rolling window.
#' @param epsilon the rescaling parameter taking value in \eqn{(0,1)}; see He
#' et al. (2021).
#' @param r a positive integer indicating the order of the transformation
#' function \eqn{g(x)=|x|^r}; see also \code{\link{gen.psi.tau.proj}}.
#' @param kmax a positive number determining the column number of the
#' projection matrix, should be larger than 0 but smaller than \eqn{p2}.
#' @param decrease a logical value. If \emph{decrease=1}, testing the
#' decrease of factor number.
#' @param method indicating the test statistic, ``ps'' for the partial-sum
#' method; others for the worst-case method.
#' @param eta a number between \eqn{[0,1)}, indicating the parameter \eqn{\eta}
#' used in the partial-sum statistic.
#' @param cv critical value; see also \code{\link{test.once.psi}}.
#' @return a list containing:
#' \item{test}{a logical value. 1 indicating the existence of change point, 0
#' indicating no change point.} \item{loc}{an integer larger than m, indicating
#' the location of change point; or \eqn{NA} when no change point is reported.}
#' @author Yong He, Xinbing Kong, Lorenzo Trapani, Long Yu
#' @references  He Y, Kong X, Trapani L, & Yu L(2021). Online change-point detection for
#' matrix-valued time series with latent two-way factor structure. \emph{arXiv
#' preprint}, arXiv:2112.13479.
#' @examples
#'
#' k1=3
#' k2=3
#' epsilon=0.05
#' Sample_T=50
#' p1=40
#' p2=20
#' kmax=8
#' r=8
#' m=p2
#'
#' # generate data
#' Y=gen.data(Sample_T,p1,p2,k1,k2,tau=0.5,change=1,pp=0.5)
#'
#' # calculate cv for "ps" with eta=0.45 and "wc"
#' cv1=getcv(0.05,method="ps",eta=0.45)
#' cv2=getcv(0.05,method="wc")
#'
#'
#' ## test with Y, projection
#' test.once.proj(Y,k1+1,m,epsilon,r,kmax,0,method="ps",eta=0.25)
#'
#'
#' test.once.proj(Y,k1+1,m,epsilon,r,kmax,0,method="ps",eta=0.45,cv1)
#'
#'
#' test.once.proj(Y,k1+1,m,epsilon,r,kmax,0,method="wc",eta=0.5,cv2)
#'
#' @export
#'
test.once.proj<-function(Y,k=1,m=20,epsilon=0.05,r=8,kmax=3,decrease=0,method="ps",eta=0.25,cv=2.3860){
  T=dim(Y)[1]
  p1=dim(Y)[2]
  p2=dim(Y)[3]
  Tm=T-m
  if(Tm>0){
    temp=log(p1)/log(m*p2)
    delta=epsilon*(temp<=0.5)+(epsilon+1-1/(2*temp))*(temp>0.5)

    psi.tau=gen.psi.tau.proj(Y,k,m,delta,r,kmax)
    if(decrease==0){
      psi.tau=psi.tau[,3]
    }else{
      psi.tau=1/psi.tau[,3]
    }

    result=test.once.psi(m,psi.tau,method,eta,cv)
    return(result)
  } else {
    stop("error: too large m")
  }
}

#' robust test of single change point for matrix-valued online time
#' series-projected version
#'
#' Based on \code{\link{test.once.proj}}, this function
#' repeats the randomized procedure multiple times and reports the majority
#' vote, thus more robust.
#'
#' See He et al. (2021).
#'
#' @param Y data, a \eqn{T\times p1\times p2} array.
#' @param k a positive integer indicating which eigenvalue to monitor.
#' \eqn{k=1} for the largest eigenvalue.
#' @param m a positive integer (\eqn{>1}) indicating the bandwidth of the
#' rolling windom.
#' @param epsilon the rescaling parameter taking value in \eqn{(0,1)}; see He
#' et al. (2021).
#' @param r a positive integer indicating the order of the transformation
#' function \eqn{g(x)=|x|^r}; see also  \code{\link{gen.psi.tau.proj}}.
#' @param kmax a positive number determining the column number of the
#' projection matrix, should be larger than 0 but smaller than \eqn{p_2}.
#' @param decrease a logical value. If \emph{decrease=1}, testing the
#' decrease of factor number.
#' @param method indicating the test statistic, ``ps'' for the partial-sum
#' method; others for the worst-case method.
#' @param eta a number between \eqn{[0,1)}, indicating the parameter \eqn{\eta}
#' used in the partial-sum statistic.
#' @param cv critical value; see also \code{\link{test.once.psi}}.
#' @param S an integer indicating the number of replications.
#' @param pr an number in \eqn{(0,1])}. The procedure reports a change point
#' only when the proportion of positive votes is over \emph{pr} in the
#' \emph{S} replications.
#' @return a list containing:
#' \item{test}{a logical value. 1 indicating the existence of change point, 0
#' indicating no change point.} \item{loc}{an integer larger than m, indicating
#' the median location of the change point among the positive votes in the
#' \emph{S} replications; or \eqn{NA} when no change point is reported.}
#' @author Yong He, Xinbing Kong, Lorenzo Trapani, Long Yu
#' @references He Y, Kong X, Trapani L, & Yu L(2021). Online change-point detection for
#' matrix-valued time series with latent two-way factor structure. \emph{arXiv
#' preprint}, arXiv:2112.13479.
#' @examples
#'
#' k1=3
#' k2=3
#' epsilon=0.05
#' Sample_T=50
#' p1=40
#' p2=20
#' kmax=8
#' r=8
#' m=p2
#'
#' # generate data
#' Y=gen.data(Sample_T,p1,p2,k1,k2,tau=0.5,change=1,pp=0.5)
#'
#' # calculate cv for "ps" with eta=0.45 and "wc"
#' cv1=getcv(0.05,method="ps",eta=0.45)
#' cv2=getcv(0.05,method="wc")
#'
#'
#' ## test with Y, projection
#' test.once.proj.robust(Y,k1+1,m,epsilon,r,kmax,0,method="ps",eta=0.25)
#'
#'
#' test.once.proj.robust(Y,k1+1,m,epsilon,r,kmax,0,method="ps",eta=0.45,cv1)
#'
#'
#' test.once.proj.robust(Y,k1+1,m,epsilon,r,kmax,0,method="wc",eta=0.5,cv2)
#'
#' @export
#'
test.once.proj.robust<-function(Y,k=1,m=20,epsilon=0.05,r=8,kmax=3,decrease=0,method="ps",eta=0.25,cv=2.3860,S=100,pr=0.75){
  result=matrix(0,S,2)
  for(i in 1:S){
    temp=test.once.proj(Y,k,m,epsilon,r,kmax,decrease,method,eta,cv)
    result[i,]=c(temp$test,temp$loc)
  }
  if(mean(result[,1])>pr){
    return(list(test=1,loc=median(result[,2],na.rm = T)))
  }else{
    return(list(test=0,loc=NA))
  }
}

#' test single change point for matrix-valued online time series -''flat''
#' version
#'
#' This function tests single change point for matrix-valued online time
#' series, under a two-way factor structure, using ''flat'' sample covariance
#' matrix.
#'
#' See He et al. (2021).
#'
#' @param Y data, a \eqn{T\times p1\times p2} array.
#' @param k a positive integer indicating which eigenvalue to monitor.
#' \eqn{k=1} for the largest eigenvalue.
#' @param m a positive integer (\eqn{>1}) indicating the bandwidth of the
#' rolling windom.
#' @param epsilon the rescaling parameter taking value in \eqn{(0,1)}; see He
#' et al. (2021).
#' @param r a positive integer indicating the order of the transformation
#' function \eqn{g(x)=|x|^r}; see also  \code{\link{gen.psi.tau.proj}}.
#' @param decrease a logical value. If \emph{decrease=1}, testing the
#' decrease of factor number.
#' @param method indicating the test statistic, ``ps'' for the partial-sum
#' method; others for the worst-case method.
#' @param eta a number between \eqn{[0,1)}, indicating the parameter \eqn{\eta}
#' used in the partial-sum statistic.
#' @param cv critical value; see also  \code{\link{test.once.psi}}.
#' @return a list containing:
#' \item{test}{a logical value. 1 indicating the existence of change point, 0
#' indicating no change point.} \item{loc}{an integer larger than m, indicating
#' the location of change point; or \eqn{NA} when no change point is reported.}
#' @author Yong He, Xinbing Kong, Lorenzo Trapani, Long Yu
#' @references He Y, Kong X, Trapani L, & Yu L(2021). Online change-point detection for
#' matrix-valued time series with latent two-way factor structure. \emph{arXiv
#' preprint}, arXiv:2112.13479.
#' @examples
#'
#' k1=3
#' k2=3
#' epsilon=0.05
#' Sample_T=50
#' p1=40
#' p2=20
#' r=8
#' m=p2
#'
#' # generate data
#' Y=gen.data(Sample_T,p1,p2,k1,k2,tau=0.5,change=1,pp=0.5)
#'
#' # calculate cv for "ps" with eta=0.45 and "wc"
#' cv1=getcv(0.05,method="ps",eta=0.45)
#' cv2=getcv(0.05,method="wc")
#'
#'
#' ## test with Y, flat version
#' test.once.flat(Y,k1+1,m,epsilon,r,0,method="ps",eta=0.25)
#'
#'
#' test.once.flat(Y,k1+1,m,epsilon,r,0,method="ps",eta=0.45,cv1)
#'
#'
#' test.once.flat(Y,k1+1,m,epsilon,r,0,method="wc",eta=0.5,cv2)
#'
#' @export
#'
test.once.flat<-function(Y,k=1,m=20,epsilon=0.05,r=8,decrease=0,method="ps",eta=0.25,cv=2.3860){
  # Y: a T by p1 by p2 array
  # decrease: if True, we are testing whether factor number decreases
  T=dim(Y)[1]
  p1=dim(Y)[2]
  p2=dim(Y)[3]
  Tm=T-m
  if(Tm>0){
    temp=log(p1)/log(m*p2)
    delta=epsilon*(temp<=0.5)+(epsilon+1-1/(2*temp))*(temp>0.5)

    psi.tau=gen.psi.tau.flat(Y,k,m,delta,r)
    if(decrease==0){
      psi.tau=psi.tau[,3]
    }else{
      psi.tau=1/psi.tau[,3]
    }

    result=test.once.psi(m,psi.tau,method,eta,cv)
    return(result)

  } else {
    stop("error: too large m")
  }
}

#' robust test of single change point for matrix-valued online time series
#' -"flat" version
#'
#' Based on \code{\link{test.once.flat}}, this function
#' repeats the randomized procedure multiple times and reports the majority
#' vote, thus more robust.
#'
#' See He et al. (2021).
#'
#' @param Y data, a \eqn{T\times p1\times p2} array.
#' @param k a positive integer indicating which eigenvalue to monitor.
#' \eqn{k=1} for the largest eigenvalue.
#' @param m a positive integer (\eqn{>1}) indicating the bandwidth of the
#' rolling windom.
#' @param epsilon the rescaling parameter taking value in \eqn{(0,1)}; see He
#' et al. (2021).
#' @param r a positive integer indicating the order of the transformation
#' function \eqn{g(x)=|x|^r}; see also  \code{\link{gen.psi.tau.proj}}.
#' @param decrease a logical value. If \emph{decrease=1}, testing the
#' decrease of factor number.
#' @param method indicating the test statistic, ``ps'' for the partial-sum
#' method; others for the worst-case method.
#' @param eta a number between \eqn{[0,1)}, indicating the parameter \eqn{\eta}
#' used in the partial-sum statistic.
#' @param cv critical value; see also  \code{\link{test.once.psi}}.
#' @param S an integer indicating the number of replications.
#' @param pr an number in \eqn{(0,1])}. The procedure reports a change point
#' only when the proportion of positive votes is over \emph{pr} in the
#' \emph{S} replications.
#' @return a list containing:
#' \item{test}{a logical value. 1 indicating the existence of change point, 0
#' indicating no change point.} \item{loc}{an integer larger than m, indicating
#' the median location of the change point among the positive votes in the
#' \emph{S} replications; or \eqn{NA} when no change point is reported.}
#' @author Yong He, Xinbing Kong, Lorenzo Trapani, Long Yu
#' @references He Y, Kong X, Trapani L, & Yu L(2021). Online change-point detection for
#' matrix-valued time series with latent two-way factor structure. \emph{arXiv
#' preprint}, arXiv:2112.13479.
#' @examples
#'
#' k1=3
#' k2=3
#' epsilon=0.05
#' Sample_T=50
#' p1=40
#' p2=20
#' kmax=8
#' r=8
#' m=p2
#'
#' # generate data
#' Y=gen.data(Sample_T,p1,p2,k1,k2,tau=0.5,change=1,pp=0.5)
#'
#' # calculate cv for "ps" with eta=0.45 and "wc"
#' cv1=getcv(0.05,method="ps",eta=0.45)
#' cv2=getcv(0.05,method="wc")
#'
#'
#' ## test with Y, flat version
#' test.once.flat.robust(Y,k1+1,m,epsilon,r,0,method="ps",eta=0.25)
#'
#'
#' test.once.flat.robust(Y,k1+1,m,epsilon,r,0,method="ps",eta=0.45,cv1)
#'
#'
#' test.once.flat.robust(Y,k1+1,m,epsilon,r,0,method="wc",eta=0.5,cv2)
#'
#' @export
#'
test.once.flat.robust<-function(Y,k=1,m=20,epsilon=0.05,r=8,decrease=0,method="ps",eta=0.25,cv=2.3860,S=100,pr=0.75){
  result=matrix(0,S,2)
  for(i in 1:S){
    temp=test.once.flat(Y,k,m,epsilon,r,decrease,method,eta,cv)
    result[i,]=c(temp$test,temp$loc)
  }
  if(mean(result[,1])>pr){
    return(list(test=1,loc=median(result[,2],na.rm = T)))
  }else{
    return(list(test=0,loc=NA))
  }
}

#' robust test of multiple change point for matrix-valued online time series
#'
#' This function tests multiple change points for matrix-valued online time
#' series, under a two-way factor structure. A change point will be reported
#' only when it's the majority vote in multiple replications. The function \code{\link{KSTP}}
#' is used to determine the initial number of factors in each regime. This function only outputs
#' the change points for row factors. For column factors, transpose the data.
#'
#'#' See empirical study in He et al. (2021).
#'
#' @param Y data, a \eqn{T\times p1\times p2} array.
#' @param k a non-negative integer indicating the initial number of factors.
#' @param m a positive integer (\eqn{>1}) indicating the bandwidth of the rolling
#' windom.
#'
#' @param epsilon1 the rescaling parameter taking value in \eqn{(0,1)}, for the
#' test of new factors or the change of loading space; see He et al. (2021).
#' @param epsilon2 the rescaling parameter taking value in \eqn{(0,1)}, for the
#' test of vanishing factors; see He et al. (2021).
#' @param r a positive integer indicating the order of the transformation
#' function \eqn{g(x)=|x|^r}; see also  \code{\link{gen.psi.tau.proj}}.
#' @param type indicates how to calculate the sample covariance. "flat" for the
#' flat version, while others for the projected version. See more details in He
#' et al. (2021).
#' @param kmax a positive number determining the column number of the projection
#' matrix, should be larger than 0 but smaller than \eqn{p_2}, required when
#' \emph{type} not being "flat".
#' @param method indicating the test statistic, ``ps'' for the partial-sum
#' method; others for the worst-case method.
#' @param eta a number between \eqn{[0,1)}, indicating the parameter \eqn{\eta}
#' used in the partial-sum statistic.
#' @param cv critical value; see also \code{\link{test.once.psi}}.
#' @param S an integer indicating the number of replications.
#' @param pr an number in \eqn{(0,1])}. The procedure reports a change point only
#' when the proportion of positive votes is over \emph{pr} in the
#' \emph{S} replications.
#' @return a matrix with two columns. The first column reports the locations
#' of change points. The second column reports the number of row factors
#' after each change point.
#' @author Yong He, Xinbing Kong, Lorenzo Trapani, Long Yu
#' @references He Y, Kong X, Trapani L, & Yu L(2021). Online change-point detection for
#' matrix-valued time series with latent two-way factor structure. \emph{arXiv
#' preprint}, arXiv:2112.13479.
#'
#' @examples
#' k1=3
#' k2=3
#' Sample_T=100
#' p1=40
#' p2=40
#' kmax=8
#' r=8
#' m=p2
#'
#' # generate data
#' Y1=gen.data(Sample_T,p1,p2,k1,k2,tau=0.5,change=1,pp=0.5)
#' Y2=gen.data(Sample_T,p1,p2,k1,k2,tau=0.5,change=0)
#' Y=array(rbind(matrix(Y1,Sample_T,p1*p2),matrix(Y2,Sample_T,p1*p2)),c(Sample_T*2,p1,p2))
#'
#' # calculate cv for "ps" with eta=0.45 and "wc"
#' cv1=getcv(0.05,method="ps",eta=0.45)
#' cv2=getcv(0.05,method="wc")
#'
#' # test with Y
#' test.multiple.robust(Y,k1,m,epsilon1=0.25,epsilon2=0.05,r,type="proj",kmax,method="ps")
#'
#' test.multiple.robust(Y,k1,m,epsilon1=0.25,epsilon2=0.05,r,type="proj",kmax,method="wc",cv=cv2)
#'
#' test.multiple.robust(Y,k1,m,epsilon1=0.25,epsilon2=0.05,r,type="flat",method="wc",cv=cv2)
#'
#' test.multiple.robust(Y,k1,m,epsilon1=0.25,epsilon2=0.05,r,type="flat",method="ps",eta=0.45,cv=cv1)
#' @export
#'
test.multiple.robust<-function(Y,k=1,m=20,epsilon1=0.25,epsilon2=0.05,r=8,kmax=4,type="proj",method="ps",eta=0.25,cv=2.3860,S=100,pr=0.75){
  result=vector()
  T=dim(Y)[1]
  p1=dim(Y)[2]
  p2=dim(Y)[3]

  while(T>m+3){
    if(k<min(p1,p2)){
      if(type=="flat"){
        re1=test.once.flat.robust(Y,k+1,m,epsilon1,r,decrease=0,method,eta,cv,S,pr)
      }else{
        re1=test.once.proj.robust(Y,k+1,m,epsilon1,r,kmax,decrease=0,method,eta,cv,S,pr)
      }
    }else{
      re1=list(test=0,loc=NA)
    }
    if(k>0){
      if(type=="flat"){
        re2=test.once.flat.robust(Y,k,m,epsilon2,r,decrease=1,method,eta,cv,S,pr)
      }else{
        re2=test.once.proj.robust(Y,k,m,epsilon2,r,kmax,decrease=1,method,eta,cv,S,pr)
      }
    }else{
      re2=list(test=0,loc=NA)
    }

    if(re1$test==1 &re2$test==1){
      loc=min(re1$loc,re2$loc)
      label=(re1$loc<re2$loc)-(re1$loc>re2$loc)
    }else if(re1$test==1 &re2$test==0){
      loc=re1$loc
      label=1
    }else if(re1$test==0 &re2$test==1){
      loc=re2$loc
      label=-1
    }else{
      break
    }
    if(label==0){
      break
    }
    if(T-loc>m){
      Y=Y[-(1:loc),,]
      T=dim(Y)[1]
      k=KSTP(Y[1:m,,],type=type,kmax=kmax,epsilon=(epsilon1+epsilon2)/2)
    }else{
      TT=dim(Y)[1]
      k=KSTP(Y[(TT-m+1):TT,,],type=type,kmax=kmax,epsilon=(epsilon1+epsilon2)/2)
      Y=Y[-(1:loc),,]
      T=dim(Y)[1]
    }

    if(is.null(nrow(result))!=1){
      result=rbind(result,c(loc+result[nrow(result),1],k))
    }else{
      result=rbind(result,c(loc,k))
    }
  }
  return(result)
}




#'determine factor number - projected
#'
#' This function determined the numbers of row and column factors
#' for matrix variate data under a two-way factor model, using a projected method.
#'
# See Yu et al. (2022).
#'
#' @references Yu L, He Y, Kong X, & Zhang X (2022). Projected estimation for
#'  large-dimensional matrix factor models. \emph{Journal of Econometrics},
#'  229(1),201-217.
#' @param Y data, a \eqn{T\times p1\times p2} array.
#' @param kmax a positive integer smaller than p2, indicating the
#' upper bound for the factor numbers, and the dimension of projection matrix.
#' @param c a non-negative but small number, to ensure the
#' denominator of the eigenvalue ratio statistics is not 0.
#'
#' @return a 2-dimensional vector containing the estimated
#' number of row and column factors, respectively.
#' @export
#'
#' @examples
#' k1=3
#' k2=3
#' Sample_T=100
#' p1=40
#' p2=20
#' kmax=8
#'
#' Y=gen.data(Sample_T,p1,p2,k1,k2,tau=0.5,change=0)
#' kpe(Y,kmax)
#'
#'

kpe<-function(Y,kmax=4,c=0){
  T=dim(Y)[1]
  p1=dim(Y)[2]
  p2=dim(Y)[3]
  d1=(1/sqrt(T*p1)+1/sqrt(T*p2)+1/p2)*c
  d2=(1/sqrt(T*p1)+1/sqrt(T*p2)+1/p1)*c

  M1=matrix(0,p1,p1)
  for(t in 1:T){
    M1=M1+Y[t,,]%*%t(Y[t,,])
  }
  R1=svds(M1,kmax,kmax,kmax)$u
  M2=matrix(0,p2,p2)
  for(t in 1:T){
    M2=M2+t(Y[t,,])%*%Y[t,,]
  }
  C1=svds(M2,kmax,kmax,kmax)$u

  k2=0
  k1=0
  iter=0
  k1_new=kmax
  k2_new=kmax
  while ((k1_new != k1 |k2_new!=k2) & iter<10) {
    k1=k1_new
    k2=k2_new
    M2=matrix(0,p2,p2)
    for(t in 1:T){
      M2=M2+t(Y[t,,])%*%R1[,1:k1]%*%t(R1[,1:k1])%*%Y[t,,]
    }
    eigval=eigen(M2)$values[1:(1+kmax)]
    k2_new=which.max(eigval[1:kmax]/(eigval[2:(1+kmax)]+d2))

    M1=matrix(0,p1,p1)
    for(t in 1:T){
      M1=M1+Y[t,,]%*%C1[,1:k2_new]%*%t(C1[,1:k2_new])%*%t(Y[t,,])
    }
    eigval=eigen(M1)$values[1:(1+kmax)]
    k1_new=which.max(eigval[1:kmax]/(eigval[2:(1+kmax)]+d1))
    iter=iter+1
  }
  return(c(k1,k2))
}


#' testing the number of row factors- without projection
#'
#' This function tests whether the number of row factors is equal or larger
#' than a given integer, under a two-way factor model, using flat version
#' of sample covariance.
#'
#' See He et al. (2023)
#' @references He Y, Kong X, Trapani L, & Yu L (2023).
#' One-way or two-way factor model for matrix sequences? \emph{Journal of Econometrics},
#' 235(2), 1981-2004.
#'
#' @param Y data, a \eqn{T\times p1\times p2} array.
#' @param k an positive integer indicating which eigenvalue to test.
#' @param alpha a number in (0,1), indicating the significance of the test.
#' @param epsilon a small positive number in (0,1), indicating the size of scaling.
#' @param r a positive number indicating the order of the power function
#' for transforming the rescaled eigenvalue.
#' @param M a large integer for the number of Gaussian variables in the randomized test.
#' @param S another large integer for the number of replications in the strong rule. Usually \eqn{M=S=T}.
#' @param fq a number in (0,0.5), controlling the threshold function of the strong rule.
#'
#' @return a logical value. 1 for "the number of row factors is smaller than k".
#' 0 for "at least k row factors exists".
#' @export
#'
#' @examples
#' k1=3
#' k2=3
#' Sample_T=100
#' p1=40
#' p2=20
#'
#' Y=gen.data(Sample_T,p1,p2,k1,k2,tau=0.5,change=0)
#' ITP_noproj(Y,k=1,M=Sample_T,S=Sample_T)
#' ITP_noproj(Y,k=4,M=Sample_T,S=Sample_T)
#'
ITP_noproj=function(Y,k=1,alpha=0.05,epsilon=0.05,r=8,M=100,S=100,fq=1/4){
  c1=qchisq((1-alpha),1)## critical value for chisq^2
  T=dim(Y)[1]
  p1=dim(Y)[2]
  p2=dim(Y)[3]
  M1=matrix(0,p1,p1)
  ## Calculate M1 (scaled)
  for(i in 1:T){
    M1=M1+Y[i,,]%*%t(Y[i,,])
  }
  M1=M1/(T*p2)


  beta=log(p1)/log(T*p2)
  delta=epsilon*(beta<0.5|beta==0.5)+(1-1/(2*beta)+epsilon)*(beta>0.5)

  ## calculate phi
  lam1=eigen(M1)$values ##  eigenvalues of M1
  lamhat=lam1/(sum(lam1)/p1)
  phihat=(abs(p1^(-delta)*lamhat))^r


  ## testing the factor number H_0: k1>=r with no projection technique
  Psihat=NULL
  for (s in 1:S){
    eta=rnorm(M/min(k,2))
    psihat1=(sqrt(phihat[k])*eta)<0.7
    psihat2=(sqrt(phihat[k])*eta)<(-0.7)
    psihat3=(sqrt(phihat[k])*eta)<2.4
    psihat4=(sqrt(phihat[k])*eta)<(-2.4)
    nuhat1=2*sum(psihat1-0.5)/sqrt(M)
    nuhat2=2*sum(psihat2-0.5)/sqrt(M)
    nuhat3=2*sum(psihat3-0.5)/sqrt(M)
    nuhat4=2*sum(psihat4-0.5)/sqrt(M)
    Psihat[s]=0.45*nuhat1^2+0.45*nuhat2^2+0.05*nuhat3^2+0.05*nuhat4^2## testing tatistics
  }
  Qhat=1/S*(sum(Psihat<=c1))

  reject_hat=1-1*(Qhat>=(1-alpha-S^(-fq)))

  return(reject_hat)
}

#' testing the number of row factors- with projection
#'
#' This function tests whether the number of row factors is equal or larger
#' than a given integer, under a two-way factor model, using projected version
#' of sample covariance.
#'
#' See He et al. (2023)
#' @references He Y, Kong X, Trapani L, & Yu L (2023).
#' One-way or two-way factor model for matrix sequences? \emph{Journal of Econometrics},
#' 235(2), 1981-2004.
#'
#' @param Y data, a \eqn{T\times p1\times p2} array.
#' @param k an positive integer indicating which eigenvalue to test.
#' @param alpha a number in (0,1), indicating the significance of the test.
#' @param kmax a positive integer smaller than p2, indicating the
#' upper bound for the factor numbers, and the dimension of projection matrix.
#' @param epsilon a small positive number in (0,1), indicating the size of scaling.
#' @param r a positive number indicating the order of the power function
#' for transforming the rescaled eigenvalue.
#' @param M a large integer for the number of Gaussian variables in the randomized test.
#' @param S another large integer for the number of replications in the strong rule. Usually \eqn{M=S=T}.
#' @param fq a number in (0,0.5), controlling the threshold function of the strong rule.
#'
#' @return a logical value. 1 for "the number of row factors is smaller than k".
#' 0 for "at least k row factors exists".
#' @export
#'
#' @examples
#' k1=3
#' k2=3
#' Sample_T=100
#' p1=40
#' p2=20
#'
#' Y=gen.data(Sample_T,p1,p2,k1,k2,tau=0.5,change=0)
#' ITP_proj(Y,k=1,M=Sample_T,S=Sample_T)
#' ITP_proj(Y,k=4,M=Sample_T,S=Sample_T)
#'
ITP_proj=function(Y,k=1,alpha=0.05,kmax=4,epsilon=0.05,r=8,M=100,S=100,fq=1/4){
  c1=qchisq((1-alpha),1)## critical value for chisq^2
  T=dim(Y)[1]
  p1=dim(Y)[2]
  p2=dim(Y)[3]

  ##  eigenvalues of projected method, k2 is set as kmax to estimate C
  M2=matrix(0,p2,p2)
  for(t in 1:T){
    M2=M2+t(Y[t,,])%*%Y[t,,]
  }
  C1=svds(M2,kmax,kmax,kmax)$u
  M1tilde=matrix(0,p1,p1)
  for(t in 1:T){
    M1tilde=M1tilde+Y[t,,]%*%C1%*%t(C1)%*%t(Y[t,,])
  }
  M1tilde=M1tilde/(T*p2)


  beta=log(p1)/log(T*p2)
  delta=epsilon*(beta<0.5|beta==0.5)+(1-1/(2*beta)+epsilon)*(beta>0.5)

  ## calculate phi
  lam1=eigen(M1tilde)$values ##  eigenvalues of M1
  lamhat=lam1/(sum(lam1)/p1)
  phihat=(abs(p1^(-delta)*lamhat))^r


  ## testing the factor number H_0: k1>=r with no projection technique
  Psihat=NULL
  for (s in 1:S){
    eta=rnorm(M/min(k,2))
    psihat1=(sqrt(phihat[k])*eta)<0.7
    psihat2=(sqrt(phihat[k])*eta)<(-0.7)
    psihat3=(sqrt(phihat[k])*eta)<2.4
    psihat4=(sqrt(phihat[k])*eta)<(-2.4)
    nuhat1=2*sum(psihat1-0.5)/sqrt(M)
    nuhat2=2*sum(psihat2-0.5)/sqrt(M)
    nuhat3=2*sum(psihat3-0.5)/sqrt(M)
    nuhat4=2*sum(psihat4-0.5)/sqrt(M)
    Psihat[s]=0.45*nuhat1^2+0.45*nuhat2^2+0.05*nuhat3^2+0.05*nuhat4^2## testing tatistics
  }
  Qhat=1/S*(sum(Psihat<=c1))

  reject_hat=1-1*(Qhat>=(1-alpha-S^(-fq)))

  return(reject_hat)
}



#' determine row factor number - test
#'
#' This function determines the number of row factors
#' under a two-way factor structure, using randomized test method.
#'
#' See He et al. (2023)
#' @references He Y, Kong X, Trapani L, & Yu L (2023).
#' One-way or two-way factor model for matrix sequences? \emph{Journal of Econometrics},
#' 235(2), 1981-2004.
#'
#' @param Y data, a \eqn{T\times p1\times p2} array.
#' @param alpha a number in (0,1), indicating the significance of the test.
#' @param type indicates how to calculate the sample covariance. "flat" for the
#' flat version, while others for the projected version.
#' @param kmax a positive integer smaller than p2, indicating the
#' upper bound for the factor numbers, and the dimension of projection matrix.
#' @param epsilon a small positive number in (0,1), indicating the size of scaling.
#' @param r a positive number indicating the order of the power function
#' for transforming the rescaled eigenvalue.
#' @param M a large integer for the number of Gaussian variables in the randomized test.
#' @param S another large integer for the number of replications in the strong rule. Usually \eqn{M=S=T}.
#' @param fq a number in (0,0.5), controlling the threshold function of the strong rule.
#'
#' @return an integer for the number of row factors. To determine the number of column
#' factors, just transpose the observation matrices.
#' @export
#'
#' @examples
#' k1=3
#' k2=3
#' Sample_T=100
#' p1=40
#' p2=20
#'
#' Y=gen.data(Sample_T,p1,p2,k1,k2,tau=0.5,change=0)
#' KSTP(Y)
#' KSTP(aperm(Y,c(1,3,2)))
#'
KSTP=function(Y,alpha=0.05,type="proj",kmax=4,epsilon=0.05,r=8,M=100,S=100,fq=1/4){

  result=vector()
  reject=0  ## if reject at some step, this will change to 1
  k.hat=kmax ## the intial estimated factor number, equal to kmax
  k=1 ## start of the testing
  if(type=="flat"){
    while(reject==0 & k<kmax+1){
      test.k=ITP_noproj(Y,k,alpha,epsilon,r,M,S,fq)
      if(test.k==1){
        reject=1 ## if reject, stop procedure and return the estimated factor number
        k.hat=k-1 ## estimated factor number
        break
      } else {
        k=k+1 ## if not reject, keep going until kmax
      }
    }
    return(k.hat)
  }else{
    while(reject==0 &k<kmax+1){
      test.k=ITP_proj(Y,k,alpha,kmax,epsilon,r,M,S,fq)
      if(test.k==1){
        reject=1 ## if reject, stop procedure and return the estimated factor number
        k.hat=k-1 ## estimated factor number
        break
      } else {
        k=k+1 ## if not reject, keep going untill kmax
      }
    }
    return(k.hat)
  }
}




#' explanatory power of factors
#'
#' This function calculates the cumulative explanatory power of the
#' leading  row factors, in terms of the explained variance, under a
#' two-way factor structure.
#'
#' @param Y data, a \eqn{T\times p1\times p2} array.
#' @param k a positive integer indicating the number of factors investigated, should be
#'  smaller than p1.
#' @param type indicates how to calculate the sample covariance. "flat" for the
#' flat version, while others for the projected version.
#' @param kmax a positive integer smaller than p2, indicating the
#' upper bound for the factor numbers, and the dimension of projection matrix.
#' @param plot a logical value. When \emph{plot=1}, a figure of the
#' cumulative explanatory power will be plotted, with x axis being the number of
#' factors, and y axis being the cumulative explained variance.
#'
#' @return a vector with k entries, corresponding to the cumulative explanatory power
#' of the leading k factors.
#' @export
#'
#' @examples
#' k1=3
#' k2=3
#' Sample_T=100
#' p1=40
#' p2=20
#'
#' Y=gen.data(Sample_T,p1,p2,k1,k2,tau=0.5,change=0)
#' var.exp(Y,k=5,plot=1)
#'
var.exp<-function(Y,k=2,type="proj",kmax=4,plot=0){
  T=dim(Y)[1]
  p1=dim(Y)[2]
  p2=dim(Y)[3]

  if(type=="flat"){
    M1=matrix(0,p1,p1)
    ## Calculate M1 (scaled)
    for(i in 1:T){
      M1=M1+Y[i,,]%*%t(Y[i,,])
    }
    M1=M1/(T*p2)
    eig=svds(M1,k,k,k)$d
    prop=cumsum(eig)/sum(diag(M1))
    if(plot==1){
      plot(prop,type="b",xlab="number",ylab="proportion")
    }
    return(prop)
  }else{
    ##  eigenvalues of projected method, k2 is set as kmax to estimate C
    M2=matrix(0,p2,p2)
    for(t in 1:T){
      M2=M2+t(Y[t,,])%*%Y[t,,]
    }
    C1=svds(M2,kmax,kmax,kmax)$u
    M1tilde=matrix(0,p1,p1)
    for(t in 1:T){
      M1tilde=M1tilde+Y[t,,]%*%C1%*%t(C1)%*%t(Y[t,,])
    }
    M1tilde=M1tilde/(T*p2)

    eig=svds(M1tilde,k,k,k)$d
    prop=cumsum(eig)/sum(diag(M1tilde))
    if(plot==1){
      plot(prop,type="b",xlab="number",ylab="proportion")
    }
    return(prop)
  }
}




