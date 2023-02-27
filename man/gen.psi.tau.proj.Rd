\name{gen.psi.tau.proj}
\alias{gen.psi.tau.proj}
\title{
calculate eigenvalue series by projected method
}
\description{
This function calculates the rolling eigenvalue series for the monitoring process, based on the projected version of sample covanriance matrix.
}
\usage{
gen.psi.tau.proj(Y, r, m, delta, p, kmax)
}
\arguments{
   \item{Y}{
  the observed \eqn{T\times p_1\times p_2} array. \eqn{T} is the sample size, \eqn{p_1} and \eqn{p_2} are the row and column dimensions, respectively.
}
  \item{r}{
  a positive integer determining which eigenvalue to monitor. \eqn{r=1} for the largest eigenvalue.
}
  \item{m}{
  a positive integer (\eqn{>1}) indicating the bandwidth of the rolling windom.
}
  \item{delta}{
  a number in \eqn{(0,1)} indicating the rescaling parameter for the eigenvalue. The default approach to calcualte delta is in the paper He et al. (2021).
}
  \item{p}{
  a positive integer indicating the order of the transformation function \eqn{g(x)=|x|^p}. Motivated by the paper, \eqn{p} should be chosen according to the moments of the data; see more details in He et al. (2021).
}
  \item{kmax}{
   a positive integer indicating the column number of the projection matrix, should be larger than 0 but smaller than \eqn{p_2}.
}
}
\details{
The rolling eigenvalue series will start at the stage \eqn{m+1}, with length \eqn{T-m}.
}
\value{
When \eqn{T>m}, the return value is a \eqn{(T-m)\times 3} matrix, containing the original, rescaled, and the transformed eigenvalue series. Otherwise, the function outputs an error.
}
\references{
He Y, Kong X, Trapani L, & Yu L(2021). Online change-point detection for matrix-valued time series with latent two-way factor structure. \emph{arXiv preprint}, arXiv:2112.13479.
}
\author{
Yong He, Xinbing Kong, Lorenzo Trapani, Long Yu
}

\examples{
## generate data
k1=3
k2=3
epsilon=0.05
Sample_T=50
p1=40
p2=20
kmax=8
p=8
m=p2

# generate data
Y=gen.data(Sample_T,p1,p2,k1,k2,tau=0.5,change=1,pp=0.3)

# calculate delta
temp=log(p1)/log(m*p2)
delta=epsilon*(temp<=0.5)+(epsilon+1-1/(2*temp))*(temp>0.5)

# calculate psi.tau
psi2=gen.psi.tau.proj(Y,k1+1,m,delta,p,kmax)

}

