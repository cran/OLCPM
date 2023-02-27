\name{test.matrix.proj}
\alias{test.matrix.proj}
\title{
test change point for matrix-valued online time series-projected version
}
\description{
This function tests change point for matrix-valued online time series, under a two-way factor structure, using projected sample covariance matrix
}
\usage{
test.matrix.proj(Y, r, m, epsilon, p, kmax, decrease = 0, method = "ps",
  
 eta = 0.25, cv = 2.3860)
}
\arguments{
  \item{Y}{
  data, a \eqn{T\times p_1\times p_2} array.
}
  \item{r}{
  a positive integer indicating which eigenvalue to monitor. \eqn{r=1} for the largest eigenvalue.
}
  \item{m}{
    a positive integer (\eqn{>1}) indicating the bandwidth of the rolling windom.
}
  \item{epsilon}{
    the rescaling parameter taking value in \eqn{(0,1)}; see He et al. (2021).
}
  \item{p}{
    a positive integer indicating the order of the transformation function \eqn{g(x)=|x|^p}. Motivated by He et al. (2021), \eqn{p} should be chosen according to the moments of the data; see more details in He et al. (2021).
}
  \item{kmax}{
     a positive number determining the column number of the projection matrix, should be larger than 0 but smaller than \eqn{p_2}.
}
  \item{decrease}{
    a logical value. If \eqn{\texttt{decrease=1}}, testing the decrease of factor number.
}
   \item{method}{
    indicating the test statistic, ``ps'' for the partial-sum method; others for the worst-case method.
}
  \item{eta}{
    a number between \eqn{[0,1)}, indicating the parameter \eqn{\eta} used in the partial-sum statistic.
}
  \item{cv}{
    critical value, related to the significance level and test statistic. For worst-case statistic, calculate cv by function get.cv. For partial-sum statistic, the default cv only works for \eqn{\eta=0.25} or \eqn{\eta=0.75}. For other \eqn{\eta}, see Horváth et al. (2004).
}
}
\details{
See He et al. (2021).
}
\value{
The return value is a list. In this list, it contains the following:
 \item{test}{a logical value. 1 indicating the existence of change point, 0 indicating no change point.}
 \item{loc}{an integer larger than m, indicating the location of change point; or NA when no change point is reported}
}
\references{
Horváth L, Hušková M, Kokoszka P, et al (2004). Monitoring changes in linear models. Journal of statistical Planning and Inference, 126(1): 225-251.

He Y, Kong X, Trapani L, & Yu L(2021). Online change-point detection for matrix-valued time series with latent two-way factor structure. \emph{arXiv preprint}, arXiv:2112.13479.
}
\author{
Yong He, Xinbing Kong, Lorenzo Trapani, Long Yu
}
\examples{
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


# calculate cv for eta=0.5 and "wc"
Tm=Sample_T-m
cv1=getcv(Tm,0.05,method="ps")
cv2=getcv(Tm,0.05,method="wc")
print(c(cv1,cv2))


## test with Y, projection
re1=test.matrix.proj(Y,k1+1,m,epsilon,p,kmax,0,method="ps",eta=0.25)
print(re1)

re2=test.matrix.proj(Y,k1+1,m,epsilon,p,kmax,0,method="ps",eta=0.75)
print(re2)

re3=test.matrix.proj(Y,k1+1,m,epsilon,p,kmax,0,method="ps",eta=0.5,cv=cv1)
print(re3)

re4=test.matrix.proj(Y,k1+1,m,epsilon,p,kmax,0,method="wc",eta=0.5,cv=cv2)
print(re4)

}
