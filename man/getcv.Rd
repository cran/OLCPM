\name{getcv}
\alias{getcv}
\title{
calculate critical values
}
\description{
This function calculates critical values for the partial-sum statistic with \eqn{\eta=0.5} or the the worst-cas statistic.
}
\usage{
getcv(Tm, alpha = 0.05, method = "ps")
}
\arguments{
  \item{Tm}{
  the length of the eigenvalue series, equal to \eqn{T-m}
}
  \item{alpha}{
  significance level of the test, taking value from \eqn{(0,1)}
}
  \item{method}{
  ``ps'' for the partial-sum staistic, others for the worst-case statistic.
}
}
\details{
For the partial-sum statistic, this function only works for \eqn{\eta=0.5}. For other values of \eqn{\eta}, the critical value should be given by the user according to Horvath and Huskova (2004).
}
\value{
The return value is a real number.
}

\references{
Horváth L, Hušková M, Kokoszka P, et al (2004). Monitoring changes in linear models. Journal of statistical Planning and Inference, 126(1): 225-251.

He Y, Kong X, Trapani L, & Yu L(2021). Online change-point detection for matrix-valued time series with latent two-way factor structure. \emph{arXiv preprint}, arXiv:2112.13479.
}
\author{
Yong He, Xinbing Kong, Lorenzo Trapani, Long Yu
}
\examples{
cv1=getcv(160,0.05,method="ps")
cv2=getcv(160,0.05,method="wc")
print(c(cv1,cv2))
}
