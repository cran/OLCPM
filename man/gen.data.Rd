\name{gen.data}
\alias{gen.data}
\title{
generate data
}
\description{
This function generates matrix-valued time series under a two-way factor structure with/without a change point.
}
\usage{
gen.data(Sample_T,p1,p2,k1,k2,tau=0.5,change=0,pp=0.3,a=0,cc=0)
}
\arguments{
  \item{Sample_T}{
  positive integer indicating the length of series
}
  \item{p1}{
  positive integer indicating the row dimension
}
  \item{p2}{
  positive integer indicating the column dimension
}
  \item{k1}{
  positive integer indicating the number of row factors
}
  \item{k2}{
  positive integer indicating the number of column factors
}
  \item{tau}{
  a real number in \eqn{(0,1)}, indicating the location of change point, i.e., (\eqn{\tau\times T})
}
  \item{change}{
  the type of change, taking 0 for no change point, taking 1 for the case that the last column of R changes, taking other values for the case that a new row factor occurs.
}
  \item{pp}{
  a number in \eqn{(0,1]}, indicating the magnitude of the break. When \eqn{\texttt{change=1}}, pp is the proportion of entries in the last column of \eqn{R} changing. When change is not equal to 0 or 1, pp is the proportion of non-zero entries in the new factor loading.
}
  \item{a}{
  a number in \eqn{[0,\text{min}(p_1,p_2))}, indicating the cross-sectional correlations of the idiosyncratic errors.
}
  \item{cc}{
  a number in \eqn{[0,1)}, indicating the AR(1) coefficient of the factor and error processes.
}
}
\details{
See the paper He et al. (2021).
}
\value{
The return value is a \eqn{T\times p_1 \times p_2} array.
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

}
