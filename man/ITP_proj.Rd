% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function.R
\name{ITP_proj}
\alias{ITP_proj}
\title{testing the number of row factors- with projection}
\usage{
ITP_proj(
  Y,
  k = 1,
  alpha = 0.05,
  kmax = 4,
  epsilon = 0.05,
  r = 8,
  M = 100,
  S = 100,
  fq = 1/4
)
}
\arguments{
\item{Y}{data, a \eqn{T\times p1\times p2} array.}

\item{k}{an positive integer indicating which eigenvalue to test.}

\item{alpha}{a number in (0,1), indicating the significance of the test.}

\item{kmax}{a positive integer smaller than p2, indicating the
upper bound for the factor numbers, and the dimension of projection matrix.}

\item{epsilon}{a small positive number in (0,1), indicating the size of scaling.}

\item{r}{a positive number indicating the order of the power function
for transforming the rescaled eigenvalue.}

\item{M}{a large integer for the number of Gaussian variables in the randomized test.}

\item{S}{another large integer for the number of replications in the strong rule. Usually \eqn{M=S=T}.}

\item{fq}{a number in (0,0.5), controlling the threshold function of the strong rule.}
}
\value{
a logical value. 1 for "the number of row factors is smaller than k".
0 for "at least k row factors exists".
}
\description{
This function tests whether the number of row factors is equal or larger
than a given integer, under a two-way factor model, using projected version
of sample covariance.
}
\details{
See He et al. (2023)
}
\examples{
k1=3
k2=3
Sample_T=100
p1=40
p2=20

Y=gen.data(Sample_T,p1,p2,k1,k2,tau=0.5,change=0)
ITP_proj(Y,k=1,M=Sample_T,S=Sample_T)
ITP_proj(Y,k=4,M=Sample_T,S=Sample_T)

}
\references{
He Y, Kong X, Trapani L, & Yu L (2023).
One-way or two-way factor model for matrix sequences? \emph{Journal of Econometrics},
235(2), 1981-2004.
}
