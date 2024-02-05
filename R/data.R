#' cv.table
#'
#' A data matrix containing approximated critical values for the scaled
#' Wiener process. A supplement to Table 1 in Horvath et al. (2004). The
#' rows are for different gamma on a grid of 50 points in [0,0.49) with step
#' size equal to 0.01, while the columns are for different alpha on a grid
#' of 50 points in [0.001,0.5) with step size equal to 0.001. All the values are
#' based on 100000 replications of simulated Wiener process on an equal grid of
#' 20000 points in [0,1], thus can be slightly different from those in Horvath et al. (2004).
#'
#'#' @references Horvath L, Huskova M, Kokoszka P, et al (2004). Monitoring
#' changes in linear models. \emph{Journal of statistical Planning and
#' Inference}, 126(1): 225-251.
#'
#' @format A matrix of size 50*500.
#' \describe{
#'   \item{gamma}{the scaling parameter}
#'   \item{alpha}{significance level (1-alpha)-th quantile}
#' }
"cv.table"
