\name{GFLassoInference}
\alias{GFLassoInference}
\docType{package}
\title{More powerful test for the graph fused lasso}

\description{
GFLassoInference is an R package for testing a difference in means between a pair of connected components estimated from the graph fused lasso. Details can be found in our manuscript (Chen et al. 2021+). In addition to the manual in R, a detailed tutorial can be found at https://yiqunchen.github.io/GFLassoInference/.}

\details{
Maintainer: Yiqun T. Chen <yiqunc@uw.edu>
The main function in the package is \code{\link{fusedlasso_inf}}. In our tutorials and experiments, we also used functions from the \code{genlasso} package and the \code{genlassoinf} package.
}

\examples{
lev1 <- 0 # mean for group 1
lev2 <- 3 # mean (absolute value) for group 2/3
sigma <- 1 # level of noise
nn <- 8 # grid size
Dmat <- genlasso::getD2d(nn, nn) # generate D matrix for the 2D fused lasso
### Create the underlying signal
A <- matrix(lev1, ncol=nn, nrow = nn)
A[1:round(nn/3),1:round(nn/3)] <- 1*lev2
A[(nn-2):(nn),(nn-2):(nn)] <- -1*lev2
### Visualize the underlying signal
lattice::levelplot(A)
set.seed(2005)
A.noisy <- A + rnorm(nn^2,mean=0,sd=sigma)
y <- c(t(A.noisy))
### Now use the fusedlasso function to obtain estimated connected components after K=13
### steps of the dual path algorithm
K = 13
complete_sol <- genlasso::fusedlasso(y=y,D=Dmat,maxsteps=K)
beta_hat <- complete_sol$beta[,K]
### estimated connected components
estimated_CC <- complete_sol$pathobjs$i
estimated_CC
### Run a test for a difference in means between estimated connected components 1 and 2
result_demo <- fusedlasso_inf(y=y, D=Dmat, c1=1, c2=2, method="K", sigma=sigma, K=K)
summary(result_demo)
}

\references{
Chen YT, Jewell SW, Witten DM. (2021+) More powerful selective inference for the graph fused lasso

Hyun S, G’Sell M, Tibshirani RJ. (2018) Exact post-selection inference for the generalized lasso path. Electron J Stat.

Tibshirani RJ, Taylor J. The solution path of the generalized lasso. Ann Stat. 2011;39(3):1335-1371. doi:10.1214/11-AOS878
}
