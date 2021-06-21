#' Summarize the  result for powerful graph fused lasso inference
#' @param object output from running spike_estimates
#' @param  ... to be passed to methods
#' @return A data frame with summarized results
#' @examples
#' lev1 <- 0 # mean for group 1
#' lev2 <- 3 # mean (absolute value) for group 2/3
#' sigma <- 1 # level of noise
#' nn <- 8 # grid size
#' Dmat <- graph2D_Dmat(nn^2) # generate D matrix for the 2D fused lasso
#' ### Create the underlying signal
#' A <- matrix(lev1, ncol=nn, nrow = nn)
#' A[1:round(nn/3),1:round(nn/3)] <- 1*lev2
#' A[(nn-2):(nn),(nn-2):(nn)] <- -1*lev2
#' set.seed(2005)
#' A.noisy <- A + rnorm(nn^2,mean=0,sd=sigma)
#' y <- do.call(c,lapply(1:nrow(A.noisy),function(irow)A.noisy[irow,]))
#' ### Run a test for a difference in means between estimated clusters 1 and 2
#' result_demo <- fusedlasso_inf(y=y, D=Dmat, c1=1, c2=2, method="K", sigma=sigma, K=13)
#' summary(result_demo)
#' @export
summary.fusedlasso_inf <- function(object, ...){
  result <- data.frame(cluster_1 = object$c1,
                       cluster_2 = object$c2,
                       diff_in_means = object$test_stats,
                       pval_c1c2 = object$Union,
                       pval_hyun = object$Hyun,
                       LCB = ifelse(is.null(object$CI_result[1]),NA,object$CI_result[1]),
                       UCB = ifelse(is.null(object$CI_result[2]),NA,object$CI_result[2]))
  return(result)
}
