# ----- main function to test equality of the means of two estimated connected components -----

#' More powerful test for  the graph fused lasso
#'
#' This functions tests the null hypothesis of no difference in means between
#' connected components \code{c1} and \code{c2} of the output of the graph fused
#' lasso solution. The ordering are numbered as per the results of the \code{fusedlasso}
#' function in the \code{genlasso} package.
#'
#' Currently, we support two different conditioning sets: conditioning set 1 is based on
#' the output after K steps dual-path algorithm; and conditioning set 2 is based on the
#' output of the after the dual-path algorithm yields c connected components in the output.
#'
#' Input:
#' @param y Numeric vector; \eqn{n} dimensional observed data
#' @param D Numeric matrix; \eqn{m} by \eqn{n} penalty matrix, i.e.,
#' the oriented incidence matrix over the underlying graph
#' @param c1,c2 Integers selecting the two connected components to test, as indexed by the results of
#' \code{genlasso::fusedlasso}.
#' @param method One of "K" or "CC", which indicates which conditioning set to use
#' @param sigma Numeric; noise standard deviation for the observed data, a non-negative number.
#' @param K Integer; number of steps to run the dual-path algorithm.
#' It must be specified if method=="K".
#' @param L Integer; the targeted number of connected components.
#' It must be specified if method=="CC".
#' @param early_stop Numeric; specify when the truncation set computation
#' should be terminated. The default is NULL, which indicates infinity.
#' @param compute_ci Logical; the default is False. Specifying whether confidence intervals for \eqn{\nu^{T}\beta}, the
#' difference in means between the two estimated connected components, should be computed.
#' @param alpha_level Numeric; parameter for the 1-\code{alpha_level} confidence interval, defeault to 0.05
#' @return Returns a list with elements:
#' \itemize{
#' \item \code{pval} the p-value in Chen et al. (2021+)
#' \item \code{truncation_set} the conditioning set of Chen et al. (2021+) stored as \code{Intervals} class
#' \item \code{test_stats} test statistics: the difference in means of two connected components
#' \item \code{beta_hat} Graph fused lasso estimates
#' \item \code{connected_comp} Estimated connected component
#' \item \code{Naive} the naive p-value using a z-test
#' \item \code{Hyun} the p-value proposed in Hyun et al. (2018)
#' \item \code{hyun_set} the conditioning set of  Hyun et al. (2018) stored as \code{Intervals} class
#' \item \code{CI_result} confidence interval of level 1-\code{alpha_level} if \code{compute_ci=TRUE}
#' }
#' @export
#'
#' @details
#' Consider the generative model \eqn{Y_j = \beta_j + \epsilon_j, \epsilon_j \sim N(0, \sigma^2). j=1,...,n}, where
#' the underlying signal \eqn{\beta} is assumed to be piecewise constant with respect to an underlying
#' graph. The fused lasso estimate minimizes the following objective function
#' \deqn{minimize_{\beta} \frac{1}{2} \sum_{j=1}^{n} ( y_j - \beta_j )^2 + \lambda \sum_{(i,j)\in E}|\beta_i-\beta_j|,}
#' where E is the edge set of the underlying graph. The solution \eqn{\hat{\beta}} can then be
#' segment into connected components; that is, the set of \eqn{\hat{\beta}} that takes on the
#' same value, and are connected in the original graph.
#'
#' Now suppose we want to test whether the means of two estimated connected components \code{c1} and \code{c2}
#' are equal; or equivalently, the null hypothesis of the form \eqn{H_{0}:  \nu^T \beta = 0} versus
#' \eqn{H_{1}:  \nu^T \beta \neq 0} for suitably chosen \eqn{\nu}.
#'
#' This function computes the following p-value:
#' \deqn{P(|\nu^T Y| \geq |\nu^T y| \; | \;  \hat{C}_1, \hat{C}_2 \in CC_K(Y),  \Pi_\nu^\perp Y  =  \Pi_\nu^\perp y)},
#' where \eqn{CC_K(Y)} is the set of estimated connected components from applying K steps of the dual path algorithm on data Y
#' , and \eqn{\Pi_\nu^\perp} is the orthogonal projection to the orthogonal complement of \eqn{\nu}.
#' In particular, the test based on this p-value controls the selective Type I error and has higher power than an existing method
#' by Hyun et al. (2018). Readers can refer to the Section 3 in Chen et al. (2021+) for more details.
#'
#' @examples
#' lev1 <- 0 # mean for group 1
#' lev2 <- 3 # mean (absolute value) for group 2/3
#' sigma <- 1 # level of noise
#' nn <- 8 # grid size
#' Dmat <- genlasso::getD2d(nn, nn) # generate D matrix for the 2D fused lasso
#' ### Create the underlying signal
#' A <- matrix(lev1, ncol=nn, nrow = nn)
#' A[1:round(nn/3),1:round(nn/3)] <- 1*lev2
#' A[(nn-2):(nn),(nn-2):(nn)] <- -1*lev2
#' ### Visualize the underlying signal
#' lattice::levelplot(A)
#' set.seed(2005)
#' A.noisy <- A + rnorm(nn^2,mean=0,sd=sigma)
#' y <- c(t(A.noisy))
#' ### Now use the fusedlasso function to obtain estimated connected components after K=13
#' ### steps of the dual path algorithm
#' K = 13
#' complete_sol <- genlasso::fusedlasso(y=y,D=Dmat,maxsteps=K)
#' beta_hat <- complete_sol$beta[,K]
#' ### estimated connected components
#' estimated_CC <- complete_sol$pathobjs$i
#' estimated_CC
#' ### Run a test for a difference in means between estimated connected components 1 and 2
#' result_demo <- fusedlasso_inf(y=y, D=Dmat, c1=1, c2=2, method="K", sigma=sigma, K=K)
#' summary(result_demo)
#'
#' @references
#' Chen YT, Jewell SW, Witten DM. (2021+) More powerful selective inference for the graph fused lasso
#' Hyun S, G’Sell M, Tibshirani RJ. (2018) Exact post-selection inference for the generalized lasso path. Electron J Stat.

fusedlasso_inf <- function(y, D, c1, c2, method, sigma, K=NULL, L=NULL, early_stop=NULL,
                           compute_ci = FALSE, alpha_level = 0.05){

  if(!method%in%c("K","CC")){stop("Method must be 'K' or 'CC'.")}
  if((method=="K")&is.null(K)){stop("Must specify K, the steps of the dual-path algorithm.")}
  if((method=="CC")&is.null(c)){stop("Must specify c, the number of connected components.")}
  if(sigma<=0){stop("Positive sigma required!")}
  # verify the dimension of y and D
  if(dim(D)[2]!=length(y)){stop("ncol(D)!=length(y), check your input!")}

  v <- rep(0, length(y))
  fused_lasso_sol <- dualpathFused_CC_indexed(
                                 y=y,
                                 D=D,
                                 v=rep(1,ncol(D)),
                                 sigma=sigma,
                                 verbose=FALSE,
                                 maxsteps = K,
                                 K_CC = L,
                                 stop_criteria = method)

  fused_lasso_mem <- fused_lasso_sol$pathobjs$membership
  n_clusters <- length(unique(fused_lasso_mem))


  if(max(c1,c2)>n_clusters){
    stop(paste0("Only ",n_clusters," connected components are estimated,
                double-check your c1 and c2 specification!"))
  }


  segment_list <- list(which(fused_lasso_mem==c1),
                       which(fused_lasso_mem==c2))
  v[segment_list[[1]]] <- 1/length(segment_list[[1]])
  v[segment_list[[2]]] <- -1/length(segment_list[[2]])

  search_result <- ComputeUnionIntervals_GFL(y=y,
                            Dmat=D,
                            v=v,
                            K=K,
                            sigma=sigma,
                            stop_criteria=method,
                            K_CC = L,
                            eta=1e-4,
                            two_sided = TRUE,
                            early_stop = early_stop,
                            segment_list = segment_list)


  names(search_result) <- c('Naive',
                            'Hyun',
                            'Union',
                            'sum_intervals',
                            'truncation_set',
                            'test_stats',
                            'sd',
                            'hyun_set',
                            "two_sided")
  search_result$pval <- search_result$Union

  K_beta_hat <- ncol(fused_lasso_sol$beta)+1

  fused_lasso_sol_output <- dualpathFused_CC_indexed(
    y=y,
    D=D,
    v=v,
    sigma=sigma,
    verbose=FALSE,
    maxsteps=K_beta_hat,
    stop_criteria = "K")

  all_beta_hat <- fused_lasso_sol_output$beta

  if(compute_ci){
    CI_result <- compute_CI(search_result$test_stats, (search_result$sd)^2, 1,
               truncation=search_result$truncation_set,
               alpha=alpha_level, steps_lim=50)
  }else{
    CI_result <- NULL
  }


  return_result <- c(search_result, list(connected_comp = fused_lasso_mem,
                                         CI_result = CI_result,
                                         c1 = c1,
                                         c2 = c2,
                                         call = match.call(),
                                         beta_hat = all_beta_hat[,ncol(all_beta_hat)]))


  class(return_result) <- "fusedlasso_inf"
  return(return_result)

  }





