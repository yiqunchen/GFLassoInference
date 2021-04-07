#' Monte Carlo significance test for GFLasso, two sided-alternatives only
#' @export
test_gfl_approx <- function(y,
                            Dmat,
                            v,
                            K_init,
                            sigma,
                            segment_list,
                            ndraws=1000,
                            stop_criteria = "K",
                            random_seed=2021){

  if(stop_criteria!="K"){
    stop("Only fixed K step inference is implemented!")
  }
  set.seed(random_seed)
  nu_T_y <- sum(v*y_vec)
  nu_norm_sq <- sum(v*v)
  scale_factor = nu_norm_sq*sigma*sigma
  # sampled proposal distribution
  phi_list <- stats::rnorm(ndraws, mean = nu_T_y, sd = sqrt(scale_factor))
  # log-sum-exp
  log_p_x_i <- stats::dnorm(phi_list, mean = 0, sd = sqrt(scale_factor), log = TRUE)
  log_q_x_i <- stats::dnorm(phi_list, mean = nu_T_y, sd = sqrt(scale_factor),log = TRUE)
  log_survives <- log_p_x_i-log_q_x_i
  # ratio_list <- p_x_i/q_x_i
  # cond 1 for the numerator
  num_condition_1 <- (abs(phi_list)>=abs(nu_T_y))
  segment_in_cluster <- rep(0, length(ndraws))

  identity_design <- diag(x=1,nrow=length(y_vec),ncol=length(y_vec))

  for(j in 1:ndraws){

    phi_j <- phi_list[j]
    y_phi_j <- as.numeric(y-(sum(v*y)-phi_j)*v/sum(v*v))

    if(stop_criteria=="K"){
      y_hat_fused_lasso <- genlasso::fusedlasso(y=y_phi_j, X=identity_design,D=Dmat, maxsteps=K_init+1)
    }

    new_cluster <- as.factor(y_hat_fused_lasso$beta[,K_init])
    segment_in_cluster[j] <- PGInference::check_segment_in_model(segment_list, new_cluster)

  }

  if(sum(segment_in_cluster)==0){
    warning("Oops - we didn't generate any samples that preserved the connected components!
            Try re-running with a larger value of ndraws.")
    return(list(stat=abs(nu_T_y), pval=NA, stderr=NA, connected_components=CC_to_test))
  }
  #  Approximate p-values
  # first truncate to only the phis that preserve the connected components
  log_survives <- log_survives[segment_in_cluster==TRUE]
  phi_list <- phi_list[segment_in_cluster==TRUE]
  # log-sum-exp to compute the p-val
  log_survives_shift <- log_survives - max(log_survives)
  props <- exp(log_survives_shift)/sum(exp(log_survives_shift))
  pval <- sum(props[abs(phi_list) >= abs(nu_T_y)])
  # variance
  var_pval <- (1 - pval)^2*sum(props[abs(phi_list) >= abs(nu_T_y)]^2) +
    pval^2*sum(props[abs(phi_list) < abs(nu_T_y)]^2)

  return(list(stat=abs(nu_T_y), pval=pval, stderr=sqrt(var_pval), connected_components=CC_to_test))

}
