#' Monte Carlo significance test for any clustering method
#'

test_gfl_approx <- function(y=y_vec,
                            Dmat=Dmat,
                            v=contrast,
                            K=K_init,
                            K_CC=CC_init,
                            stop_criteria = stop_criteria,
                            sigma=sigma,
                            two_sided=TRUE,
                            segment_list=CC_to_test,
                            ndraws=1000,
                            random_seed=2021){
  set.seed(random_seed)
  nu_T_y <- sum(contrast*y_vec)
  nu_norm_sq <- sum(contrast*contrast)
  # sampled proposal distribution
  phi_list <- stats::rnorm(ndraws, mean = nu_T_y, sd = sqrt(scale_factor))


  for(j in 1:ndraws){
    phi_j <- phi_list[j]

    test_1 <- check_segment_GFL(y, v, phi,
                                  Dmat, K, sigma,
                                  c_1 = sum(v*y),
                                  c_2 = sum(v*v),
                                  old_cluster,
                                  comparison,
                                  stop_criteria,
                                  K_CC = NULL,
                                  K_lambda = NULL,
                                  segment_list = NULL)


  }


}
