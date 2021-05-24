### This function computes iteratively the partition of intervals
### selected by the Gamma matrix with perturbed y
### The result is a union of intervals with the same selected
### components as y (observed value)
### Input:
### @y: observed data
### @v: contrast vector
### @f0: object
### @end_tolerance: allowed tolerance in the difference between two
### adjacent intervals; default to be 1e-6
### @eta: proposed length of perturbation for new y
### @max_union_size: how many intervals to initialize on each end
### Output:
### truncation_list: a list of intervals (each item in the form (vlo, vup))
#' @export
ComputeUnionIntervals_GFL <- function(y,
                                      Dmat,
                                      v,
                                      K,
                                      sigma,
                                      stop_criteria,
                                      K_CC = NULL,
                                      K_lambda = NULL,
                                      eta=1e-4,
                                      comparison = 'CC',
                                      two_sided = TRUE,
                                      end_tolerance=1e-6,
                                      t=0.5,
                                      max_union_size = 100000,
                                      early_stop = sqrt(sum(v*v*sigma))*20,
                                      segment_list = NULL,
                                      complete_k_step = TRUE,
                                      intervals_delta = .Machine[[ "double.eps" ]]^0.5){

  if (comparison!='CC'){
    stop('Comparisons other than CC is not available!\n')
  }

  # initial choice of k matters
  f0 <- dualpathFused_CC_indexed(y,Dmat,v = v,
                                 sigma=sigma,
                                 verbose=FALSE,
                                 maxsteps = K,
                                 K_CC = K_CC,
                                 K_lambda = K_lambda,
                                 stop_criteria = stop_criteria)

  K_cc <- f0$pathobjs$q
  K_lambda <- f0$lambda[length(f0$lambda)]
  K <- length(f0$action)
  # we never really need the old graph?
  cluster_0 <- f0$pathobjs$membership

  temp_vlo <- unlist(f0$LS_list[names(f0$LS_list)=='vlo'])
  temp_vup <- unlist(f0$LS_list[names(f0$LS_list)=='vup'])

  initail_line_seg <- list(max(temp_vlo),
                           min(temp_vup),
                           max(unlist(f0$LS_list[names(f0$LS_list)=='z'])),
                           max(unlist(f0$LS_list[names(f0$LS_list)=='sd'])))

  names(initail_line_seg) <- c('vlo','vup','z','sd')
  cat(initail_line_seg$z,'z',initail_line_seg$sd,'sd','\n')


  if(two_sided){
    naive_p_val <- naive.two.sided.pval(z = initail_line_seg$z,
                                        mean = 0,
                                        sd = initail_line_seg$sd)
    hyun_p_val <- two.sided.tnorm.pval(z = initail_line_seg$z,
                             mean = 0,
                             sd = initail_line_seg$sd,
                             a = initail_line_seg$vlo,
                             b = initail_line_seg$vup,
                             bits = NULL)

  }else{
    naive_p_val <- naive.one.sided.pval(z = initail_line_seg$z,
                                        mean = 0,
                                        sd = initail_line_seg$sd)
    hyun_p_val <- tnorm.surv(z = initail_line_seg$z,
                             mean = 0,
                             sd = initail_line_seg$sd,
                             a = initail_line_seg$vlo,
                             b = initail_line_seg$vup,
                             bits = NULL)
  }
  #cat(naive_p_val,'pNaive, 2 sided','\n')
  #cat(hyun_p_val,'pHyun, 2 sided','\n')


  line_search_results <- GetUnionIntervals_GFL(y, v,
                                               eta, t,
                                               Dmat, K = K,
                                               K_CC = K_CC,
                                               K_lambda = K_lambda,
                                               stop_criteria = stop_criteria,
                                               sigma = sigma,
                                               cluster_0 = cluster_0,
                                               end_tolerance = end_tolerance,
                                               max_union_size = max_union_size,
                                               initail_line_seg = initail_line_seg,
                                               comparison = comparison,
                                               early_stop = early_stop,
                                               segment_list = segment_list,
                                               complete_k_step = complete_k_step)

  # post-processing the intervals
  # 1. select the relevant intervals for the positive direction
  pos_dir_union <- line_search_results[[2]][which(line_search_results[[1]])]
  # 2. format the union correctly
  pos_dir_union <- lapply(pos_dir_union,function(x)c(min(x),max(x))) # numerical stability for intervals
  pos_intervals <- intervals::Intervals(
    matrix(unlist(pos_dir_union), ncol = 2, byrow = TRUE))
  # 3. take the union - up to machine epsilon accuracy
  pos_intervals_union <- intervals::interval_union(intervals::expand(pos_intervals,
                                                                     intervals_delta,
                                                                     type="relative"))
  # 4. do the same for the negative direction
  neg_dir_union <- line_search_results[[4]][which(line_search_results[[3]])]
  neg_dir_union <- lapply(neg_dir_union,function(x)c(min(x),max(x))) # numerical stability for intervals
  # 2. format the union correctly
  neg_intervals <- intervals::Intervals(
    matrix(unlist(neg_dir_union), ncol = 2, byrow = TRUE))
  # 3. take the union - up to machine epsilon accuracy
  neg_intervals_union <- intervals::interval_union(intervals::expand(neg_intervals,
                                                                     intervals_delta,
                                                                     type="relative"))

  truncation_set <- intervals::interval_union(pos_intervals_union,
                                              neg_intervals_union)

  num_neg_segments <- sum(unlist(lapply(line_search_results[[4]],
                                        function(x)!is.null(x))))

  num_pos_segments <- sum(unlist(lapply(line_search_results[[2]],
                                        function(x)!is.null(x))))

  sum_intervals <- (num_neg_segments+num_pos_segments) #just auxiliary info
  test_stats <- sum(v*y)

  # formatting truncation set so that we can use safe pvals
  truncation_set_df <- data.frame(truncation_set)
  colnames(truncation_set_df) <- c("min_mean","max_mean")
  truncation_set_df$contained <- 1
  # end of formatting

  union_cond_p_val <- calc_p_value_safer(truncation_set_df,
                                         test_stats, sum(v*v),
                                         sigma^2, mu = 0, two_sided=two_sided)

  #cat(union_cond_p_val,'p Union, 2-sided', '\n')

  # hyun_set = list(c(initail_line_seg$vlo,
  #                   initail_line_seg$vup))

  hyun_set <- intervals::Intervals(matrix(c(initail_line_seg$vlo,
                                   initail_line_seg$vup), ncol = 2))

  p_val_list <- list(naive_p_val,
                     hyun_p_val,
                     union_cond_p_val,
                     sum_intervals,
                     truncation_set,
                     test_stats,
                     sigma*sqrt(sum(v*v)),
                     hyun_set,
                     two_sided)

  names(p_val_list) <- c('Naive','Hyun','Union','sum_intervals',
                         'truncation_set','test_stats','sd',
                         'hyun_set',
                         "two_sided")
  return(p_val_list)

}


### This function computes iteratively the partition of intervals
### selected by the Gamma matrix with perturbed y
### The result is a union of intervals with the same selected
### components as y (observed value)

GetUnionIntervals_GFL <- function(y, v, eta,
                                  t,Dmat,
                                  K,K_CC = NULL,
                                  K_lambda = NULL,
                                  stop_criteria,
                                  sigma,
                                  cluster_0,
                                  end_tolerance,
                                  max_union_size,
                                  initail_line_seg,
                                  comparison,
                                  early_stop=NULL,
                                  segment_list = NULL,
                                  complete_k_step = TRUE){

  ### TODO: better way to store CC info
  neg_dir_union <- vector('list', max_union_size)
  pos_dir_union <- vector('list', max_union_size)
  pos_dir_graph_cc <- vector('logical', max_union_size)
  neg_dir_graph_cc <- vector('logical', max_union_size)
  pos_step_taken  <- vector('numeric', max_union_size)
  neg_step_taken <- vector('numeric', max_union_size)
  pos_dir_graph_cc[1] <- TRUE
  neg_dir_graph_cc[1] <- TRUE
  pos_dir_union[[1]] <- c(initail_line_seg$vlo, initail_line_seg$vup)
  neg_dir_union[[1]] <- c(initail_line_seg$vlo, initail_line_seg$vup)

  pos_counter <- 1
  neg_counter <- 1

  c_1 <- sum(v*y)
  c_2 <- sum(v*v)

  pos_upper_limit <- initail_line_seg$vup
  neg_lower_limit <- initail_line_seg$vlo

  # large number to prevent overflow
  if(is.null(early_stop)){
    is.practically.finite <- function(x, upper_bound = 2147483647){
      return(is.finite(x)&(abs(x)<=upper_bound))
    }

  }else{
    is.practically.finite <- function(x, upper_bound = early_stop){
      return(is.finite(x)&(abs(x)<=upper_bound))
    }
  }

  # we first start by looking at the positive direction
  while(is.practically.finite(pos_upper_limit)&(pos_counter<max_union_size)){
    # first proposing
    eta_increase <- eta
    phi <- pos_upper_limit+eta_increase
    temp_seg <- check_segment_GFL(y, v,
                                  phi = phi,
                                  Dmat, K,
                                  sigma,
                                  K_CC = K_CC,
                                  K_lambda = K_lambda,
                                  stop_criteria = stop_criteria,
                                  old_cluster = cluster_0,
                                  comparison = comparison,
                                  segment_list = segment_list)

    # if we overshoot at the left endpoint -> make the increase smaller
    # index 1 is vlo and 2 is vup TODO: named list
    # no abs needed -> note that for early stopping it's possible to have vlo extended
    # beyind last vup

    ## abs don't matter
    while((temp_seg[[1]]$vlo-pos_dir_union[[pos_counter]][2])>=end_tolerance){
      eta_increase <- eta_increase*t
      if(eta_increase<=1e-6){
        cat('current phi', phi,'\n')
        print('eta too small, not gonna work,\n')
        break
      }
      phi <- pos_upper_limit+eta_increase
      temp_seg <- check_segment_GFL(y, v, phi = phi,
                                    Dmat, K, sigma,
                                    K_CC = K_CC,
                                    K_lambda = K_lambda,
                                    stop_criteria = stop_criteria,
                                    old_cluster = cluster_0,
                                    comparison = comparison,
                                    segment_list = segment_list)
    }

    # update stopping criteria

    pos_upper_limit <- temp_seg[[1]]$vup
    pos_counter <- pos_counter+1
    # store data into the pre-allocated lists
    pos_dir_graph_cc[pos_counter] <- temp_seg[[2]]
    pos_step_taken[pos_counter] <- temp_seg[[3]]
    pos_dir_union[[pos_counter]] <- c(temp_seg[[1]]$vlo, temp_seg[[1]]$vup)
    cat('test sanity check, pos ', pos_dir_union[[pos_counter]],'\n')
  }

  # we then proceed to the negative direction
  while(is.practically.finite(neg_lower_limit)&(neg_counter<max_union_size)){

    eta_increase <- eta

    # first proposing eta decrease
    phi <- neg_lower_limit-eta_increase
    temp_seg <- check_segment_GFL(y, v, phi = phi, Dmat,
                                  K, sigma,
                                  K_CC = K_CC,
                                  K_lambda = K_lambda,
                                  stop_criteria = stop_criteria,
                                  old_cluster = cluster_0,
                                  comparison = comparison,
                                  segment_list = segment_list)

    # if we overshoot at the left endpoint -> make the increase smaller
    # index 1 is vlo and 2 is vup TODO: named list
    while(abs(temp_seg[[1]]$vup-neg_dir_union[[neg_counter]][1])>=end_tolerance){
      eta_increase <- eta_increase*t
      if(eta_increase<=1e-6){
        cat('current phi', phi,'\n')
        print('eta too small, not gonna work,\n')
        break
      }
      phi <- neg_lower_limit-eta_increase
      temp_seg <- check_segment_GFL(y, v, phi = phi, Dmat,
                                    K, sigma,
                                    K_CC = K_CC,
                                    K_lambda = K_lambda,
                                    stop_criteria = stop_criteria,
                                    old_cluster = cluster_0,
                                    comparison = comparison,
                                    segment_list = segment_list)
    }
    # update stopping criteria
    neg_lower_limit <- temp_seg[[1]]$vlo
    #cat('neg_lower_limit',neg_lower_limit,'\n')
    neg_counter <- neg_counter+1
    # store data into the pre-allocated lists
    neg_dir_graph_cc[neg_counter] <- temp_seg[[2]]
    neg_step_taken[neg_counter] <- temp_seg[[3]]
    neg_dir_union[[neg_counter]] <- c(temp_seg[[1]]$vlo, temp_seg[[1]]$vup)
    cat('test sanity check, neg ', neg_dir_union[[neg_counter]],'\n')
  }

  if(!is.null(early_stop)){
    neg_dir_union[[neg_counter+1]] <- c(-Inf,neg_dir_union[[neg_counter]][[1]])
    pos_dir_union[[pos_counter+1]] <- c(pos_dir_union[[pos_counter]][[2]],Inf)
    neg_dir_graph_cc[neg_counter+1] <- FALSE
    pos_dir_graph_cc[pos_counter+1] <- FALSE
  }
  list_to_return <- list(pos_dir_graph_cc,pos_dir_union,neg_dir_graph_cc,neg_dir_union,pos_step_taken,neg_step_taken)

  names(list_to_return) <- c('pos_dir_graph_cc','pos_dir_union',
                             'neg_dir_graph_cc','neg_dir_union',
                             'pos_step_taken','neg_step_taken')
  return(list_to_return)
}


### This function computes S_hyun for the perturbation phi
### and decides whether if it should be included in the final set
check_segment_GFL <- function(y, v, phi,
                              Dmat, K, sigma,
                              c_1 = sum(v*y),
                              c_2 = sum(v*v),
                              old_cluster,
                              comparison,
                              stop_criteria,
                              K_CC = NULL,
                              K_lambda = NULL,
                              segment_list = NULL){

  #### y_obs(phi) = y_obs - (v^Ty_obs-phi)*v/||v||_2^2
  # upper bound set up be the max int for now
  new_y <- as.numeric(y-(sum(v*y)-phi)*v/sum(v*v)) # y'(phi)

  new_f0 <- dualpathFused_CC_indexed(new_y,
                                     Dmat, v=v,
                                     sigma=sigma,
                                     K_CC = K_CC,
                                     K_lambda = K_lambda,
                                     stop_criteria = stop_criteria,
                                     verbose=FALSE,
                                     maxsteps = K,
                                     segment_list = NULL)

  if(comparison == 'CC'){
    graph_equal <- check_segment_in_model(segment_list,
                                          new_f0$pathobjs$membership)
  }

  steps_taken <- max(new_f0$pathobjs$k,length(new_f0$action))
  # side info, not critical for inference

  temp_vlo <- unlist(new_f0$LS_list[names(new_f0$LS_list)=='vlo'])
  temp_vup <- unlist(new_f0$LS_list[names(new_f0$LS_list)=='vup'])

  curr_line_seg <- list(max(temp_vlo),
                        min(temp_vup))
  names(curr_line_seg) <- c('vlo','vup')
  return(list(curr_line_seg,graph_equal,steps_taken))

}



