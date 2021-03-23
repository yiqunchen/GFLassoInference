# We compute a solution path of the fused lasso dual problem:
#
# \hat{u}(\lambda) =
# \argmin_u \|y - D^T u\|_2^2 \rm{s.t.} \|\u\|_\infty \leq \lambda
#
# where D is the incidence matrix of a given graph.
#
# Note: the df estimates at each lambda_k can be thought of as the df
# for all solutions corresponding to lambda in (lambda_k,lambda_{k-1}),
# the open interval to the *right* of the current lambda_k.

### modified to compute post selection events 
# input contrast vector v


#### indexed by number of CC

dualpathFused_CC_indexed <- function(y, D, v, sigma, K_CC = NULL, K_lambda = NULL,
                                     approx=FALSE, maxsteps=2000, minlam=0, maxsteps_factor = 100,
                             rtol=1e-7, btol=1e-7, cdtol = 1e-9,verbose=FALSE,
                             object=NULL, segment_list = NULL, stop_criteria = 'K') {

  Ds = NULL # for saving purposes
  
  if (!stop_criteria%in%c('K','CC','lambda')){
    stop('Criteria can only be K, CC, or lambda\n')
  }
  
  if((stop_criteria =='CC')&(is.null(K_CC))){
    stop('must specify K_CC when stopping with K\n')
  }
  
  if((stop_criteria =='lambda')&(is.null(K_lambda))){
    stop('must specify K_lambda when stopping with lambda\n')
  }
  
  if(stop_criteria%in%c('lambda','CC')){
    maxsteps <- maxsteps*maxsteps_factor
  }
  
  # If we are starting a new path
  #cat('segment_list',segment_list[[1]],'\n')
  if (is.null(object)) {
    m = nrow(D)
    n = ncol(D)
    
    # a list to keep track of all line segments
    LS_list <- list() 
    
    
    # Find the minimum 2-norm solution, using some linear algebra
    # tricks and a little bit of graph theory
    L = abs(crossprod(D))
    diag(L) = 0
    gr = graph.adjacency(L,mode="upper") # Underlying graph
    cl = clusters(gr)
    q = cl$no                            # Number of clusters
    i = cl$membership                    # Cluster membership
    x = f = numeric(n)  # solving linear system efficiently
    xv = numeric(n) # solving linear system with v efficiently
    
    # For efficiency, don't loop over singletons
    tab = tabulate(i)
    oo = which(tab[i]==1)
    if (length(oo)>0) {
      f[oo] = y[oo]
    }
    
    
    # Same for groups with two elements (doubletons?)
    oi = order(i)
    oo = which(tab[i][oi]==2)
    if (length(oo)>0) {
      mm = colMeans(matrix(y[oi][oo],nrow=2))
      f[oi][oo] = rep(mm,each=2)
      # added for [post selection]
      mmv = colMeans(matrix(v[oi][oo],nrow=2))
      ii = oo[Seq(1,length(oo),by=2)]
      x[oi][ii] = y[oi][ii] - mm # 
      # added for [post selection]
      xv[oi][ii] = v[oi][ii] - mmv
    }
    
    # Now all groups with at least three elements
    cs = cumsum(tab)
    grps = which(tab>2)
    for (j in grps) {
      oo = oi[Seq(cs[j]-tab[j]+1,cs[j])]
      yj = y[oo]
      # added for [post selection]
      vj = v[oo]
      f[oo] = mean(yj)
      # individual Laplacian
      Lj = crossprod(Matrix(D[,oo[-1]],sparse=TRUE))
      x[oo][-1] = as.numeric(solve(Lj,as.numeric(yj-mean(yj))[-1]))
      # # added for [post selection]
      xv[oo][-1] = as.numeric(solve(Lj,as.numeric(vj-mean(vj))[-1]))
    }
    
    uhat = as.numeric(D%*%x)     # Dual solution
    uvhat = as.numeric(D%*%xv)     # Dual solution with v replacing y
    betahat = f                  # Primal solution
    ihit = which.max(abs(uhat))  # Hitting coordinate
    hit = abs(uhat[ihit])        # Critical lambda
    s = Sign(uhat[ihit])         # Sign
    
    
    if (verbose) {
      cat(sprintf("1. lambda=%.3f, adding coordinate %i, |B|=%i...",
                  hit,ihit,1))
    }
    
    # Now iteratively find the new dual solution, and
    # the next critical point
    
    # Things to keep track of, and return at the end
    buf = min(maxsteps,1000)
    lams = numeric(buf)        # Critical lambdas
    h = logical(buf)           # Hit or leave?
    df = numeric(buf)          # Degrees of freedom
    action = numeric(buf)      # Action taken
    # cat('hit', hit,'\n')
    lams[1] = hit
    h[1] = TRUE
    df[1] = q
    action[1] = ihit # first action is hit
    
    u = matrix(0,m,buf)      # Dual solutions
    beta = matrix(0,n,buf)   # Primal solutions
    u[,1] = uhat
    beta[,1] = betahat
    
    # Update our graph
    e = which(D[ihit,]!=0)
    gr[e[1],e[2]] = 0             # Delete edge
    newcl = subcomponent(gr,e[1]) # New cluster
    oldcl = which(i==i[e[1]])     # Old cluster
    # If these two clusters aren't the same, update
    # the memberships
    if (length(newcl)!=length(oldcl) || any(sort(newcl)!=sort(oldcl))) {
      i[newcl] = q+1
      q = q+1
    }
    
    
    # rows to add, for first hitting time (no need for sign-eligibility--just add all sign pairs)
    
    # added for [post selection]
    Gy_init <- t(t(c(s*uhat[ihit]-uhat,s*uhat[ihit]+uhat)))
    Gv_init <- t(t(c(s*uvhat[ihit]-uvhat, s*uvhat[ihit]+uvhat)))
    
    LS_t <- PolyInt_G_free(y = y, v = v, Gy = Gy_init,
                           Gv = Gv_init, sigma, 
                           tol = 1e-6, bits=NULL)
    LS_list <- c(LS_list, LS_t) # add the truncated line segment
    
    # Other things to keep track of
    r = 1                      # Size of boundary set
    B = ihit                   # Boundary set
    I = Seq(1,m)[-ihit]        # Interior set
    D1 = D[-ihit,,drop=FALSE]  # Matrix D[I,]
    D2 = D[ihit,,drop=FALSE]   # Matrix D[B,]
    k = 2                      # What step are we at?




  }
  
  # If iterating an already started path
  else {
    # Grab variables needed to construct the path
    lambda = NULL
    for (j in 1:length(object)) {
      if (names(object)[j] != "pathobjs") {
        assign(names(object)[j], object[[j]])
      }
    }
    for (j in 1:length(object$pathobjs)) {
      assign(names(object$pathobjs)[j], object$pathobjs[[j]])
    }
    lams = lambda
  }
  
  tryCatch({
    while (k<=maxsteps && lams[k-1]>=minlam) {

               # check if we can just end early - early stopping by segment 
      if (stop_criteria =='CC'){
        if(length(unique(i))>=K_CC){
          #cat('step ',k ,'yields right CC\n')
          break
        }
      }
      
      if (stop_criteria =='lambda'){
        #cat('lambda', K_lambda,'\n')
        if(lams[1]<=K_lambda){
          #cat('step ',k ,'yields right lambda\n')
          break
        }
      }
      
      if  (stop_criteria =='K'){
      if (!is.null(segment_list)){
        #cat('segment_list',segment_list[[1]],'\n')
        early_stop <- check_segment_in_model(segment_list, i)
        if (early_stop){
          # cat(paste0('early stopping at step ',k,'\n'))
          break
        }
      }
      }
      
      ##########
      # Check if we've reached the end of the buffer
      if (k > length(lams)) {
        buf = length(lams)
        lams = c(lams,numeric(buf))
        action = c(action,numeric(buf)) # action
        h = c(h,logical(buf))
        df = c(df,numeric(buf))
        u = cbind(u,matrix(0,m,buf))
        beta = cbind(beta,matrix(0,n,buf))
      }
      
      ##########
      Ds = as.numeric(t(D2)%*%s)
      
      # If the interior is empty, then nothing will hit
      if (r==m) {
        fa = y
        fa_v = v
        fb = Ds
        hit = 0
      }
      
      # Otherwise, find the next hitting time
      else {
        xa = xb = numeric(n)
        fa = fb = numeric(n)
        xa_v = fa_v =  numeric(n) # final same a - except substitute y with v
        # For efficiency, don't loop over singletons
        tab = tabulate(i)
        oo = which(tab[i]==1)
        if (length(oo)>0) {
          fa[oo] = y[oo]
          fb[oo] = Ds[oo]
          fa_v[oo] = v[oo]
        }
        
        # Same for groups with two elements (doubletons?)
        oi = order(i)
        oo = which(tab[i][oi]==2)
        if (length(oo)>0) {
          ma = colMeans(matrix(y[oi][oo],nrow=2))
          mb = colMeans(matrix(Ds[oi][oo],nrow=2))
          # added for [post selection]
          ma_v = colMeans(matrix(v[oi][oo],nrow=2))
          fa[oi][oo] = rep(ma,each=2)
          fb[oi][oo] = rep(mb,each=2)
          fa_v[oi][oo] = rep(ma_v,each=2)
          ii = oo[Seq(1,length(oo),by=2)]
          xa[oi][ii] = y[oi][ii] - ma
          xb[oi][ii] = Ds[oi][ii] - mb
          # added for [post selection]
          xa_v[oi][ii] = v[oi][ii] - ma_v
        }
        
        # Now all groups with at least three elements
        cs = cumsum(tab)
        grps = which(tab>2)
        for (j in grps) {
          oo = oi[Seq(cs[j]-tab[j]+1,cs[j])]
          yj = y[oo]
          vj = v[oo]
          Dsj = Ds[oo]
          fa[oo] = mean(yj)
          fb[oo] = mean(Dsj)
          fa_v[oo] = mean(vj)
          Lj = crossprod(Matrix(D1[,oo[-1]],sparse=TRUE))
          xa[oo][-1] = as.numeric(solve(Lj,as.numeric(yj-mean(yj))[-1]))
          xb[oo][-1] = as.numeric(solve(Lj,as.numeric(Dsj-mean(Dsj))[-1]))
          # added for [post selection]
          xa_v[oo][-1] = as.numeric(solve(Lj,as.numeric(vj-mean(vj))[-1]))
        }
        
        a = as.numeric(D1%*%xa)
        # added for [post selection]
        a_v = as.numeric(D1%*%xa_v) 
        b = as.numeric(D1%*%xb)
        shits = Sign(a)
        hits = a/(b+shits)
        
        # Make sure none of the hitting times are larger
        # than the current lambda (precision issue)
        hits[hits>lams[k-1]+btol] = 0
        hits[hits>lams[k-1]] = lams[k-1]
        
        ihit = which.max(hits)
        hit = hits[ihit]
        shit = shits[ihit]
        
        
        a_over_b_ihit <- a[ihit]/(shit+b[ihit])
        a_v_over_b_ihit <- a_v[ihit]/(shit+b[ihit])
        ### added for [post selection]
        # adding gamma matrices for the hitting time
        Gy_hit <- t(t(c(shits*a, (a[ihit]/(shit+b[ihit]))-a/(shits+b) )))
        Gv_hit <- t(t(c(shits*a_v,  (a_v[ihit]/(shit+b[ihit]))-a_v/(shits+b) )))
        LS_t <- PolyInt_G_free(y = y, v = v, Gy = Gy_hit,
                               Gv = Gv_hit, sigma, 
                               tol = 1e-6, bits=NULL)
        LS_list <- c(LS_list, LS_t) # add the truncated line segment
      }
      
      ##########
      # If nothing is on the boundary, then nothing will leave
      # Also, skip this if we are in "approx" mode
      if (r==0 || approx) {
        leave = 0
      }
      
      # Otherwise, find the next leaving time
      else {
        c = as.numeric(s*(D2%*%fa))
        c_v = as.numeric(s*(D2%*%fa_v)) # substitute y with v in c's definition
        c[abs(c)<cdtol] = 0 # round small error
        d = as.numeric(s*(D2%*%fb))
        leaves = c/d
        
        # c must be negative
        leaves[c>=0|d>=0] = 0
        
        # Make sure none of the leaving times are larger
        # than the current lambda (precision issue)
        leaves[leaves>lams[k-1]+btol] = 0
        leaves[leaves>lams[k-1]] = lams[k-1]
        
        
        ileave = which.max(leaves)
        leave = leaves[ileave]
        
        # identify which on boundary set are c<0 and d<0
        C_i = (c < 0)
        D_i = (d < 0)
        C_i[!B] = D_i[!B] = FALSE
        
        CDi = (C_i & D_i)
        CDi[ileave] = FALSE
        CDind = which(CDi) 
        
        
        ### added for [post selection]
        # adding gamma matrices for the leaving time
        c_vec <- c # disambigous names...
        d_vec <- d
        
        if(d_vec[ileave]!=0){
          c_over_d_ileave <- leave#c_vec[ileave]/d_vec[ileave]
          c_v_over_d_ileave <- c_v[ileave]/d_vec[ileave]
        }else{
          c_over_d_ileave <- NULL
          c_v_over_d_ileave <- NULL
        }
        
        c_over_d <- c_vec[CDind]/d_vec[CDind]
        
        
        c_v_over_d <- c_v[CDind]/d_vec[CDind]
        
        Gy_leave_1 <- t(t(c( (-1)*c_vec[C_i&D_i], c_vec[(!C_i)&D_i] )))
        Gv_leave_1 <- t(t(c( (-1)*c_v[C_i&D_i], c_v[(!C_i)&D_i]  )))
        if (length(Gy_leave_1)>0&(d_vec[ileave]!=0)){
          LS_t <- PolyInt_G_free(y = y, v = v, Gy = Gy_leave_1,
                                 Gv = Gv_leave_1, sigma, 
                                 tol = 1e-6, bits=NULL)
          LS_list <- c(LS_list, LS_t) # add the truncated line segment
          
          Gy_leave_2 <- t(t(c(c_over_d_ileave-c_over_d)))
          Gv_leave_2 <- t(t(c(c_v_over_d_ileave-c_v_over_d)))
          
          LS_t <- PolyInt_G_free(y = y, v = v, Gy = Gy_leave_2,
                                 Gv = Gv_leave_2, sigma, 
                                 tol = 1e-6, bits=NULL)
          LS_list <- c(LS_list, LS_t) # add the truncated line segment
        }
        # add rows to get the leaving time characterization
        
        
        
      }
      
      ##########
      # Stop if the next critical point is negative
      if (hit<=0 && leave<=0) break
      
      # If a hitting time comes next
      if (hit > leave) {
        # Record the critical lambda and properties
        lams[k] = hit
        action[k] = I[ihit]
        h[k] = TRUE
        df[k] = q
        uhat = numeric(m)
        uhat[B] = hit*s
        uhat[I] = a-hit*b
        betahat = fa-hit*fb
        
        # Update our graph
        e = which(D1[ihit,]!=0)
        gr[e[1],e[2]] = 0             # Delete edge
        newcl = subcomponent(gr,e[1]) # New cluster
        oldcl = which(i==i[e[1]])     # Old cluster
        # If these two clusters aren't the same, update
        # the memberships
        if (length(newcl)!=length(oldcl) || any(sort(newcl)!=sort(oldcl))) {
          i[newcl] = q+1
          q = q+1
        }
        
        if((!is.null(c_over_d_ileave))&(leave!=0)){
          ##### add to gamma!; this is if hit comes next 
          Gy_hit_geq_leave <- t(t(1*c(a_over_b_ihit-c_over_d_ileave)))
          Gv_hit_geq_leave <- t(t(1*c(a_v_over_b_ihit-c_v_over_d_ileave)))
          #cat(paste0('a_over_b_ihit',a_over_b_ihit,'hit',hit,'c_over_d_ileave',
          #c_over_d_ileave,'leave',leave,'hit > leave, 1 ,\n'))
          LS_t <- PolyInt_G_free(y = y, v = v, Gy = Gy_hit_geq_leave,
                                 Gv = Gv_hit_geq_leave, sigma, 
                                 tol = 1e-6, bits=NULL)
          #cat('hit > leave, 2 ,\n')
          LS_list <- c(LS_list, LS_t) # add the truncated line segment
        }else{
          # else nothing
        }
        
        
        # Update all other variables
        r = r+1
        B = c(B,I[ihit])
        I = I[-ihit]
        s = c(s,shit)
        D2 = rbind(D2,D1[ihit,])
        D1 = D1[-ihit,,drop=FALSE]
        
        if (verbose) {
          cat(sprintf("\n%i. lambda=%.3f, adding coordinate %i, |B|=%i...",
                      k,hit,B[r],r))
        }
      }
      
      # Otherwise a leaving time comes next
      else {
        # Record the critical lambda and properties
        lams[k] = leave
        action[k] = -B[ileave]
        h[k] = FALSE
        df[k] = q
        uhat = numeric(m)
        uhat[B] = leave*s
        uhat[I] = a-leave*b
        betahat = fa-leave*fb
        
        if(d_vec[ileave]!=0&(!is.null(c_over_d_ileave))){
          ##### add to gamma!; this is if hit comes next 
          Gy_leave_geq_hit <- t(t(-1*c(a_over_b_ihit-c_over_d_ileave)))
          Gv_leave_geq_hit <- t(t(-1*c(a_v_over_b_ihit-c_v_over_d_ileave)))
          
          #cat(paste0('a_over_b_ihit',a_over_b_ihit,'hit',hit,'c_over_d_ileave',
          #  c_over_d_ileave,'leave',leave,'hit < leave, 1 ,\n'))
          
          LS_t <- PolyInt_G_free(y = y, v = v, Gy = Gy_leave_geq_hit,
                                 Gv = Gv_leave_geq_hit, sigma, 
                                 tol = 1e-6, bits=NULL)
          #cat('hit < leave, 2 ,\n')
          LS_list <- c(LS_list, LS_t) # add the truncated line segment
        }else{
          
        }
        
        
        # Update our graph
        e = which(D2[ileave,]!=0)
        gr[e[1],e[2]] = 1             # Add edge
        newcl = subcomponent(gr,e[1]) # New cluster
        oldcl = which(i==i[e[1]])     # Old cluster
        # If these two clusters aren't the same, update
        # the memberships
        if (length(newcl)!=length(oldcl) || !all(sort(newcl)==sort(oldcl))) {
          newno = i[e[2]]
          oldno = i[e[1]]
          i[oldcl] = newno
          i[i>oldno] = i[i>oldno]-1
          q = q-1
        }
        

        # Update all other variables
        r = r-1
        I = c(I,B[ileave])
        B = B[-ileave]
        s = s[-ileave]
        D1 = rbind(D1,D2[ileave,])
        D2 = D2[-ileave,,drop=FALSE]
        
        
        if (verbose) {
          cat(sprintf("\n%i. lambda=%.3f, deleting coordinate %i, |B|=%i...",
                      k,leave,I[m-r],r))
        }
      }
      
      
      
      # check if we can just end early - early stopping by segment 
      if (stop_criteria =='CC'){
        if(length(unique(i))>=K_CC){
          #cat('step ',k ,'yields right CC\n')
          action = action[action!=0]
          break
        }
      }
      
      if (stop_criteria =='lambda'){
        #cat('lambda', K_lambda,'\n')
        if(lams[k]<=K_lambda){
          #cat('step ',k ,'yields right lambda\n')
          action = action[action!=0]
          break
        }
      }
      
      if (stop_criteria =='K'){
      if (!is.null(segment_list)){
        #cat('segment_list',segment_list[[1]],'\n')
        early_stop <- check_segment_in_model(segment_list, i)
        if (early_stop){
          # cat(paste0('early stopping at step ',k,'\n'))
          action = action[action!=0]
          break
        }
      }
      }
      
      u[,k] = uhat
      beta[,k] = betahat
      
      # Step counter
      k = k+1
    }
  }, error = function(err) {
    err$message = paste(err$message,"\n(Path computation has been terminated;",
                        " partial path is being returned.)",sep="")
    warning(err)})
  
  # Trim
  lams = lams[Seq(1,k-1)]
  h = h[Seq(1,k-1)]
  df = df[Seq(1,k-1)]
  u = u[,Seq(1,k-1),drop=FALSE]
  beta = beta[,Seq(1,k-1),drop=FALSE]
  
  # If we reached the maximum number of steps
  if (k>maxsteps) {
    if (verbose) {
      cat(sprintf("\nReached the maximum number of steps (%i),",maxsteps))
      cat(" skipping the rest of the path.")
    }
    completepath = FALSE
  }
  
  # If we reached the minimum lambda
  else if (lams[k-1]<minlam) {
    if (verbose) {
      cat(sprintf("\nReached the minimum lambda (%.3f),",minlam))
      cat(" skipping the rest of the path.")
    }
    completepath = FALSE
  }
  
  # Otherwise, note that we completed the path
  else completepath = TRUE
  
  # The least squares solution (lambda=0)
  bls = y
  if (verbose) cat("\n")
  
  # Save needed elements for continuing the path
  pathobjs = list(type="fused",r=r, B=B, I=I, Q1=NA, approx=approx,
                  Q2=NA, k=k, df=df, D1=D1, D2=D2, Ds=Ds, ihit=ihit, m=m, n=n, q=q,
                  h=h, q0=NA, rtol=rtol, btol=btol, s=s, y=y, gr=gr, i=i, membership = i)
  
  colnames(u) = as.character(round(lams,3))
  colnames(beta) = as.character(round(lams,3))
  
  states = get.states(action)
  action = action[action!=0]
  
  return(list(lambda=lams,beta=beta,u=u,hit=h,df=df,y=y,states=states,
              action = action,
              completepath=completepath,bls=bls,pathobjs=pathobjs,
              LS_list=LS_list))
}
