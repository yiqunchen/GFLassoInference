dualpathSvd2_YC <- function(y, D, v, sigma = 1, approx=FALSE, maxsteps=2000, minlam=0,
                            rtol=1e-7, btol=1e-7, verbose=FALSE, object=NULL,
                            ctol=1e-10, cdtol=1e-4, do.gc=F){

  # Error checking
  stopifnot(ncol(D) == length(y))

  nk = 0
  ss = list() # list of ss
  LS_list <- list() # keep all the truncated line seg

  # If we are starting a new path
  if (is.null(object)) {
    m = nrow(D)
    n = ncol(D)

    # Initialize Gamma matrix (work in progress)
    # Gammat = Matrix(NA,nrow= maxsteps*ncol(D)*4 ,ncol=n)

    # Compute the dual solution at infinity, and
    # find the first critical point
    D_sparse_T <- Matrix(t(D),sparse=T)
    In = diag(1,n)
    sv = sparse_svdsolve(D_sparse_T,y,rtol)
    uhat = as.numeric(sv$x[,1])        # Dual solution
    q = sv$q                           # Rank of D

    ihit = which.max(abs(uhat))   # Hitting coordinate
    hit = abs(uhat[ihit])         # Critical lambda
    s = Sign(uhat[ihit])          # Sign
    k = 1
    ss[[k]] = s

    if (verbose) {
      cat(sprintf("1. lambda=%.3f, adding coordinate %i, |B|=%i...",
                  hit,ihit,1))
    }

    # Now iteratively find the new dual solution, and
    # the next critical point

    # Things to keep track of, and return at the end
    buf = min(maxsteps,1000)
    u = matrix(0,m,buf)        # Dual solutions
    lams = numeric(buf)        # Critical lambdas
    h = logical(buf)           # Hit or leave?
    df = numeric(buf)          # Degrees of freedom
    action = numeric(buf)      # Action taken
    upol = c()                 # Constant in polyhedral constraint

    lams[1] = hit
    action[1] = ihit
    h[1] = TRUE
    df[1] = n-q
    u[,1] = uhat


    # add rows to Gamma
    tDinv <- sv$gen_inv  #ginv(as.matrix(t(D)))

    # rows to add, for first hitting time (no need for sign-eligibility--just add all sign pairs)
    #Gammat = diag(Sign(uhat)  %*% tDinv)
    M = matrix(s*tDinv[ihit,], nrow(tDinv[-ihit,]), n, byrow=TRUE)
    Gammat = rbind(M + tDinv[-ihit,],
                   M - tDinv[-ihit,])


    #nk = nrow(Gammat)
    # first gamma !
    LS_t <- PolyInt(y = y, G = Gammat, v = v, sigma = sigma)
    LS_list <- c(LS_list, LS_t) # add the truncated line segment

    # Other things to keep track of, but not return
    r = 1                      # Size of boundary set
    B = ihit                   # Boundary set
    I = Seq(1,m)[-ihit]        # Interior set
    Ds = D[ihit,]*s            # Vector t(D[B,])%*%s
    D1 = D[-ihit,,drop=FALSE]  # Matrix D[I,]
    D2 = D[ihit,,drop=FALSE]   # Matrix D[B,]
    k = 2                      # What step are we at?

    D1_sparse_T <- Matrix(t(D1),sparse=T)

  } else {
    # If iterating an already started path
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

      if(do.gc) gc()
      if(verbose){ cat('\n'); show.glutton(environment(),4)}

      ##########
      # Check if we've reached the end of the buffer
      if (k > length(lams)) {
        buf = length(lams)
        lams = c(lams,numeric(buf))
        action = c(action,numeric(buf))
        h = c(h,logical(buf))
        df = c(df,numeric(buf))
        u = cbind(u,matrix(0,m,buf))
      }

      ##########
      # If the interior is empty, then nothing will hit
      if (r==m) {
        a = b = numeric(0)
        hit = 0
        q = 0
      }
      # Otherwise, find the next hitting time
      else {
        In = diag(1,n)
        sv = sparse_svdsolve(D1_sparse_T,cbind(y,Ds),rtol)
        #sv = svdsolve(t(D1),cbind(y,Ds,In),rtol)
        a = as.numeric(sv$x[,1])  # formerly a = as.numeric(D3 %*% y)
        b = as.numeric(sv$x[,2])
        D3 =  sv$gen_inv#ginv(as.matrix(t(D1)))# formerly as.matrix(sv$x[,3:(n+2)])
        # meant to be pseudoinverse of t(D[-I,])

        q = sv$q
        shits = Sign(a)
        hits = a/(b+shits);


        # Make sure none of the hitting times are larger
        # than the current lambda (precision issue)
        hits[hits>lams[k-1]+btol] = 0
        hits[hits>lams[k-1]] = lams[k-1]

        ihit = which.max(hits)
        hit = hits[ihit]
        shit = shits[ihit]

        # Gamma Matrix!
        # rows to add, for viable hitting signs:
        tDinv = D3
        rows.to.add = (if(length(shits)>1){
          do.call(rbind, lapply(1:length(shits), function(ii){shits[ii] * tDinv[ii,]  }))
        } else {
          shits* tDinv
        })
        # add line segments

        #Gammat = rbind(Gammat, rows.to.add)
        LS_t <- PolyInt(y = y, G = rows.to.add, v = v, sigma = sigma)
        LS_list <- c(LS_list, LS_t) # add the truncated line segment

        # rows to add, for hitting event: (This is just dividing each row of tDinv by corresponding element of tDinv%*%Ds+shit)
        A = asrowmat(D3/(b+shits)) # tDinv / as.numeric(tDinv %*% Ds + shits)
        if(nrow(A)!=1){
          nleft = nrow(A[-ihit,])
          if(is.null(nleft)) nleft = 1
          M = matrix(A[ihit,], nrow = nleft, ncol = n, byrow = TRUE)
          #Gammat = rbind(Gammat, M - A[-ihit,])
          LS_t <- PolyInt(y = y, G =  M - A[-ihit,], v = v, sigma = sigma)
          LS_list <- c(LS_list, LS_t) # add the truncated line segment
        }
      }

      ##########
      # If nothing is on the boundary, then nothing will leave
      # Also, skip this if we are in "approx" mode
      if (r==0 || approx) {
        leave = 0
      }

      # Otherwise, find the next leaving time
      else {
        c = as.matrix(s*(D2%*%(y-D1_sparse_T%*%a)))
        d = as.matrix(s*(D2%*%(Ds-D1_sparse_T%*%b)))


        # round small values of c to zero (but not d)
        #cdtol = 1E-10
        c[abs(c) < cdtol] = 0

        # get leave times
        leaves = c/d

        # identify which on boundary set are c<0 and d<0
        Ci = (c < 0)
        Di = (d < 0)
        Ci[!B] = Di[!B] = FALSE
        CDi = (Ci & Di)

        # c and d must be negative at all coordinates to be considered
        leaves[c>=0|d>=0] = 0

        # Make sure none of the leaving times are larger
        # than the current lambda (precision issue)
        super.lambda = leaves>lams[k-1]+btol
        leaves[leaves>lams[k-1]+btol] = 0
        leaves[leaves>lams[k-1]] = lams[k-1]

        # If a variable just entered, then make sure it
        # cannot leave (added from lasso.R)
        if (action[k-1]>0) leaves[r] = 0

        # index of leaving coordinate
        ileave = which.max(leaves)
        leave = leaves[ileave]

        # Gamma Matrix!!
        # rows to add, for leaving event:
        if(dim(D1)[1]==0) {
          D1 = rep(0,ncol(D1))
          D1_sparse_T <- Matrix(t(D1),sparse=T)
        } # temporarily added because of
        # dimension problem in next line,
        # at last step of algorithm
        gmat = s*(D2%*%(In - D1_sparse_T%*%D3)) # coefficient matrix to c

        # close-to-zero replacement is hard-coded in
        gmat[abs(c)<cdtol,] = rep(0,ncol(gmat))

        # we still want to see that gmat&%y ~= c
        if(!(max(gmat%*%y-c) < ctol)) print(max(gmat%*%y-c))

        gd = gmat / as.numeric(d)


        # hard coding in the zero replacements, to make it identical with lea
        gd[c>=0,] = rep(0, ncol(gd))        # re-doing zero-replacement in the g/d matrix
        gd[super.lambda,] = rep(0,ncol(gd)) # re-doing larger-than-lambda-replacement
        #gd[closeto.lambda,] = rep(0,ncol(gd)) # re-doing close-to-lambda-replacement (not sure how to make the i'th leaving time gd[i,]%*%y == lam[k-1] properly;
        # solving to get some gd[i,] is like finding a non-unique solution to an overdetermined system; because such a gd[i,] is not unique,
        # how do I know that adding this row to Gamma won't do weird and mysterious things?)

        #if( (length(Di)!=0) & (which(closeto.lambda) %in% which(Di))) print("closeto.lambda replacement SHOULD have happenned (but didn't).")

        # add rows that ensure c<0 #(only in )
        #Gammat <- rbind(Gammat, gmat[Ci&Di,]*(-1),
        #               gmat[(!Ci)&Di,])

        LS_t <- PolyInt(y = y, G = gmat[Ci&Di,]*(-1), v = v, sigma = sigma)
        LS_list <- c(LS_list, LS_t) # add the truncated line segment
        LS_t <- PolyInt(y = y, G = gmat[(!Ci)&Di,], v = v, sigma = sigma)
        LS_list <- c(LS_list, LS_t) # add the truncated line segment
        #        print("after leave eligibility (c<0)")
        #        print(nrow(Gammat))

        # get rid of NA rows in Gammat (temporary fix)
        #Gammat <- na.omit(Gammat) # faster; r naive fix
        #missing.rows = apply(Gammat, 1, function(row) any(is.na(row)))
        #if(sum(missing.rows)>=1){ Gammat <- Gammat[-which(missing.rows),] }

        # add rows for maximizer
        CDi = (Ci & Di)
        CDi[ileave] = FALSE
        CDind = which(CDi)

        #Gammat <- rbind(Gammat, gd[rep(ileave,length(CDind)),] - gd[CDind,])
        #        print("after leave times")
        #        print(nrow(Gammat))
        LS_t <- PolyInt(y = y, G = gd[rep(ileave,length(CDind)),] - gd[CDind,], v = v, sigma = sigma)
        LS_list <- c(LS_list, LS_t) # add the truncated line segment
      }
      ##########
      # Stop if the next critical point is negative

      if (hit<=0 && leave<=0) {break}

      # If a hitting time comes next
      if (hit > leave) {

        # Record the critical lambda and solution
        lams[k] = hit
        action[k] = I[ihit]
        h[k] = TRUE
        df[k] = n-q
        uhat = numeric(m)
        uhat[B] = hit*s
        uhat[I] = a-hit*b
        u[,k] = uhat


        # add row to Gamma to characterize the hit coming next
        if(!approx)  {
          # this is literally h_k - l_k > 0
          #Gammat = rbind(Gammat,  A[ihit,] - gd[ileave,])
          LS_t <- PolyInt(y = y, G =  A[ihit,] - gd[ileave,], v = v, sigma = sigma)
          LS_list <- c(LS_list, LS_t) # add the truncated line segment
          #          print("after hit vs leave (hit wins)")
          #          print(nrow(Gammat))
        }

        #nk = c(nk,nrow(Gammat))

        # Update all of the variables
        r = r+1
        B = c(B,I[ihit])
        I = I[-ihit]
        Ds = Ds + D1[ihit,]*shit
        s = c(s,shit)
        D2 = rbind(D2,D1[ihit,])
        D1 = D1[-ihit,,drop=FALSE]
        D1_sparse_T <- Matrix(t(D1),sparse=T)
        ss[[k]] = s
        if (verbose) {
          cat(sprintf("\n%i. lambda=%.3f, adding coordinate %i, |B|=%i...",
                      k,hit,B[r],r))
        }
      }

      # Otherwise a leaving time comes next
      else {
        # Record the critical lambda and solution
        lams[k] = leave
        action[k] = -B[ileave]
        h[k] = FALSE
        df[k] = n-q
        uhat = numeric(m)
        uhat[B] = leave*s
        uhat[I] = a-leave*b
        u[,k] = uhat


        # add row to Gamma to characterize the leave coming next
        if(!approx)  {
          #Gammat = rbind(Gammat, - A[ihit,] + gd[ileave,])
          LS_t <- PolyInt(y = y, G =  - A[ihit,] + gd[ileave,], v = v, sigma = sigma)
          LS_list <- c(LS_list, LS_t)

        }
        # nk = c(nk,nrow(Gammat))

        # Update all of the variables
        r = r-1
        I = c(I,B[ileave])
        B = B[-ileave]
        Ds = Ds - D2[ileave,]*s[ileave]
        s = s[-ileave]
        D1 = rbind(D1,D2[ileave,])
        D1_sparse_T <- Matrix(t(D1),sparse=T)
        D2 = D2[-ileave,,drop=FALSE]
        ss[[k]] = s

        if (verbose) {
          cat(sprintf("\n%i. lambda=%.3f, deleting coordinate %i, |B|=%i...",
                      k,leave,I[m-r],r))
        }

      }

      # Step counter + other stuff
      k = k+1
      # resetting tDinv
      #tDinv = t(as.matrix(rep(NA,length(tDinv)))[,-1,drop=F])
    }
  }
  , error = function(err) {
    err$message = paste(err$message,"\n(Path computation has been terminated;",
                        " partial path is being returned.)",sep="")
    warning(err)})

  # Trim
  lams = lams[Seq(1,k-1)]
  h = h[Seq(1,k-1)]
  df = df[Seq(1,k-1),drop=FALSE]
  u = u[,Seq(1,k-1),drop=FALSE]

  # Save needed elements for continuing the path
  pathobjs = list(type="svd", r=r, B=B, I=I, approx=approx,
                  k=k, df=df, D1=D1, D2=D2, Ds=Ds, ihit=ihit, m=m, n=n, h=h,
                  rtol=rtol, btol=btol, s=s, y=y)
  # If we reached the maximum number of steps
  if (k>maxsteps) {
    if (verbose) {
      cat(sprintf("\nReached the maximum number of steps (%i),",maxsteps))
      cat(" skipping the rest of the path.")
    }
    completepath = FALSE
  } else if (lams[k-1]<minlam) {
    # If we reached the minimum lambda
    if (verbose) {
      cat(sprintf("\nReached the minimum lambda (%.3f),",minlam))
      cat(" skipping the rest of the path.")
    }
    completepath = FALSE
  } else {
    # Otherwise, note that we completed the path
    completepath = TRUE
  }
  if (verbose) cat("\n")

  colnames(u) = as.character(round(lams,3))

  beta <-  apply(t(D)%*%u,2,function(column){y-column})
  ss = c(NA,ss)
  states = get.states(action)

  return(list(lambda=lams,beta=beta,u=u,hit=h,df=df,y=y,ss=ss,states=states,
              completepath=completepath,bls=y,pathobjs=pathobjs,
              action=action, LS_list=LS_list))
}




