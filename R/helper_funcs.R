
# Computes x = A^+ * b using an SVD
# (slow but stable)
svdsolve <- function(A,b,rtol) {
  s = svd(A)
  di = s$d
  ii = di>rtol
  di[ii] = 1/di[ii]
  di[!ii] = 0
  return(list(x=s$v%*%(di*(t(s$u)%*%b)),q=sum(ii)))
}

# Computes x = A^+ * b using a sparse SVD
# (a tiny bit faster than svdsolve)s
sparse_svdsolve <- function(A,b, rtol=1e-7) {
  s <- sparsesvd(A, rank=0L, tol=1e-7, kappa=1e-6)
  di <- s$d
  ii <- (di>rtol)
  di[!ii] <- 0
  gen_inv <- s$v %*% (di * t(s$u))
  return(list(x=s$v%*%(di*(t(s$u)%*%b)),q=sum(ii),s=s,
              gen_inv=gen_inv))
}


#' Utility function for creating the penalty matrix D that corresponds to a chain graph of length n
#' @keywords internal
#' @param n: length of the chain graph
#' @return An (n-1)x(n) penalty matrix D
#' @examples
#' n <- 10
#' D_1d <- dual1d_Dmat(n)
dual1d_Dmat = function(m){
  D = matrix(0, nrow = m-1, ncol = m)
  for(ii in 1:(m-1)){
    D[ii,ii] = -1
    D[ii,ii+1] = 1
  }
  return(D)
}


# Makes a D matrix for a matrix that is stacked as rows
#' @keywords internal
graph2D_Dmat = function(m){

  m0 = sqrt(m)
  stopifnot( m0 == round(m0) )
  D = matrix(0, nrow = 2*(m0-1)*m0, ncol = m)

  # connect all within chunks
  for(jj in 1:m0){
    chunk = 1:m0 + (jj-1)*m0
    D[(m0-1)*(jj-1) + 1:(m0-1), chunk] = dual1d_Dmat(m0)
  }
  sofar = m0*(m0-1)
  # connect between all adjacent chunks
  for(jj in 1:(m-m0)){
    newrow = rep(0,m)
    newrow[jj]    = -1
    newrow[jj+m0] =  1
    D[sofar+jj,] = newrow
  }
  return(D)
}

# takes in adjacency matrix (either produced from igraph,
#                            or a matrix with zero entries
#                            in non-adjacent edges and 1 in adjacent edges)
# and produces a (sparse) D matrix for usage in the generalized lasso
##' @import Matrix
##' @export
getDmat.from.adjmat = function(adjmat, sparseMatrix=FALSE){
  n = ncol(adjmat)
  if(sparseMatrix){
    Dmat = Matrix(0,nrow=n^2,ncol=n, sparse=TRUE)
  } else {
    Dmat = matrix(0,nrow=n^2,ncol=n)
  }
  count = 1
  for(jj in 1:n){
    for(kk in jj:n){
      if(adjmat[jj,kk]==1){ # there was an error !=1
        Dmat[count,jj] = 1
        Dmat[count,kk] = -1
        count = count+1
      }
    }
  }
  Dmat = Dmat[1:(count-1),]
  return(Dmat)
}


convert_segment <- function(membership_vec){
  unique_member <- unique(membership_vec)
  segment_list <- vector('list',length = length(unique_member))
  for (i in seq_along(unique_member)){
    segment_list[[i]] <- which(membership_vec==i)
  }
  return(segment_list)
}


convert_cc_to_list <- function(cc_vec){
  unique_cc <- sort(unique(cc_vec))
  cc_list <- vector('list', length = length(unique_cc))
  for (i in seq_along(unique_cc)){
    cc_list[[i]] <- which(cc_vec==unique_cc[i])
  }
  return(cc_list)
}

# Function to take in an action object  ex) f0 = dualpathSvd2(..); action.obj = f0$action;
# and returns a _list_ of the states after each step in algorithm, starting with NA as the first state.
get.states = function(action.obj){

  # Obtain a list of actions at each step.
  actionlists = sapply(1:length(action.obj),function(ii) action.obj[1:ii] )

  # Helper function to extract the final state after running through (a list of) actions
  get.final.state = function(actionlist){
    if(length(actionlist)==1) return(actionlist)
    to.delete = c()
    for(abs.coord in unique(abs(actionlist))){
      all.coord.actions = which(abs(actionlist)==abs.coord)
      if( length(all.coord.actions) > 1 ){
        if( length(all.coord.actions) %%2 ==1){
          to.delete = c(to.delete, all.coord.actions[1:(length(all.coord.actions)-1)])
        } else {
          to.delete = c(to.delete, all.coord.actions)
        }
      }
    }
    if(!is.null(to.delete)){
      return(actionlist[-to.delete])
    } else {
      return(actionlist)
    }
  }
  states = lapply(actionlists, get.final.state)
  states = c(NA,states)

  return(states)
}

# Returns the sign of x, with Sign(0)=1.
Sign <- function(x) {
  return(-1+2*(x>=0))
}


# Returns a sequence of integers from a to b if a <= b,
# otherwise nothing. You have no idea how important this
# function is...

Seq = function(a,b,by) {
  if (a<=b) return(seq(a,b,by=by))
  else return(integer(0))
}

#' check_segment_in_model; util function
#' @keywords internal
#' @details
#' This model takes in
#' 1. new_cluster new cluster assignment
#' 2. old_segment a list segment (represented by the positions of the nodes),
#' and checks if the list of segments of interest is in the new graph
#' (assume same ordering of nodes)
check_segment_in_model <- function(old_segment, new_cluster){

  old_segment_length <- length(old_segment)

  for (i in c(1:old_segment_length)){
    curr_seg <- old_segment[[i]]
    curr_elem <- curr_seg[1]
    new_name <- new_cluster[curr_elem]
    new_seg_nodes <- sort(which(new_cluster == new_name), decreasing = FALSE)
    sorted_seg_nodes <- sort(curr_seg,decreasing = FALSE)
    if (length(sorted_seg_nodes)!=length(new_seg_nodes)){
      return(FALSE)
    }else{
      if (!all(new_seg_nodes==sorted_seg_nodes)){
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

getadjmat.from.Dmat = function(Dmat){

  Dmat = rbind(Dmat)
  n = ncol(Dmat)
  adjmat = matrix(0,n,n)

  for(jj in 1:nrow(Dmat)){
    inds = which(Dmat[jj,]!=0)
    adjmat[inds[1],inds[2]] = 1
    adjmat[inds[2],inds[1]] = 1
  }
  return(adjmat)
}

GetGamma_k = function(obj, condition.step){
  ### get k_th (condition.step) Gammat matrix for conditioning polyhedron
  ### obj: output from dualpathSvd2()
  ### condition.step: k
  if(length(obj$action) < condition.step){
    stop("\n You must ask for polyhedron from less than ", length(obj$action), " steps! \n")
  }
  myGammat <- obj$Gamma[1:obj$nk[condition.step],]
  return(myGammat)
}

GetBeta_k = function(obj, condition.step){
  ### get the estimated beta(mean in the denoising case)
  ### at step k
  ### obj: output from dualpathSvd2()
  ### condition.step: k
  if(length(obj$action) < condition.step){
    stop("\n You must ask for polyhedron from less than ", length(obj$action), " steps! \n")
  }
  myBeta <- obj$beta[,condition.step]
  return(myBeta)
}

GetGraph_k = function(obj, condition.step, Dmat){
  ### get the estimated graph at step k
  ### obj: output from dualpathSvd2()
  ### condition.step: k
  ### Dmat: the graph incidence (penalty) matrix
  if(length(obj$action) < (condition.step-1)){
    stop("\n You must ask for polyhedron from less than ", length(obj$action), " steps! \n")
  }
  myStates <- get.states(obj$action)[[condition.step+1]]
  new_adjancey <- Dmat[-myStates,]
  myGraph <- graph_from_adjacency_matrix(getadjmat.from.Dmat(new_adjancey), mode="undirected")
  return(myGraph)
}

test.cc.equality <- function(g1, g2){
  ### compares if two connected components
  ### are the same - assume common indices
  ### Input: igraph object g1 and g2
  ### Output: logical vector
  cc_g1 <- components(g1)
  cc_g2 <- components(g2)
  # number of cc has to be equal
  if(cc_g1$no!=cc_g1$no){
    return(FALSE)
  }
  # since membership ordering is determined
  # by the vertex ordering; sufficient to check
  # the equality in membership vectors
  if(all(cc_g1$membership==cc_g2$membership)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}


# Check correctness of polyhedron:
polyhedron.checks.out = function(y0, G, u=ifelse(is.null(nrow(G)), 0,rep(0,nrow(G))),
                                 tol = 1e-7, throw.error=F){
  if(is.null(nrow(G))){
    return(NULL)
  }
  all.nonneg = all(G%*%y0-u>=-1*tol)
  if(all.nonneg){
    return(all.nonneg)
  } else {
    if(throw.error){
      print(paste0("min margin",min(G%*%y0),"\n"))
      print(paste0("max margin", max(G%*%y0),"\n"))
      stop("failed Gy >= u test")
    } else {
      print(paste0("min margin",min(G%*%y0),"\n"))
      print(paste0("max margin", max(G%*%y0),"\n"))
      print("failed Gy >= u test")
    }
    return(all.nonneg)
  }
}

PolyInt <- function(y, G, v, sigma, u=ifelse(is.null(nrow(G)), 0,rep(0,nrow(G))),
                    tol = 1e-6, bits=NULL){
  ### Compute the interval [vlo, vup] defined by y:Gy<=0
  ### [vlow, vup] = v^ty|Gy<=u
  ### Input: y (observed data), G (Gamma matrix), u (0 vector)
  ### v (contrast vector), sigma (variance; assumed to be known)
  ### bits (optional argument for mpfr; precise calculation)
  ### Output: an ordered tuple [vlow, vup]
  z = sum(v*y)
  # 2 norm square
  vv = sum(v^2)
  sd = sigma*sqrt(vv)
  if(is.null(nrow(G))){
    return(list(vlo=-Inf,vup=+Inf,sd=sd,z=z))
  }
  # rho as defined in hyun et al.
  rho = G %*% v / vv
  rho[abs(rho)<=tol] <- 0 # floating point round off
  vec = (u - G %*% y + rho*z) / rho
  vlo = suppressWarnings(max(vec[rho>0])) #remove suppress warnings
  vup = suppressWarnings(min(vec[rho<0])) #remove suppress warnings
  if(vlo>vup){
    stop('Error: vlow should be smaller than vup, invalid interval!')
  }
  return(list(vlo=vlo,vup=vup,sd=sd,z=z))
}

PolyInt_G_free <- function(y, v, Gy, Gv, sigma,
                           tol = 1e-6, gy_tol = 1e-9, bits=NULL){
  ### compute the same values as polyint but instead pre-input Gy and Gv instead!
  ### if we can computr Gy and Gv efficiently

  ### Compute the interval [vlo, vup] defined by y:Gy<=0
  ### [vlow, vup] = v^ty|Gy<=u
  ### Input: y (observed data), G (Gamma matrix), u (0 vector)
  ### v (contrast vector), sigma (variance; assumed to be known)
  ### bits (optional argument for mpfr; precise calculation)
  ### Output: an ordered tuple [vlow, vup]
  Gy[abs(Gy)<=gy_tol] <- 0 # machine precision

  if(!all(Gy>=0)){
    #
    cat('bad gy',Gy,'\n')
  }
  u=ifelse(is.null(Gy), 0,rep(0,length(Gy)))

  z = sum(v*y)
  # 2 norm square
  vv = sum(v^2)
  sd = sigma*sqrt(vv)
  # rho as defined in hyun et al.
  rho = Gv / vv # rho should be an mpfr object
  rho[abs(rho)<=tol] <- 0 # floating point round off
  vec = (u - Gy + rho*z) / rho
  vlo = suppressWarnings(max(vec[rho>0])) #remove suppress warnings
  vup = suppressWarnings(min(vec[rho<0])) #remove suppress warnings
  if(is.null(vlo)|is.null(vup)){
    vlo = -Inf
    vup = Inf
  }
  if(vlo>vup){
    stop('Error: vlow should be smaller than vup, invalid interval!')
  }
  return(list(vlo=vlo,vup=vup,sd=sd,z=z))
}

