# Auxiliary Functions
# * aux_dist               : fast distance computation using 'Rfast::Dist'
# * aux_kernel_standard    : given a distance object, construct Gaussian kernel
# * aux_compression        : kernel compression
# * aux_entropyrule        : von Neumann entropy for stopping criterion
#   aux_entropyrule_markov : why not from transition matrix ?
# * aux_fastpca            : use 'rARPACK' 
# * aux_yariv_local        : alternative to Local PCA by Yariv's algorithm
# * aux_binarynetwork      : network input
# * aux_effective          : compute effective resistance on a binary graph
# * aux_pinv               : pseudo-inverse


# aux_dist ----------------------------------------------------------------
#' @keywords internal
#' @noRd
aux_dist <- function(X){
  return(stats::as.dist(cpp_distance(X)))
}

# * aux_kernel_standard ---------------------------------------------------
#' @keywords internal
#' @noRd
aux_kernel_standard <- function(D, nbdk, alpha){
  # input
  x = as.matrix(D)
  n = base::nrow(x)
  
  # k-th nearest distance
  nndist = rep(0,n)
  for (i in 1:n){
    tgt = as.vector(x[i,])
    nndist[i] = tgt[order(tgt)[nbdk+1]]
  }
  
  # Build Kernel Matrix
  matK = array(1,c(n,n))
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      term1 = exp(-((x[i,j]/nndist[i])^(alpha)))
      term2 = exp(-((x[i,j]/nndist[j])^(alpha)))
      matK[i,j] <- matK[j,i] <- 0.5*(term1+term2)
    }
  }
  return(matK)
}

# * aux_compression -------------------------------------------------------
#' @keywords internal
#' @noRd
aux_compression <- function(P, K, M, clustering){
  # preprocessing
  M = round(M)      # number of clusters
  N = base::nrow(P) # number of observations
  if (length(clustering)!=N){
    stop("* aux_compression : wrong clustering label.")
  }
  
  # extract information
  list_indices = list()
  for (j in 1:M){
    list_indices[[j]] = which(clustering==j)
  }
  
  # construct P_NM
  P_NM = array(0,c(N,M))
  for (i in 1:N){
    for (j in 1:M){
      idj = as.vector(list_indices[[j]])
      P_NM[i,j] = base::sum(as.vector(P[i,idj]))
    }
  }
  
  # construct P_MN
  # P_MN = array(0,c(M,N))
  # for (j in 1:M){
  #   for (i in 1:N){
  #     idj  = as.vector(list_indices[[j]])
  #     summ = 0
  #     for (kk in 1:length(idj)){
  #       k = idj[kk]
  #       Qjk_term1 = base::sum(as.vector(K[k,]))
  #       Qjk_term2 = base::sum(K[idj,])
  #       Qjk = (Qjk_term1/Qjk_term2)
  #       summ = summ + (Qjk*P[k,i])
  #     }
  #     P_MN[j,i] = summ
  #   }
  # }
  P_MN = aux_compression_Q(K, M, N, list_indices)%*%P
  
  # return 
  output = list()
  output$P_MN = P_MN
  output$P_NM = P_NM
  return(output)
}
#' @keywords internal
#' @noRd
aux_compression_Q <- function(K, M, N, list_indices){
  Q = array(0,c(M,N))
  rowsumK  = base::rowSums(K)
  groupsum = rep(0, M)
  for (j in 1:M){
    groupsum[j] = base::sum(K[list_indices[[j]],])
  }
  
  for (j in 1:M){
    for (k in 1:N){
      Q[j,k] = rowsumK[k]/groupsum[j]
      # Q[j,k] = rowsumK[k]/base::sum(K[list_indices[[j]],])  
    }
  }
  padding = array(0,c(M,N))
  for (i in 1:M){
    padding[i,list_indices[[i]]]=1
  }
  return(Q*padding)
}


# * aux_entropyrule -------------------------------------------------------
#' @keywords internal
#' @noRd
aux_entropyrule <- function(K){
  vecDhalfinv = 1/sqrt(rowSums(K))
  matA = diag(vecDhalfinv)%*%K%*%diag(vecDhalfinv)
  
  eigA  = eigen(matA)$values
  eigA  = eigA[(eigA>0)]
  vec.t = 1:1000
  vec.H = rep(0,1000)
  for (i in 1:1000){
    eig.t  = (eigA^i) + (1e-7) # modified for zero-padding
    eig.t  = eig.t/base::sum(eig.t)
    term.t = -base::sum(eig.t*base::log(eig.t))
    if (is.na(term.t)){
      vec.t = vec.t[1:(i-1)]
      vec.H = vec.H[1:(i-1)]
      break
    } else {
      vec.H[i] = term.t
    }
  }
  
  fun_detect = utils::getFromNamespace("hidden_knee_clamped", "maotai")
  opt.t      = round(fun_detect(vec.t, vec.H))
  return(opt.t)
}
#' @keywords internal
#' @noRd
aux_entropyrule_markov <- function(P){
  
  eigP = base::Re(eigen(P)$values)
  eigP  = eigP[(eigP>0)]
  vec.t = 1:1000
  vec.H = rep(0,1000)
  for (i in 1:1000){
    eig.t  = (eigP^i) + (1e-7) # modified for zero-padding
    eig.t  = eig.t/base::sum(eig.t)
    term.t = -base::sum(eig.t*base::log(eig.t))
    if (is.na(term.t)){
      vec.t = vec.t[1:(i-1)]
      vec.H = vec.H[1:(i-1)]
      break
    } else {
      vec.H[i] = term.t
    }
  }
  
  fun_detect = utils::getFromNamespace("hidden_knee_clamped", "maotai")
  opt.t      = round(fun_detect(vec.t, vec.H))
  return(opt.t)
}


# * aux_fastpca -----------------------------------------------------------
#' @export
aux_fastpca <- function(data, ndim=2){
  centered   = as.matrix(base::scale(data, center=TRUE, scale=FALSE))
  projection = rARPACK::svds(centered, k=round(ndim))$v
  return(centered%*%projection)
}



# * aux_yariv_local -------------------------------------------------------
#' @export
aux_yariv_local <- function(rowdata, point, ndim, printer=TRUE){
  # setup
  d = round(ndim)
  R = t(rowdata)
  
  n = base::nrow(R) # ambient dimension
  N = base::ncol(R)
  U = rARPACK::eigs_sym(stats::cov(rowdata), k=round(ndim))$vectors # (n x d)
  
  r      = as.vector(point)
  q_old  = as.vector(point)
  Rtilde = array(0,c(n,N))
  
  # iteration
  maxiter = 100
  abstol  = 1e-7
  for (it in 1:maxiter){
    # compute R_tilde
    for (i in 1:N){
      Rtilde[,i] = as.vector(R[,i])-q_old
    }
    # X & Xtilde
    Xtilde = cbind(rep(1,N), t(Rtilde)%*%U)
    # solve for alpha
    alpha = base::solve(t(Xtilde)%*%Xtilde, t(Xtilde)%*%t(Rtilde))
    # update q
    qtilde = q_old + as.vector(alpha[1,])
    # QR Decomposition
    Q = qr.Q(qr(t(alpha[2:nrow(alpha),])))
    # update
    U     = Q
    q_new = qtilde + as.vector(U%*%t(U)%*%(r-qtilde))
    q_inc = sqrt(sum((q_old-q_new)^2))
    q_old = q_new
    if (q_inc < abstol){
      if (printer){
        print(paste0("* MLS : iteration ",it, " termination!."))
      }
      break
    }
    if (printer){
      print(paste0("* MLS : iteration ",it, " complete: inc=",round(q_inc,4)))  
    }
  }
  
  # return
  output = list()
  output$location = q_old
  output$subspace = U
  return(output)
}
#' @export
aux_projection_distance <- function(projcube){
  return(cpp_distance_proj(projcube))
}



# * aux_binarynetwork -----------------------------------------------------
#' @keywords internal
#' @noRd
aux_binarynetwork <- function(input){
  if (inherits(input, "igraph")){
    output = round(as.matrix(igraph::as_adjacency_matrix(input)))
  } else {
    output = round(as.matrix(input))
  }
  uvector = unique(as.vector(output))
  if (length(uvector)!=2){
    stop("* binary network only.")
  }
  if (!all(uvector==c(0,1))){
    stop("* binary : {0,1} values.")
  }
  return(output)
}


# * aux_pinv               : pseudo-inverse -------------------------------
#' @keywords internal
#' @noRd
aux_pinv <- function(A){
  svdA      = base::svd(A)
  tolerance = (.Machine$double.eps)*max(c(nrow(A),ncol(A)))*as.double(max(svdA$d))
  
  idxcut    = which(svdA$d <= tolerance)
  invDvec   = (1/svdA$d)
  invDvec[idxcut] = 0
  
  output = (svdA$v%*%diag(invDvec)%*%t(svdA$u))
  return(output)
}

# * aux_effective ---------------------------------------------------------
#' @export
aux_effective <- function(A){ # input matrix
  # parameters and pre-compute L
  if (inherits(A, "dgCMatrix")){
    A = as.matrix(A)
  }
  N = nrow(A)
  diag(A) = rep(0, as.integer(N))
  L = base::diag(base::rowSums(A)) - A # graph laplacian
  
  # 1. compute Q
  tmp = diag(N)
  tmp[,1] = rep(1,N)/sqrt(N)
  qrq = base::qr.Q(base::qr(tmp))
  Q = t(qrq[,2:N])
  
  # 2. compute \bar{L}
  Lbar = Q%*%L%*%t(Q)
  
  # 3. solve Lyapunov equation
  Sigma = maotai::lyapunov(Lbar, diag(rep(1,N-1)))
  
  # 4. compute X
  X = 2*(t(Q)%*%Sigma%*%Q)
  
  # 5. compute 
  out.cpp = cpp_effective(X)
  if (!isSymmetric(out.cpp)){
    out.cpp = (out.cpp + t(out.cpp))/2
  }
  return(out.cpp)
}
#' @export
aux_effectivesym <- function(A){
  n = nrow(A)
  L = diag(rowSums(A))-A
  Linv = aux_pinv(L)
  
  output = array(0,c(n,n))
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      est = rep(0,n)
      est[i] = 1
      est[j] = -1
      
      output[i,j] <- output[j,i] <-sum(as.vector(Linv%*%est)*est)
    }
  }
  return(output)
}
# star10 = igraph::make_full_citation_graph(30, directed=FALSE)
# matt10 = as.matrix(igraph::as_adjacency_matrix(star10))
# 
# aa = aux_effective(matt10)
# bb = aux_effectivesym(matt10)