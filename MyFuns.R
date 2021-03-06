rm(list = ls())
#####
if (!require('ape')) install.packages('ape')
require('ape')
if (!require('geiger')) install.packages('geiger')
require(geiger)
if (!require('adephylo')) install.packages('adephylo')
require(adephylo)

if (!require('genlasso')) install.packages('genlasso')
require('genlasso')
if (!require('truncnorm')) install.packages('truncnorm')
require(truncnorm)

#####

# plot_tree: function used to plot my tree
#' @param tree the tree to be plotted
#' @param cls the associated clusters
#' @return a list with clustered tips

plot_tree <- function(tree, cls){
  colored_edge <- NULL
  colored_tip <- NULL
  tip_color <- NULL
  
  colo <- rainbow(length(cls))
  for (i in 1:length(cls)){
    colored_tip[[i]] <- tips(tree, cls[i])
    colored_edge <- c(colored_edge, which.edge(tree,colored_tip[[i]]))
    tip_color <- c(tip_color,rep(colo[i],length(which.edge(tree,colored_tip[[i]]))))
  }
  colr <- rep("darkgrey", nrow(tree$edge)) #colors of unassociated edges
  colr[colored_edge] <- tip_color #colors of associated edges
  plot(phy,show.tip.label = F,edge.color = colr)
  nodelabels(cex = 0.5,frame = 'none',adj = c(0,0.5))
  tiplabels(cex = 0.5,frame = 'none',adj = c(0,0.5))
  return(colored_tip)
}

######

# ln_param: generate mean and vcov matrix for logistic normal
#' @param tree the phenogenetic tree structure
#' @param seed random seed - used when generate mean vector
#' @hyper_mu hyperparameter mean when generate mean
#' @hyper_sig hyperparameter sd when generate mean
#' @return a list (mu,Sig)
ln_param <- function(tree, seed = 1012, hyper_mu = 0, hyper_sig = 1){
  p <- length(tree$tip.label) #number of Tips
  
  set.seed(seed)
  mu <- rnorm(p, hyper_mu, hyper_sig)
  
  patr.dis <- cophenetic(tree)
  temp <- 1 - patr.dis/max(patr.dis) #1^t * 1 - C
  Sig <- max(-min(eigen(temp)$values),0) * diag(p) + temp
  
  return(list(mu = mu,Sig = Sig))
}

#####

# gen_OTU : generate OTU data
#' @param ln_par parameter of logistic normal dist'n list(mu,Sig)
#' @param n number of subjects
#' @return OTU compositional data
gen_OTU <- function(ln_par, n){
  require(MASS)
  W <- mvrnorm(n,ln_par$mu,ln_par$Sig)
  X <- t(apply(exp(W),1,FUN = function(x) x/sum(x)))
  return(X)
}

#####

# gen_gma : generate coefficients with given proportion of variability
# in response explained by each factor
#' @param pi gene sqrt(prob=0)
#' @param p1 proportion explained by gene
#' 
#' @param ln_par parameter of logistic normal dist'n
#' @param N number of simulation times to estimate OTU var
#' @param p2 proportion explained by microbiome
#' @param tree the phylogenetic tree
#' @param cls the associated clusters (denoted by internal node)
#' 
#' @param sig2 given variance of error term (as baseline)
#' 
#' @param p3 proportion explained by the interaction
#' @return gma : regression coefficients
gen_gma <- function(pi = 0.5, p1 = 0.02, ln_par, N = 10000, p2 = 0.05, tree, cls, sig2 = 1, p3 = 0){
  p_t = p1 + p2 + p3
  var_G <- 2 * pi * (1-pi)
  gma_G <- sqrt(p1*sig2/((1-p_t)*var_G))
  
  # set.seed(1012)
  X <- gen_OTU(ln_par,N)
  Z <- matrix(NA, nrow = nrow(X),ncol = length(cls))
  for (i in 1:length(cls)){
    tip_i <- tips(tree,cls[i])
    Z[,i] <- rowSums(X[,tip_i])
  }
  # define relative ratio of the OTU effects
  # ratio <- 1:length(cls) 
  if(length(cls) %% 2) ratio <- c(setdiff(-floor(length(cls)/2):floor(length(cls)/2),0), ceiling(length(cls)/2)) else 
    ratio <- setdiff(-floor(length(cls)/2):floor(length(cls)/2),0)
  var_M <- var(as.numeric(Z %*% ratio))
  w <- sqrt(p2*sig2/((1-p_t)*var_M))
  gma_M <- w * as.vector(ratio)
  
  #interaction
  G <- sample(x = c(0,1,2), size = N, replace = TRUE, prob = c(pi^2,2*pi*(1-pi),(1-pi)^2))
  GZ <- G * Z
  var_I <- var(as.numeric(GZ %*% ratio))
  ww <- sqrt(p3*sig2/((1-p_t)*var_I))
  gma_I <- ww * as.vector(ratio)
  
  return(list(gma_G = gma_G, gma_M = gma_M, gma_I = gma_I))
}

#####

# true_beta : return the true coef vector
#' @param tree the tree structure
#' @param cls the associated tips
#' @param gma the true coeff
#' @return beta : the true coeff with the same length as the ncol. of X (# OTUs)

true_beta <- function(tree,cls,gma){
  Ntips <- length(tree$tip.label)
  true_beta <- rep(0, Ntips)
  for (i in 1:length(cls)){
    asso_tips <- tips(tree,cls[i])
    true_beta[asso_tips] <- gma$gma_M[i]
  }
  return(true_beta)
}

######

# gen_dat : generate data (G,X,Z,y)
#' @param n the number of subjects
#' @param pi the gene prob.
#' @param ln_par - logistic normal parameters, a list with mean and cov matrix
#' @param gma regression coefficients: list(gma_G,gma_M)
#' @param tree the tree
#' @param cls the associated clusters
#' @param sig s.d. of the error term. Default to be 1.
#' @return A list
gen_dat <- function(n, pi = 0.5, ln_par, gma, tree, cls, sig = 1){
  G <- sample(x = c(0,1,2), size = n, replace = TRUE, prob = c(pi^2,2*pi*(1-pi),(1-pi)^2))
  X <- gen_OTU(ln_par,n)
  
  Z <- matrix(NA, nrow = nrow(X), ncol = length(cls))
  for (i in 1:length(cls)){
    tip_i <- tips(tree,cls[i])
    Z[,i] <- rowSums(X[,tip_i])
  }
  
  GZ <- G * Z
  eps <- rnorm(n, 0, sig)
  
  y <- G * gma$gma_G + Z %*% gma$gma_M + GZ %*% gma$gma_I + eps
  y <- c(y)
  return(list(G = G, X = X, Z = Z, y = y))
}

#####

# gen_D: generate the fusion penalty matrix from tree structure
#' @param tree the phenogenetic tree
#' @param m the power added to the weights
#' @param weight c('avarage','max','pair') how to set weights for interior node
#' @param type c('myown','wang1') how to integrate fusion pattern in D
#' @return list(D,W) a list of the penality matrix and weights

gen_D <- function(tree, m = 1, weight = 'max', type = 'myown'){
  require(geiger)
  require(phylobase)
  
  p <- length(tree$tip.label)
  N <- tree$Nnode
  Dis <- dist.nodes(tree); Dis <- Dis/max(Dis)
  
  W_v <- numeric(N)
  if (type == 'myown'){
    D <- NULL
    for (v in 1:N){
      Lv <- tips(tree, p+v)
      n_Lv <- length(Lv)
      d_v <- rep(0,p); d_v[Lv] <- 1
      BrLen_v <- Dis[p+v,Lv]
      w_v <- ifelse(weight=='max',1/max(BrLen_v)^m,1/mean(BrLen_v)^m)
      #w_v <- 1/max(BrLen_v)^m
      W_v[v] <- w_v
      for (j in Lv){
        e_j <- rep(0,p); e_j[j] <- 1
        D <- rbind(D, (e_j - d_v/n_Lv)*w_v/n_Lv)
      }
    }
  } else {
    y <- phylo4(tree)
    D <- matrix(nrow = N, ncol = p)
    for (v in 1:N){
      child_v <- as.vector(children(y,p+v))
      Lv_1 <- tips(tree,child_v[1])
      Lv_2 <- tips(tree,child_v[2])
      dv_1 <- rep(0,p); dv_1[Lv_1] <- 1/length(Lv_1)
      dv_2 <- rep(0,p); dv_2[Lv_2] <- 1/length(Lv_2)
      dv <- dv_1 - dv_2
      
      BrLen_v <- Dis[p+v,c(Lv_1,Lv_2)]
      w_v <- ifelse(weight=='max',1/max(BrLen_v)^m,1/Dis[child_v[1],child_v[2]]^m)
      W_v[v] <- w_v
      
      D[v,] <- dv * w_v
    }
  }
  return(list(D = D,Weights = W_v))
}

#####

# data_prep : generate the data and penalty matrix for penalized regression
#' @param Data the data containing G, X, y
#' @param DW the list containing the penalty matrix for fusion pattern & weights
#' @param alp the weight assigned on sparsity
#' @param ref the reference OTU level
#' @param normlz whether to normalize the patristic distances or not
#' @return G_cen, X_cen, X_cen1, y_cen, D1
data_prep <- function(Data, DW, ref, alp = 0, normlz = F){
  G_cen <- Data$G - mean(Data$G)
  y_cen <- Data$y - mean(Data$y)
  X_cen <- t(t(Data$X) - colMeans(Data$X))
  X_cen1 <- X_cen[,-ref]
  Z_cen <- t(t(Data$Z) - colMeans(Data$Z))
  D <- DW$D[,-ref]
  if (normlz) normalizer <- c(DW$Weights,NCOL(X_cen1)) else normalizer <- c(1,1)
  if (alp > 0 & alp < 1) {
    D1 <- rbind((1-alp)*D/normalizer[1], alp*diag(NCOL(X_cen1))/normalizer[2])
  } else if (alp == 1) {
    D1 <- diag(NCOL(X_cen1))
  } else D1 <- D
  return(list(y_cen = y_cen, X_cen = X_cen, G_cen = G_cen,
              X_cen1 = X_cen1, Z_cen = Z_cen, D1 = D1))
}



#######

# select_OTU : a function to select and fuse associated OTUs
#' @param y,X,D The prepared data for penalized regression
#' @param c_bic c('min','psi')
#' @param dtol the tolerant level for decreasing in BIC
#' @param svd whether to use SVD in genlasso
#' @return genlasso object, stop.index, bic values at each step and sig2
select_OTU <- function(y, X, D, maxsteps = 2000, svd = F){
  require(genlasso)
  # run genlasso
  g_lasso <- genlasso(y, X, D, maxsteps = maxsteps, svd = svd)
  # estimate sigma^2
  sig2 <- sum((g_lasso$y - X %*% g_lasso$bls)^2)/(NROW(X) - NCOL(X))
  # compute BIC
  bic <- colSums((g_lasso$y - g_lasso$fit)^2) + log(length(g_lasso$y)) * g_lasso$df * sig2
  bic <- as.vector(bic)
  # suggest which step to stop
  index <- which.min(bic)
  # add new elements to genlasso 
  g_lasso$stop.index = index
  g_lasso$bic = bic
  g_lasso$sig2 = sig2
  return(g_lasso)
}

#####

#' gen_select: a function implementing the generalized lasso
#' using SVD, which also generate the affine contraints in the meantime
#' @param y the response (n)
#' @param X the designe matrix (n,p), assumed to have full col rank
#' @param D the penalty matrix (m,p)
#' @param rtol singular values below rtol are viewed as 0; default to be sqrt(.Machine$double.eps)
#' @param maxsteps maximum number of steps, default to be 0
#' 
#' @return U : matrix of dual variables (m, maxsteps)
#' @return lambda : vector of knot values (maxsteps)
#' @return fit : the fitted response
#' @return beta : the optimal primal variable
#' @return df : nullity(D2[-B,])
#' @return h : whether hit or leave at each step
#' @return bls : the least squres fit
#' @return bic : the bic at each step
#' @return stop.index : the stopping step suggested by BIC
#' @return Gama, dd : The affine constraints: Gama %*% y <= dd
#' @return sig2 : the estimate of the error variances


gen_select <- function(y, X, D, rtol = 1e-7, btol = 1e-7, maxsteps = 2000){
  svd_D <- function(A, b, rtol = 1e-7){
    # A function to calculate inv(t(D2[B_c,]))
    # A = t(D2[B_c,])
    # inv(A) %*% b
    # nullity(D2[B_c,]) = n - rank(t(D2[B_c,]))
    svd.A <- svd(A)
    d <- svd.A$d
    r <- sum(d>rtol)
    d <- d[1:r]
    u <- svd.A$u[,1:r]; v <- svd.A$v[,1:r] # D_k <- u %*% d %*% t(v) is the condensed SVD
    
    x <- v %*% ((t(u)/d) %*% b)
    inv <- v %*% (t(u)/d)
    Proj.D_k <- u %*% t(u)
    return(list(x = x, inv = inv, Prj = Proj.D_k, r = r))
  }
  # some basic parameters
  n <- length(y)
  m <- nrow(D)
  p <- ncol(D)
  
  # do some transformations
  # reformulate the problem with general X as a signal approx. problem
  x <- svd(X)
  if (min(x$d) < rtol) stop('X does not have full col. rank!')
  y2 <- as.numeric(x$u %*% (t(x$u) %*% y))
  X_inv <- x$v %*% (t(x$u)/x$d)
  D2 <- D %*% X_inv
  bls <- y2
  sig2 <- sum((y - y2)^2)/(n-p)
  
  # Intialize things to keep track of & return with
  U <- matrix(nrow = m, ncol = maxsteps) # optimal dual variables at each step
  lambs <- numeric(maxsteps) # knot values at each step
  B <- NULL # boundary set
  B_c <- seq(1:m) # interior set
  S <- NULL # signs of coordinates in B
  bic <- numeric(maxsteps) # adjusted bic at each step
  df <- numeric(maxsteps) # df of the fit, i.e. nullity(D2[B_c,])
  h <- logical(maxsteps) # whether hit or not
  
  
  Gama <- NULL; dd <- NULL
  
  ###########Begin!###########
  
  ##First step
  temp <- svd_D(t(D2), y2)
  U[,1] <- temp$x
  i_hit <- which.max(abs(U[,1]))
  r_hit <- sign(U[i_hit,1])
  lambs[1] <- abs(U[i_hit,1])
  h[1] <- T
  df[1] <- n - temp$r
  bic[1] <- round(sum((temp$Prj %*% y2)^2) + log(n) * sig2 * df[1], 3)
  B <- c(B,B_c[i_hit]) # must hit
  
  change <- list(new = D2[B_c[i_hit],], Prj = temp$Prj)
  
  B_c <- B_c[-i_hit]
  S <- c(S,r_hit)
  Ds <- D2[i_hit,] * S # vector t(D[B,]) %*% S
  D_1 <- D2[-i_hit,,drop = F] # matrix D[B_c,]
  D_2 <- D2[i_hit,,drop = F] # matrix D[B,]

  Gama_add <- cbind(t(temp$inv[-i_hit,]) - r_hit * temp$inv[i_hit,],
                    -t(temp$inv[-i_hit,]) - r_hit * temp$inv[i_hit,])
  dd_add <- rep(0,2*(m-1))
  if (sum(t(Gama_add) %*% y2 > dd_add + btol)) stop('Impossible1~')
  Gama <- rbind(Gama, t(Gama_add))
  dd <- c(dd, dd_add)
  
  
  
  YELLOW <- F
  # flag <- 0
  k <- 2
  
  while(k <= maxsteps & lambs[k-1] > 0){
    
    temp <- svd_D(t(D_1), cbind(y2,Ds), rtol = 1e-7)
    
    df[k] <- n - temp$r
    
    # decide when to quit loop based on BIC
    if (df[k] == df[k-1]) {
      bic[k] <- bic[k-1]
    } else {
      bic[k] <- round(sum((temp$Prj %*% y2)^2) + log(n) * sig2 * df[k], 3)
      if (bic[k] > bic[k-1]) {
        # bic is increasing
        if (h[k-1]) {
          vec <- as.vector((diag(n) - temp$Prj) %*% change$new)
          norm.vec <- vec/sqrt(sum(vec^2))
          Gama_add <- rbind(norm.vec,-norm.vec)
          dd_add <- rep(sqrt(sig2 * log(n)),2)
        } else {
          vec <- as.vector((diag(n) - change$Prj) %*% change$new)
          norm.vec <- vec/sqrt(sum(vec^2))
          Gama_add <- -sign(c(norm.vec %*% y2)) * norm.vec
          dd_add <- -sqrt(sig2 * log(n))
        }
        if (YELLOW == T) {
          stop.index <- k - 1
          break
        } else YELLOW <- T
      } else {
        # bic is decreasing
        if (h[k-1]) {
          vec <- as.vector((diag(n) - temp$Prj) %*% change$new)
          norm.vec <- vec/sqrt(sum(vec^2))
          Gama_add <- -sign(c(norm.vec %*% y2)) * norm.vec
          dd_add <- -sqrt(sig2 * log(n))
        } else {
          vec <- as.vector((diag(n) - change$Prj) %*% change$new)
          norm.vec <- vec/sqrt(sum(vec^2))
          Gama_add <- rbind(norm.vec,-norm.vec)
          dd_add <- rep(sqrt(sig2 * log(n)),2)
        }
        YELLOW <- F
       # flag <- k
      }
    # if (sum(Gama_add %*% y2 > dd_add + btol)) stop('Impossible2~')
    # Gama <- rbind(Gama, Gama_add)
    # dd <- c(dd,dd_add)
    }
    
    
    
    
    # hitting times
    a = as.numeric(temp$x[,1])
    b = as.numeric(temp$x[,2])
    
    R <- sign(a)
    Gama_add <- -R * temp$inv
    dd_add <- rep(0, nrow(Gama_add))
    if (sum(Gama_add %*% y2 > dd_add + btol)) stop('Impossible2+~')
    Gama <- rbind(Gama, Gama_add)
    dd <- c(dd, dd_add)
    
    hits <- a/(R+b)
    hits[hits > lambs[k-1] + btol] <- 0
    hits[hits > lambs[k-1]] <- lambs[k-1]
    
    i_hit <- which.max(hits)
    hit <- hits[i_hit]
    r_hit <- R[i_hit]
    
    Gama_add <- t(temp$inv[-i_hit,,drop = F]/(R[-i_hit] + b[-i_hit])) - temp$inv[i_hit,]/(r_hit + b[i_hit])
    dd_add <- rep(0, nrow(t(Gama_add)))
    if (sum(t(Gama_add) %*% y2 > dd_add + btol)) stop('Impossible3~')
    Gama <- rbind(Gama, t(Gama_add))
    dd <- c(dd, dd_add)
    
    
    # leaving times
    c <- S * (D_2 %*% (y2 - t(D_1) %*% a))
    d <- S * (D_2 %*% (Ds - t(D_1) %*% b))
    
    
    I <- which(c < -btol & d < -btol)
    I_c <- which(c > -btol & d < -btol)
    
    C.mat <- (S * D_2 %*% (diag(n) - temp$Prj))
    Gama_add <- rbind(C.mat[I,], -C.mat[I_c,])
    dd_add <- rep(0, nrow(Gama_add))
    if (sum(Gama_add %*% y2 > dd_add + btol)) stop('Impossible4~')
    Gama <- rbind(Gama, Gama_add)
    dd <- c(dd,dd_add)
    
    leaves <- c/d
    leaves[c > -btol] <- 0
    leaves[leaves > lambs[k-1] + btol] <- 0
    leaves[leaves > lambs[k-1]] = lambs[k-1]
    
    i_leave <- which.max(leaves)
    leave <- leaves[i_leave]
    
    if (leave){
      I_i <- setdiff(I,i_leave)
      D.vec <- C.mat %*% Ds
      Gama_add <- t(C.mat[I_i,,drop = F]/D.vec[I_i]) - C.mat[i_leave,]/D.vec[i_leave]
      dd_add <- rep(0, nrow(t(Gama_add)))
      if (sum(t(Gama_add) %*% y2 > dd_add + btol)) stop('Impossible5~')
      Gama <- rbind(Gama, t(Gama_add))
      dd <- c(dd, dd_add)
      
      Gama_add <- -temp$inv[i_hit,]/(r_hit + b[i_hit]) + C.mat[i_leave,]/D.vec[i_leave] #if hit
      dd_add <- 0
    }
    
    # compare between hit and leave
    if (hit > leave) {
      lambs[k] <- hit
      h[k] <- T
      
      uhat <- numeric(m)
      uhat[B] <- hit*S
      uhat[B_c] <- a - hit*b
      U[,k] <- uhat
      
      change <- list(new = D_1[i_hit,], Prj = temp$Prj)
      
      B <- c(B,B_c[i_hit])
      B_c <- B_c[-i_hit]
      S <- c(S, r_hit)
      
      Ds <- Ds + D_1[i_hit,] * r_hit
      D_2 <- rbind(D_2,D_1[i_hit,])
      D_1 <- D_1[-i_hit,,drop = F]
      
      if (leave) {
        if (sum(Gama_add %*% y2 > dd_add + btol)) stop('Impossible6~')
        Gama <- rbind(Gama, Gama_add)
        dd <- c(dd,dd_add)
      }
    } else {
      lambs[k] <- leave
      h[k] <- F
      uhat <- numeric(m)
      uhat[B] <- leave*S
      uhat[B_c] <- a - leave*b
      U[,k] <- uhat
      
      change <- list(new = D_2[i_leave,], Prj = temp$Prj)
      
      B_c <- c(B_c, B[i_leave])
      B <- B[-i_leave]
      Ds <- Ds - D_2[i_leave,] * S[i_leave]
      S <- S[-i_leave]
      
      D_1 <- rbind(D_1,D_2[i_leave,])
      D_2 <- D_2[-i_leave,,drop = F]
      
      if (leave) {
        Gama_add <- -Gama_add
        if (sum(Gama_add %*% y2 > dd_add + btol)) stop('Impossible7~')
        Gama <- rbind(Gama, Gama_add)
        dd <- c(dd,dd_add)
      }
    }
    k <- k + 1
    # if (k>=400) break
    if (!k%%30) print(k)
  }

  #trim
  lambs <- lambs[1:(k-1)]
  U <- U[,1:(k-1), drop = F]
  df <- df[1:(k-1)]
  h <- h[1:(k-1)]
  bic <- bic[1:(k-1)]
  
  fit <- y2 - t(D2) %*% U
  beta <- X_inv %*% fit
  gama <- Gama %*% x$u %*% t(x$u)
  
  return(list(lambda = lambs, beta = beta, fit = fit, sig2 = sig2,
              U = U, df = df, h = h, bls = bls, bic = bic,
              stop.index = stop.index, Gama = Gama, gama = gama, dd = dd))
}

######
# Another version without generating gamma matrix
gen_select2 <- function(y, X, D, rtol = 1e-7, btol = 1e-7, maxsteps = 2000){
  svd_D <- function(A, b, rtol = 1e-7){
    # A function to calculate inv(t(D2[B_c,]))
    # A = t(D2[B_c,])
    # inv(A) %*% b
    # nullity(D2[B_c,]) = n - rank(t(D2[B_c,]))
    svd.A <- svd(A)
    d <- svd.A$d
    r <- sum(d>rtol)
    d <- d[1:r]
    u <- svd.A$u[,1:r]; v <- svd.A$v[,1:r] # D_k <- u %*% d %*% t(v) is the condensed SVD
    
    x <- v %*% (1/d * (t(u) %*% b))
    # inv <- v %*% (t(u)/d)
    Proj.D_k <- u %*% t(u)
    return(list(x = x, Prj = Proj.D_k, r = r))
  }
  # some basic parameters
  n <- length(y)
  m <- nrow(D)
  p <- ncol(D)
  
  # do some transformations
  # reformulate the problem with general X as a signal approx. problem
  x <- svd(X)
  if (min(x$d) < rtol) stop('X does not have full col. rank!')
  y2 <- as.numeric(x$u %*% (t(x$u) %*% y))
  X_inv <- x$v %*% (t(x$u)/x$d)
  D2 <- D %*% X_inv
  bls <- c(y2)
  sig2 <- sum((y - y2)^2)/(n-p)
  
  # Intialize things to keep track of & return with
  U <- matrix(nrow = m, ncol = maxsteps) # optimal dual variables at each step
  lambs <- numeric(maxsteps) # knot values at each step
  B <- NULL # boundary set
  B_c <- seq(1:m) # interior set
  S <- NULL # signs of coordinates in B
  bic <- numeric(maxsteps) # adjusted bic at each step
  # bic_c2 <- numeric(maxsteps)
  bic_n <- numeric(maxsteps)
  df <- numeric(maxsteps) # df of the fit, i.e. nullity(D2[B_c,])
  h <- logical(maxsteps) # whether hit or not

  
  ###########Begin!###########
  
  ##First step
  temp <- svd_D(t(D2), y2)
  U[,1] <- temp$x
  i_hit <- which.max(abs(U[,1]))
  r_hit <- sign(U[i_hit,1])
  lambs[1] <- abs(U[i_hit,1])
  h[1] <- T
  df[1] <- n - temp$r
  # bic[1] <- round(sum((temp$Prj %*% y2)^2) + log(n) * sig2 * df[1], 3)
  # bic_c2[1] <- round(sum((temp$Prj %*% y2)^2) + log(n) * sig2 * df[1], 3)
  bic_n[1] <- round(sum((temp$Prj %*% y2)^2) + log(n) * sig2 * df[1], 3)
  B <- c(B,B_c[i_hit]) # must hit
  
  # change <- list(new = D2[B_c[i_hit],], Prj = temp$Prj)
  
  B_c <- B_c[-i_hit]
  S <- c(S,r_hit)
  Ds <- D2[i_hit,] * S # vector t(D[B,]) %*% S
  D_1 <- D2[-i_hit,,drop = F] # matrix D[B_c,]
  D_2 <- D2[i_hit,,drop = F] # matrix D[B,]
  
  
  
  
  YELLOW <- F
  flag <- 1
  k <- 2
  
  while(k <= maxsteps & lambs[k-1] > 0){
    
    temp <- svd_D(t(D_1), cbind(y2,Ds), rtol = 1e-7)
    
    df[k] <- n - temp$r
    
    # decide when to quit loop based on BIC
    # if (df[k] == df[k-1]) {
      # bic[k] <- bic[k-1]
    # } else {
      # bic[k] <- round(sum((temp$Prj %*% y2)^2) + log(n) * sig2 * df[k], 3)
      # if (bic[k] > bic[k-1]) {
      #   # bic is increasing
      #   if (YELLOW == T) {
      #     stop.index <- flag
      #     # break
      #   } else YELLOW <- T
      # } else {
      #   # bic is decreasing
      #   YELLOW <- F
      #   flag <- k
      # }
    # }
    
    
    # hitting times
    a = as.numeric(temp$x[,1])
    b = as.numeric(temp$x[,2])
    
    R <- sign(a)
    
    hits <- a/(R+b)
    hits[hits > lambs[k-1] + btol] <- 0
    hits[hits > lambs[k-1]] <- lambs[k-1]
    
    i_hit <- which.max(hits)
    hit <- hits[i_hit]
    r_hit <- R[i_hit]
    
    
    
    # leaving times
    c <- S * (D_2 %*% (y2 - t(D_1) %*% a))
    d <- S * (D_2 %*% (Ds - t(D_1) %*% b))
    
    leaves <- c/d
    leaves[c > -btol] <- 0
    leaves[leaves > lambs[k-1] + btol] <- 0
    leaves[leaves > lambs[k-1]] <- lambs[k-1]
    
    i_leave <- which.max(leaves)
    leave <- leaves[i_leave]
    
    # compare between hit and leave
    if (hit > leave) {
      lambs[k] <- hit
      # bic_c2[k] <- round(sum((temp$Prj %*% y2 + hit * (diag(n) - temp$Prj) %*% Ds)^2) + log(n) * sig2 * df[k], 3)
      
      h[k] <- T
      
      uhat <- numeric(m)
      uhat[B] <- hit*S
      uhat[B_c] <- a - hit*b
      U[,k] <- uhat
      
      bic_n[k] <- round(sum((t(D2) %*% uhat)^2) + log(n) * sig2 * df[k],3)
      
      B <- c(B,B_c[i_hit])
      B_c <- B_c[-i_hit]
      Ds <- Ds + D_1[i_hit,] * r_hit
      S <- c(S, r_hit)
      
      D_2 <- rbind(D_2,D_1[i_hit,])
      D_1 <- D_1[-i_hit,,drop = F]
      
    } else {
      lambs[k] <- leave
      # bic_c2[k] <- round(sum((temp$Prj %*% y2 + leave * (diag(n) - temp$Prj) %*% Ds)^2) + log(n) * sig2 * df[k], 3)
      
      h[k] <- F
      uhat <- numeric(m)
      uhat[B] <- leave*S
      uhat[B_c] <- a - leave*b
      U[,k] <- uhat
      
      bic_n[k] <- round(sum((t(D2) %*% uhat)^2) + log(n) * sig2 * df[k],3)
      
      B_c <- c(B_c, B[i_leave])
      B <- B[-i_leave]
      Ds <- Ds - D_2[i_leave,] * S[i_leave]
      S <- S[-i_leave]
      
      D_1 <- rbind(D_1,D_2[i_leave,])
      D_2 <- D_2[-i_leave,,drop = F]
    }
    # decide whether to stop
    if (bic_n[k] > bic_n[k-1]) {
      # bic is increasing
      if (YELLOW == T) {
        # stop.index <- flag
        # break
      } else YELLOW <- T
    } else if (bic_n[k] < bic_n[k-1]){
      # bic is decreasing
      YELLOW <- F
      flag <- k
    }
    
    
    k <- k + 1
    if (!k%%30) print(k)
  }
  
  #trim
  lambs <- lambs[1:(k-1)]
  U <- U[,1:(k-1), drop = F]
  df <- df[1:(k-1)]
  h <- h[1:(k-1)]
  # bic <- bic[1:(k-1)]
  # bic_c2 <- bic_c2[1:(k-1)]
  bic_n <- bic_n[1:(k-1)]
  stop.index <- which.min(bic_n)
  
  fit <- y2 - t(D2) %*% U
  beta <- X_inv %*% fit
  
  # bic_c <- colSums((c(y2) - fit)^2) + log(n) * sig2 * df
  
  return(list(lambda = lambs, beta = beta, fit = fit, sig2 = sig2, bic_n = bic_n,
              U = U, df = df, h = h, bls = bls, bic = bic, stop.index = stop.index))
}

######

#' a function to calculate the IC criterion
#' 
#' @param obj a genlasso object containing beta, fit, df, bls
#' @param X the design matrix
#' @param y the response vector
#' @param type c('aic','bic','gic'), corresponding to 2, log(n), log(log(n))*log(n)
#' @return IC the total IC calculated from y, X, beta, df
IC <- function(obj, X, y, type = 'bic'){
  if (type == 'bic'){
    IC <- colSums((y - X %*% obj$beta)^2) + log(n) * obj$sig2 * obj$df
  } else if (type == 'aic'){
    IC <- colSums((y - X %*% obj$beta)^2) + 2 * obj$sig2 * obj$df
  } else if (type == 'gic'){
    IC <- colSums((y - X %*% obj$beta)^2) + log(n) * log(log(n)) * obj$sig2 * obj$df
  } else {
    stop('Not a valid IC type!')
  }
  return(IC)
}

######

# form_new : form new covariates according to estimated beta's
#' @param X_cen the (centered) design matrix (input of genlasso)
#' @param beta_esti the estimated beta
#' @param ref the reference level (whose coeff is 0)
#' @param rdig the variance of beta_esti below this level (1e-rdig) would be treated as the same
#' @return fused_OTU : a list with each element corresponds to the index of OTUs being fused
#' @return X_new: the reformed OTUs
form_new <- function(X_cen, beta_esti, rdig = 6, ref = 1){
  approx_beta <- round(beta_esti, rdig)
  level_beta <- unique(approx_beta)
  
  ref_new <- which(level_beta == approx_beta[ref])
  X1 <- sapply(level_beta, FUN = function(x) rowSums(X_cen[,approx_beta == x,drop = F]))
  X_new <- X1[,-ref_new]
  
  fused_OTU <- lapply(level_beta, FUN = function(x) sort(which(approx_beta == x)))
  fused_OTU[[ref_new]] <- NULL
  return(list(fused_OTU = fused_OTU, X_new = X_new))
}





######

# assess_fuse : Assess the model's performance in estimating fusion pattern
# Compare between the beta_est and beta_tru in their fusion pattern
#' @param tree the tree structure
#' @param beta_est the full estimated coefficient vector
#' @param beta_tru the true coefficient vector
#' @param ftol the tolerent level for whether being 0
#' @return FPR : False positive rate w.r.t internal nodes, positive means being fused together;
#'               This is analogous to Type I error with H_0: not fused together, i.e. falsely fused/really not fused
#' @return FNR : False negative rate w.r.t. internal nodes, analogous to Type II error, falsely not fused/really fused
assess_fuse <- function(tree, beta_esti, beta_true, ftol = 1e-6){
  Ntips <- length(tree$tip.label)
  
  nFPR_d <- nFPR_n <- 0
  nFNR_d <- nFNR_n <- 0
  for (i in 1:tree$Nnode){
    int_node <- Ntips + i
    tips_i <- tips(tree, int_node)
    range_tru <- diff(range(beta_true[tips_i])) > ftol
    range_est <- diff(range(beta_esti[tips_i])) > ftol
    if (range_tru) {
      # really not fused
      nFPR_d <- nFPR_d + 1
      if (!range_est) nFPR_n <- nFPR_n + 1
    } else {
      nFNR_d <- nFNR_d + 1
      if (range_est) nFNR_n <- nFNR_n + 1
    }
  }
  
  lFPR_d <- lFPR_n <- 0
  lFNR_d <- lFNR_n <- 0
  # parwise difference
  diff_tru <- abs(outer(beta_true, beta_true, '-'))
  diff_est <- abs(outer(beta_esti, beta_esti, '-'))
  leaf_tru <- diff_tru[upper.tri(diff_tru)]
  leaf_est <- diff_est[upper.tri(diff_est)]
  # compute FPR and FNR
  lFP <- leaf_tru > ftol
  lFPR <- sum(leaf_est[lFP] <= ftol) / sum(lFP)
  lFN <- !lFP
  lFNR <- sum(leaf_est[lFN] > ftol) / sum(lFN)
  return(list(nFPR = nFPR_n/nFPR_d, nFNR = nFNR_n/nFNR_d, lFPR = lFPR, lFNR = lFNR))
}

######

# assess_sparse : Assess the model's performance in estimating sparsity pattern
# Compare between the beta_est and beta_tru in their sparsity pattern
#' @param beta_est the estimated coefficient vector
#' @param beta_tru the true coefficient vector
#' @param ftol the tolerent level for whether being 0
#' @param ref the reference level
#' @return FPR : False positive rate w.r.t internal nodes, positive means being zero;
#'               This is analogous to Type I error with H_0: beta_i != 0, i.e. falsely zero/really nonzero
#' @return FNR : False negative rate w.r.t. internal nodes, analogous to Type II error, falsely nonzero/really zero
assess_sparse <- function(beta_esti, beta_true, ftol = 1e-6){
  zeros <- abs(beta_true) <= ftol
  nonzeros <- !zeros
  FPR <- sum(abs(beta_esti[nonzeros]) <= ftol)/sum(nonzeros)
  FNR <- sum(abs(beta_esti[zeros]) > ftol)/sum(zeros)
  return(list(FPR = FPR, FNR = FNR))
}



#####

# plot_beta_bic : plot the bic and beta's
#' @param beta_true the true coefficient vector
#' @param beta_esti the estimated coefficient vector
#' @param bic the reference level
plot_beta_bic <- function(beta_true, beta_esti, bic){
  opar <- par(no.readonly = T)
  par(mfrow = c(1,2))
  plot(bic,main = 'The BIC plot')
  matplot(cbind(beta_esti,beta_true), col = c('black','red'), pch = c(1,2),
          ylab = 'beta',xlab = 'index', main = 'The coefficients')
  legend(30, -0.5, c('Est.','True'), pch = c(1,2),col = c('black','red'))
  par(opar)
}


#####

# esti_beta : recover the full coef vec from the reference-level omitted version

esti_beta <- function(beta_est, ref = 1) append(beta_est, 0, after = ref-1)

#####

# get_dir : get the direction we're interested in
#' @param y_cen the response
#' @param X_new the reformed OTU
#' @param G_cen the centrered genetic factor
#' @param sig2 the sigma squared
#' @param type c('single','interact')

#' @return y_norm, eta_norm, test: if type == "interact", then eta_norm would be a single vector;
#'          o./w. a matrix with each row corresponds to a OTU factor.
#'          deg : the degree of freedom.

get_dir <- function(y_cen, X_new, G_cen, type = 'interact'){
  X_new1 <- cbind(G_cen,X_new)
  # y_norm <- c(y_cen/sqrt(sig2))
  if (type == 'interact'){
    Intact <- G_cen * X_new # interaction terms
    Intact_cen <- t(t(Intact) - colMeans(Intact)) # get each columns centered
    rsd <- lsfit(X_new1, Intact_cen, intercept = F)$residuals # modulo the information over main effects
    prj <- lsfit(rsd, y_cen, intercept = F) # projection onto interaction after modulo the main effects
    eta <- c(rsd %*% prj$coef) 
    # eta_norm <- eta / sqrt(sum(eta^2))
    deg <- prj$qr$rank
    test <- c(y_cen %*% eta)
    return(list(eta = eta, test = test, deg = deg))
  } else {
    eta <- ginv(X_new1)[2:ncol(X_new1),]
    # eta_norm <- apply(eta, 1, FUN = function(x) x/sqrt(sum(x^2)))
    test <- c(eta %*% y_cen)
    return(list(eta = eta, test = test))
  }
}


######

# get_bound : get the upper and lower bound along a certain direction
#' @param y the response
#' @param eta the direction we're interested in
#' @param obj the objective; new_res$fused_OTU
#' @param sd the initial learning rate
#' @param stop stopping criterion; default to be 0.01 (relative change in estimated bound)
#' @param max the max bound*sig expored, otherwise, set to be Inf
#' @param upper the sign of the direction is positive or negative; default to be positive
#' @return bnd the bnd for the test stat

get_bound <- function(y, eta, sd, stop = 0.01, max = 10, upper = T, verbose = F){
  dir <- eta/norm.v(eta)
  bnd <- 0
  s <- ifelse(upper, 1, -1)
  lr <- sd
  y_new <- y
  while (abs(bnd) < max * sd){
    y_new <- y_new + s * lr * dir
    
    atmpt <- select_OTU(y_new, X_cen1, D1)
    atmpt_est <- esti_beta(atmpt$beta[,atmpt$stop.index])
    if (verbose == T) plot_beta_bic(beta_true, atmpt_est, atmpt$bic)
    atmpt_res <- form_new(model.data$X_cen, atmpt_est)
    
    if (setequal(atmpt_res$fused_OTU,new_res$fused_OTU)) {
      bnd <- bnd + s * lr
    } else {
      y_new <- y_new - s * lr * dir
      lr <- lr/2
      if (lr < stop * sd) return(bnd)
    }
  }
  bnd <- ifelse(upper, Inf, -Inf)
  return(bnd)
}

######

# trunc_test : run the truncated test based on the results of get_bound
#' @param test the test statistic
#' @param bnd A vector with two elements c(lower, upper)
#' @param type can be chosen from c("norm", "chi")
#' @param deg the degree of freedom for chi test
#' @param side one-sided or two-sided test c('one','two')
#' @return p-value

trunc_test <- function(test, bnd, type = 'chi', deg = NULL, sig2 = NULL, side = 'one'){
  if (type == 'chi'){
    Test <- test / sig2
    p_val <- diff(pchisq(c(Test,bnd[2]^2), df = deg))/diff(pchisq(bnd^2, df = deg))
  } else {
    temp1 <- ptruncnorm(test, bnd[1], bnd[2], 0, sig)
    if (side == 'one'){
      p_val <- min(temp1, 1 - temp1)
    } else p_val <- 2*min(temp1, 1 - temp1)
  }
  return(p_val)
}

######

# fool_test : use all the OTUs without selecting the OTUs
#' @param data an object created by function data_prep
#' @return p_values

fool_test <- function(data){
  y_cen <- data$y_cen
  G_cen <- data$G_cen
  X_cen1 <- data$X_cen1
  GZ_cen <- cbind(G_cen, X_cen1)
  
  fit <- lm(y_cen ~ G_cen * X_cen1)
  aov_res <- anova(fit)
  p_val <- aov_res$`Pr(>F)`[3]
  return(p_val)
}

#####

# oracle_test : use the real associated OTU to do the testing
#' @param data an object created by function data_prep
#' @return p_values

oracle_test <- function(data){
  y_cen <- data$y_cen
  G_cen <- data$G_cen
  Z_cen <- data$Z_cen
  GZ_cen <- cbind(G_cen, Z_cen)
  
  fit <- lm(y_cen ~ G_cen * Z_cen)
  aov_res <- anova(fit)
  p_val <- aov_res$`Pr(>F)`[3]
  return(p_val)
}

# naive_test : use the selected OTU to do the testing but without truncation
#' @param test the testing statistic
#' @return p_values

naive_test <- function(test, sig2 = NULL, deg = NULL){
  Test <- test / sig2
  p_val <- 1- pchisq(Test, df = deg)
  return(p_val)
}





#####

# test_norm : post-selection inference for some parameter eta^t * mu
#            Assumptions: y \sim N(mu,sig2*I) and the selection event {y: Gama * y <= dd}
#' @param y observed response
#' @param Gama used for selection event
#' @param dd used for selection event
#' @param eta the linear combination we're interested in
#' @param sig2 the estimate of error variance
#' @param type c('test','CI')
#' @param lsig the significance level of C.I., default to be 0.05
#' @return p-value if type == 'test'; OR a confidence interval if type == 'CI'
test_norm <- function(y, Gama, dd, eta, sig2, type = 'test', lsig = 0.05, btol = 1e-7){
  # compute the lower and upper boundary
  y_norm <- c(y)/sqrt(sig2) # with sd = 1
  eta_norm <- eta / norm.v(eta)
  alp <- c(Gama %*% eta_norm)
  temp <- (dd - Gama %*% y_norm + alp * c(eta_norm %*% y_norm)) / alp
  V_up <- min(temp[alp > btol])
  V_lo <- max(temp[alp < -btol])
  if (V_lo >= V_up) return(99)
  # testing whether eta^t mu == 0
  require(truncnorm)
  if (type == 'test'){
    temp1 <- ptruncnorm(c(eta_norm %*% y_norm), V_lo, V_up, 0, 1)
    Test.stat <- 2 * min(temp1, 1 - temp1)
    return(Test.stat)
  } 
  # else {
  #   mu_tnorm <- function(mu) ptruncnorm(eta %*% y, V_lo, V_up, mu, sqrt(sig2 * sum(eta^2)))
  #   U.search <- V_up # min(V_up, eta %*% y + 5 * sqrt(sig2))
  #   L.search <- V_lo # min(V_lo, eta %*% y - 3 * sqrt(sig2 * sum(eta^2)))
  #   CI_up <- uniroot(function(x) mu_tnorm(x) - lsig/2, c(L.search,U.search))$root
  #   CI_lo <- uniroot(function(x) mu_tnorm(x) - (1-lsig/2), c(L.search,U.search))$root
  #   return(c(CI_lo, CI_up))
  # }
}



######

test_norm2 <- function(y, eta, sig2){
  y <- c(y)
  Test.stat <- c(y %*% eta) / sqrt(sig2 * sum(eta^2))
  p_val <- pnorm(Test.stat, 0, 1)
  return(2 * min(p_val, 1 - p_val))
}


######

# test_chi : post selection inference for a bunch of predictors; generally, under the same 
# assumptions of test_chi, we're interested in whether P * mu = 0 where P is a projection
#' @param y the observed response
#' @param Gama
#' @param dd used for selection event
#' @param P a projection matrix
#' @param sig2 the estimate of error variance
#' @return the p-value

test_chi <- function(y, Gama, dd, P, sig2, btol = 1e-7){
  y_norm <- c(y)/sqrt(sig2)
  r <- sum(diag(P))
  R <- P %*% y_norm
  R.norm <- norm.v(R)
  u <- R/R.norm
  v <- y_norm - R
  Test.stat <- R.norm
  
  alp <- Gama %*% u
  temp <- (dd - Gama %*% v)/alp
  V_up <- min(temp[alp > btol])
  V_lo <- max(c(temp[alp < -btol],0))
  if (Test.stat > V_up | Test.stat < V_lo) return(99)
  p_val <- diff(pchisq(c(Test.stat^2,V_up^2), df = r))/diff(pchisq(c(V_lo^2,V_up^2), df = r))
  return(p_val)
}


#####

test_chi2 <- function(y, P, sig2){
  y <- c(y)
  r <- sum(diag(P))
  R <- c(P %*% y)
  R.norm2 <- sum(R^2)
  Test.stat <- R.norm2/sig2

  p_val <- 1 - pchisq(Test.stat, df = r)
  return(p_val)
}


# proj.mat calculates the projection matrix onto the columns of a matrix
proj.mat <- function(A, rtol = 1e-7){
  svd.A <- svd(A)
  d <- svd.A$d
  r <- sum(d>rtol)
  u <- svd.A$u[,1:r,drop = F]
  return(u %*% t(u))
}


norm.v <- function(x) sqrt(sum(x^2))


