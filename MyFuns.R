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
#' @return gma : regression coefficients
gen_gma <- function(pi = 0.5, p1 = 0.02, ln_par, N = 10000, p2 = 0.05, tree, cls, sig2 = 1){
  var_G <- 2 * pi * (1-pi)
  gma_G <- sqrt(p1*sig2/((1-p1-p2)*var_G))
  
  # set.seed(1012)
  X <- gen_OTU(ln_par,N)
  Z <- matrix(NA, nrow = nrow(X),ncol = length(cls))
  for (i in 1:length(cls)){
    tip_i <- tips(tree,cls[i])
    Z[,i] <- rowSums(X[,tip_i])
  }
  # define relative ratio of the OTU effects
  # ratio <- 1:length(cls) 
  if(length(cls) %% 2) ratio <- c(setdiff(-floor(length(cls)/2):floor(length(cls)/2),0), ceiling(length(cls)/2)) else ratio <- setdiff(-floor(length(cls)/2):floor(length(cls)/2),0)
  var_M <- var(as.numeric(Z %*% ratio))
  w <- sqrt(p2*sig2/((1-p1-p2)*var_M))
  gma_M <- w * as.vector(ratio)
  
  return(list(gma_G = gma_G, gma_M = gma_M))
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
  Z <- matrix(NA, nrow = nrow(X),ncol = length(cls))
  for (i in 1:length(cls)){
    tip_i <- tips(tree,cls[i])
    Z[,i] <- rowSums(X[,tip_i])
  }
  eps <- rnorm(n, 0, sig)
  
  y <- G * gma$gma_G + Z %*% gma$gma_M + eps
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
  D <- DW$D[,-ref]
  if (normlz) normalizer <- c(DW$Weights,NCOL(X_cen1)) else normalizer <- c(1,1)
  if (alp > 0 & alp < 1) {
    D1 <- rbind((1-alp)*D/normalizer[1], alp*diag(NCOL(X_cen1))/normalizer[2])
  } else if (alp == 1) {
    D1 <- diag(NCOL(X_cen1))
  } else D1 <- D
  return(list(y_cen = y_cen, X_cen = X_cen, G_cen = G_cen,
              X_cen1 = X_cen1, D1 = D1))
}



#######

# select_OTU : a function to select and fuse associated OTUs
#' @param y,X,D The prepared data for penalized regression
#' @param c_bic c('min','psi')
#' @param dtol the tolerant level for decreasing in BIC
#' @param svd whether to use SVD in genlasso
#' @return genlasso object, stop index and bic values at each step
select_OTU <- function(y, X, D, c_bic = 'min', dtol = 0.1, svd = F){
  require(genlasso)
  # run genlasso
  g_lasso <- genlasso(y, X, D, svd = svd)
  # estimate sigma^2
  sig2 <- sum((g_lasso$y - X %*% g_lasso$bls)^2)/(NROW(X) - NCOL(X))
  # compute BIC
  bic <- colSums((g_lasso$y - g_lasso$fit)^2) + log(length(g_lasso$y)) * g_lasso$df * sig2
  bic <- as.vector(bic)
  # suggest which step to stop
  bic_1 <- (bic[-1] - bic[-length(bic)]) > dtol
  index <- ifelse(c_bic == 'min', which.min(bic), min(which(bic_1[-length(bic_1)] & bic_1[-1])))
  return(list(g_lasso = g_lasso, stop_index = index, bic = bic))
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
  
  change <- list(new = D2[B_c[i_hit],], Prj = temp$Prj, s = r_hit)
  
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
  flag <- 0
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
          stop.index <- flag
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
        flag <- k - 1
      }
    if (sum(Gama_add %*% y2 > dd_add + btol)) stop('Impossible2~')
    Gama <- rbind(Gama, Gama_add)
    dd <- c(dd,dd_add)
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
      
      change <- list(new = D_1[i_hit,], Prj = temp$Prj, s = r_hit)
      
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
      
      change <- list(new = D_2[i_leave,], Prj = temp$Prj, s = S[i_leave])
      
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
    # if (k==43) stop('check')
    if (!k%%30) print(k)
  }

  #trim
  lambs <- lambs[1:stop.index]
  U <- U[,1:stop.index, drop = F]
  df <- df[1:stop.index]
  h <- h[1:stop.index]
  bic <- bic[1:stop.index]
  
  fit <- y2 - t(D2) %*% U
  beta <- X_inv %*% fit
  gama <- Gama %*% x$u %*% t(x$u)
  
  return(list(lambda = lambs, beta = beta, fit = fit, sig2 = sig2,
              U = U, df = df, h = h, bls = bls, bic = bic,
              stop.index = stop.index, Gama = Gama, gama = gama, dd = dd))
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
#' @return FPR : False positive rate w.r.t internal nodes, positive means being nonzero;
#'               This is analogous to Type I error with H_0: beta_i = 0, i.e. falsely nonzero/really zero
#' @return FNR : False negative rate w.r.t. internal nodes, analogous to Type II error, falsely zero/really nonzero
assess_sparse <- function(beta_esti, beta_true, ftol = 1e-6){
  zeros <- abs(beta_true) <= ftol
  nonzeros <- !zeros
  FPR <- sum(abs(beta_esti[zeros]) > ftol)/sum(zeros)
  FNR <- sum(abs(beta_esti[nonzeros]) <= ftol)/sum(nonzeros)
  return(list(FPR = FPR, FNR = FNR))
}



#####

# plot_beta : plot the bic and beta's
#' @param beta_true the true coefficient vector
#' @param beta_esti the estimated coefficient vector
#' @param bic the reference level
plot_beta <- function(beta_true, beta_esti, bic){
  opar <- par(no.readonly = T)
  par(mfrow = c(1,2))
  plot(bic,main = 'The BIC plot')
  matplot(cbind(beta_esti,beta_true), col = c('black','red'), pch = c(1,2),
          ylab = 'beta',xlab = 'index', main = 'The coefficients')
  legend(30, 1.5, c('Est.','True'), pch = c(1,2),col = c('black','red'))
  par(opar)
}

#####
# esti_beta : recover the full coef vec from the reference-level omitted version
esti_beta <- function(beta_est, ref) append(beta_est, 0, after = ref-1)


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
  alp <- c(Gama %*% eta) / sum(eta^2)
  temp <- (dd - Gama %*% y + alp * c(eta %*% y)) / alp
  V_up <- min(temp[alp > btol])
  V_lo <- max(temp[alp < -btol])
  if (V_lo >= V_up) stop('Lower bound is larger than the upper bound!')
  # testing whether eta^t mu == 0
  require(truncnorm)
  if (type == 'test'){
    temp1 <- ptruncnorm(eta %*% y, V_lo, V_up, 0, sqrt(sig2 * sum(eta^2)))
    Test.stat <- 2 * min(temp1, 1 - temp1)
    return(Test.stat)
  } else {
    mu_tnorm <- function(mu) ptruncnorm(eta %*% y, V_lo, V_up, mu, sqrt(sig2 * sum(eta^2)))
    U.search <- V_up # min(V_up, eta %*% y + 5 * sqrt(sig2))
    L.search <- V_lo # min(V_lo, eta %*% y - 3 * sqrt(sig2 * sum(eta^2)))
    CI_up <- uniroot(function(x) mu_tnorm(x) - lsig/2, c(L.search,U.search))$root
    CI_lo <- uniroot(function(x) mu_tnorm(x) - (1-lsig/2), c(L.search,U.search))$root
    return(c(CI_lo, CI_up))
  }
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
  r <- sum(diag(P))
  R <- P %*% y
  R.norm <- sqrt(sum(R^2))
  u <- R/R.norm
  v <- y - R
  Test.stat <- R.norm/sqrt(sig2)
  
  alp <- Gama %*% u * sqrt(sig2)
  temp <- (dd - Gama %*% v)/alp
  V_up <- min(temp[alp > btol])
  V_lo <- max(c(temp[alp < -btol],0))
  if (Test.stat > V_up | Test.stat < V_lo) return('Impossible!')
  p_val <- 1 - diff(pchisq(c(Test.stat^2,V_up^2), df = r))/diff(pchisq(c(V_lo^2,V_up^2), df = r))
  return(p_val)
}


#####

# proj.mat calculates the projection matrix onto the columns of a matrix
proj.mat <- function(A, rtol = 1e-7){
  svd.A <- svd(A)
  d <- svd.A$d
  r <- sum(d>rtol)
  u <- svd.A$u[,1:r,drop = F]
  return(u %*% t(u))
}





