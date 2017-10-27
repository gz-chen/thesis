#' gen_select: a function implementing the generalized lasso
#' using SVD, which also generate the affine contraints in the meantime
#' @param y the response (n)
#' @param X the designe matrix (n,p), assumed to have full col rank
#' @param D the penalty matrix (m,p)
#' @param tol the tolerant level for knot value
#' @param rtol singular values below rtol are viewed as 0; default to be sqrt(.Machine$double.eps)
#' @param maxsteps maximum number of steps, default to be 0
#' @param sig2 estimate of error variance; used for bic computing.
#' @return U : matrix of dual variables (m, maxsteps)
#' @return lambda : vector of knot values (maxsteps)
#' @return fit :
#' @return Gama, ddd : The affine constraints: Gama %*% y <= ddd

sig2 <- 1.48915

# some basic parameters
n <- length(y)
m <- nrow(D1)
p <- ncol(D1)

X <- X_cen1
y <- y_cen

# do some transformations
# reformulate the problem with general X as a signal approx. problem
x <- svd(X)

y2 <- as.numeric(x$u %*% (t(x$u) %*% y))

X_inv <- x$v %*% (t(x$u)/x$d)
D2 <- D1 %*% X_inv

# Intialize things to keep track of & return with
U <- matrix(nrow = m, ncol = maxsteps) # optimal dual variables at each step
lambs <- numeric(maxsteps) # knot values at each step
B <- NULL # boundary set
B_c <- seq(1:m) # interior set
S <- NULL # signs of coordinates in B
bic <- numeric(maxsteps) # adjusted bic at each step
df <- numeric(maxsteps) # df of the fit, i.e. nullity(D2[B_c,])
hh <- logical(maxsteps) # whether hit or not


Gama <- NULL; ddd <- NULL

svd_D <- function(A, b, rtol = 1e-7){
  # A function to calculate inv(t(D2[B_c,]))
  # A = t(D2[B_c,])
  # inv(A) %*% b
  # nullity(D2[B_c,]) = n - rank(t(D2[B_c,]))
  svd.A <- svd(A)
  d <- svd.A$d
  r <- sum(d>rtol)
  d <- d[1:r]
  di <- 1/d
  u <- svd.A$u[,1:r]; v <- svd.A$v[,1:r] # D_k <- u %*% d %*% t(v) is the condensed SVD

  # x <- v %*% ((t(u)/d) %*% b)
  x <- v %*% ( di * t(u) %*% b)
  return(list(x = x, r = r))
}

# svdsolve <- function(A,b,rtol = 1e-7) {
#   s = svd(A)
#   di = s$d
#   ii = di>rtol
#   di[ii] = 1/di[ii]
#   di[!ii] = 0
#   return(list(x=s$v %*% (di* (t(s$u) %*% b)),r=sum(ii)))
# }

##First step
temp <- svd_D(t(D2), y2)
U[,1] <- temp$x
i_hit <- which.max(abs(U[,1]))
r_hit <- sign(U[i_hit,1])
lambs[1] <- abs(U[i_hit,1])
hh[1] <- T
df[1] <- n - temp$r
# bic[1] <- sum((temp$Prj %*% y2)^2) + log(n) * sig2 * df[1]
B <- c(B,B_c[i_hit]) # must hit
B_c <- B_c[-i_hit]
S <- c(S,r_hit)
Ds <- D2[i_hit,] * S # vector t(D[B,]) %*% S
D_1 <- D2[-i_hit,,drop = F] # matrix D[B_c,]
D_2 <- D2[i_hit,,drop = F] # matrix D[B,]



k <- 2

while(k <= maxsteps & lambs[k-1] > 0){
  # see below
  temp <- svd_D(t(D_1), cbind(y2,Ds), rtol = 1e-7)

  df[k] <- n - temp$r

  # hitting times
  a = as.numeric(temp$x[,1])
  b = as.numeric(temp$x[,2])
  # a <- as.vector(temp$inv %*% y2)
  # b <- as.vector(temp$inv %*% Ds)

  R <- sign(a)

  hits <- a/(R+b)

  hits[hits > lambs[k-1] + rtol] <- 0
  hits[hits > lambs[k-1]] <- lambs[k-1]

  i_hit <- which.max(hits)
  hit <- hits[i_hit]
  r_hit <- R[i_hit]

  # leaving times
  c <- S * (D_2 %*% (y2 - t(D_1) %*% a))
  # c1 <- S * (D_2 %*% (y2 - t(D_1) %*% aa))
  d <- S * (D_2 %*% (Ds - t(D_1) %*% b))

  leaves <- c/d
  leaves[c >= 0] <- 0

  leaves[leaves > lambs[k-1] + rtol] <- 0
  leaves[leaves > lambs[k-1]] = lambs[k-1]

  i_leave <- which.max(leaves)
  leave <- leaves[i_leave]


  if (hit > leave) {
    lambs[k] <- hit
    hh[k] <- T

    uhat <- numeric(m)
    uhat[B] <- hit*S
    uhat[B_c] <- a - hit*b
    U[,k] <- uhat

    B <- c(B,B_c[i_hit])
    B_c <- B_c[-i_hit]
    S <- c(S, r_hit)

    Ds <- Ds + D_1[i_hit,] * r_hit
    D_2 <- rbind(D_2,D_1[i_hit,])
    D_1 <- D_1[-i_hit,,drop = F]


  } else {
    lambs[k] <- leave
    hh[k] <- F
    uhat <- numeric(m)
    uhat[B] <- leave*S
    uhat[B_c] <- a - leave*b
    U[,k] <- uhat


    B_c <- c(B_c, B[i_leave])
    B <- B[-i_leave]
    Ds <- Ds - D_2[i_leave,] * S[i_leave]
    S <- S[-i_leave]

    D_1 <- rbind(D_1,D_2[i_leave,])
    D_2 <- D_2[-i_leave,,drop = F]

  }
  if (k == 100) break
  k <- k + 1
}
# my version