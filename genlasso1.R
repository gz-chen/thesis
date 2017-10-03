
svdsolve <- function(A,b,rtol) {
  s = svd(A)
  di = s$d
  ii = di>rtol
  di[ii] = 1/di[ii]
  di[!ii] = 0
  return(list(x=s$v %*% (di* (t(s$u) %*% b)),q=sum(ii)))
}
  
  X <- X_cen1; y <- y_cen
  x <- svd(X)
  # if (min(x$d) < rtol) stop('X does not have full col. rank!')
  y2 <- as.numeric(x$u %*% (t(x$u) %*% y))
  X_inv <- x$v %*% (t(x$u)/x$d)
  D2 <- D1 %*% X_inv
   
  y = y2; D = D2; approx=FALSE; maxsteps=2000; minlam=0;
  rtol=1e-7; btol=1e-7; verbose=FALSE; object=NULL;
  
  # If we are starting a new path
  if (is.null(object)) {
    m = nrow(D)
    n = ncol(D)
    
    # Compute the dual solution at infinity, and
    # find the first critical point
    sv = svdsolve(t(D),y,rtol)    # SVD solver
    uhat = sv$x                   # Dual solution
    q = sv$q                      # Rank of D
    ihit = which.max(abs(uhat))   # Hitting coordinate
    hit = abs(uhat[ihit])         # Critical lambda
    s = sign(uhat[ihit])          # sign
    
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
    
    lams[1] = hit
    h[1] = TRUE
    df[1] = n-q
    u[,1] = uhat
    
    # Other things to keep track of, but not return
    r = 1                      # Size of boundary set
    B = ihit                   # Boundary set
    I = seq(1,m)[-ihit]        # Interior set
    Ds = D[ihit,]*s            # Vector t(D[B,])%*%s
    DD_1 = D[-ihit,,drop=FALSE]  # Matrix D[I,]
    DD_2 = D[ihit,,drop=FALSE]   # Matrix D[B,]
    k = 2                      # What step are we at?
  } else { 
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
  

    while (k<=maxsteps && lams[k-1]>=minlam) {
      ##########
      # Check if we've reached the end of the buffer
      if (k > length(lams)) {
        buf = length(lams)
        lams = c(lams,numeric(buf))
        h = c(h,logical(buf))
        df = c(df,numeric(buf))
        u = cbind(u,matrix(0,m,buf))
      }
      
      ##########
      # If the interior is empty, then nothing will hit
      if (r==m) {
        aa = bb = numeric(0)
        hit = 0
        q = 0
      } else {
        sv = svdsolve(t(DD_1),cbind(y,Ds),rtol)
        aa = as.numeric(sv$x[,1])
        bb = as.numeric(sv$x[,2])
        q = sv$q
        shit_s = sign(aa)
        hit_s = aa/(bb+shit_s);
        
        # Make sure none of the hitting times are larger
        # than the current lambda (precision issue)
        hit_s[hit_s>lams[k-1]+btol] = 0
        hit_s[hit_s>lams[k-1]] = lams[k-1]
        
        ihit = which.max(hit_s)
        hit = hit_s[ihit]
        shit = shit_s[ihit]
      }
      
      ##########
      # If nothing is on the boundary, then nothing will leave
      # Also, skip this if we are in "approx" mode
      if (r==0 || approx) {
        leave = 0
      } else {
        cc = s*(DD_2%*%(y-t(DD_1)%*%aa))
        dd = s*(DD_2%*%(Ds-t(DD_1)%*%bb))
        leave_s = cc/dd
        
        # c must be negative
        leave_s[cc >= 0] = 0
        
        # Make sure none of the leaving times are larger
        # than the current lambda (precision issue)
        leave_s[leave_s>lams[k-1]+btol] = 0
        leave_s[leave_s>lams[k-1]] = lams[k-1]
        
        ileave = which.max(leave_s)
        leave = leave_s[ileave]
      }
      
      ##########
      # Stop if the next critical point is negative
      if (hit<=0 && leave<=0) break
      
      # If a hitting time comes next
      if (hit > leave) {
        # Record the critical lambda and solution
        lams[k] = hit
        h[k] = TRUE
        df[k] = n-q
        uhat = numeric(m)
        uhat[B] = hit*s
        uhat[I] = aa-hit*bb
        u[,k] = uhat
        
        # Update all of the variables
        r = r+1
        B = c(B,I[ihit])
        I = I[-ihit]
        Ds = Ds + DD_1[ihit,]*shit
        s = c(s,shit)
        DD_2 = rbind(DD_2,DD_1[ihit,])
        DD_1 = DD_1[-ihit,,drop=FALSE]
        
        if (verbose) {
          cat(sprintf("\n%i. lambda=%.3f, adding coordinate %i, |B|=%i...",
                      k,hit,B[r],r))
        }
      } else {
        # Record the critical lambda and solution
        lams[k] = leave
        h[k] = FALSE
        df[k] = n-q
        uhat = numeric(m)
        uhat[B] = leave*s
        uhat[I] = aa-leave*bb
        u[,k] = uhat
        
        # Update all of the variables
        r = r-1
        I = c(I,B[ileave])
        B = B[-ileave]
        Ds = Ds - DD_2[ileave,]*s[ileave]
        s = s[-ileave]
        DD_1 = rbind(DD_1,DD_2[ileave,])
        DD_2 = DD_2[-ileave,,drop=FALSE]
        
        
        if (verbose) {
          cat(sprintf("\n%i. lambda=%.3f, deleting coordinate %i, |B|=%i...",
                      k,leave,I[m-r],r))
        }
      }
      if (k == 200) break
      # Step counter
      k = k+1
    }
# origin version

