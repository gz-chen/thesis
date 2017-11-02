# this.dir <- dirname(parent.frame(2)$ofile)
# setwd(this.dir)
setwd('/Users/gzchen/Documents/GitHub/thesis')
source('MyFuns.R')


RES <- NULL
P3 <- c(0,0.01,0.03)

for (ip in 1:length(P3)){
  {
    set.seed(1012)
    p <- 50
    phy <- rtree(p)
    phy$tip.label <- 1:p
    
    ######Define clusters & plot tree########
    
    #define the associated clusters by the highest internal node
    CLS <- c(56,75)
    
    asso_tips <- plot_tree(phy, cls = CLS)
    
    ############################################
    ################generate simu data##########
    
    ln_par <- ln_param(phy)
    # generate parameters used for generating OTUs
    
    
    gma <- gen_gma(ln_par = ln_par, p1 = 0.02, p2 = 0.05, tree = phy, cls = CLS, p3 = P3[ip])
    # generate coefficients with given effect size
    
    
    beta_true <- true_beta(phy, CLS, gma)
    # recover full vector
    
    
    #######generate penalty matrix & Data prep#######
    
    # generate penalty matrix
    DW <- gen_D(phy, m = 2, weight = 'max', type = 'myown')
    # data prep. for pen. reg.
    ref <- 1
  }
  
  {
    ################check selective type I error#####
    NN <- 100
    correct <- 0
    chi_s<- NULL
    bounds <- NULL
  }
  
  # set.seed(as.numeric(Sys.time()))
  set.seed(15)
  for (i in 1:NN) {
    {
      # data
      Data <- gen_dat(n = 500, ln_par = ln_par, gma = gma, tree = phy, cls = CLS, sig = 1)
      model.data <- data_prep(Data, DW, ref, alp = 0.26, normlz = F)
      
      {
        G_cen <- model.data$G_cen
        y_cen <- model.data$y_cen
        X_cen1 <- model.data$X_cen1
        D1 <- model.data$D1
      }
      
      # variabl selection
      res <- select_OTU(y_cen, X_cen1, D1)
      
      # model result and assessment
      beta_esti <- esti_beta(res$beta[,res$stop.index], ref)
      
      #plot beta
      plot_beta_bic(beta_true, beta_esti, res$bic)
      
      fuse_ass <- assess_fuse(phy, beta_esti, beta_true)
      sparse_ass <- assess_sparse(beta_esti, beta_true)
    }
    
    
    
    if (!fuse_ass$nFPR & !sparse_ass$FPR) {
      # form new covariates
      new_res <- form_new(model.data$X_cen, beta_esti)
      X_new <- new_res$X_new
      X_new1 <- cbind(G_cen, X_new)
      
      # new testing scheme for interaction
      norms <- get_dir(y_cen, X_new, G_cen, type = 'interact')
      
      lr <- sqrt(norms$deg * res$sig2) 
      V_up <- get_bound(y_cen, norms$eta, sd = lr, max = 10, upper = T, verbose = T)
      V_lo <- get_bound(y_cen, norms$eta, sd = lr, max = 10, upper = F, verbose = T)
      
      std <- sqrt(res$sig2)
      org <- norm.v(norms$eta)
      V_lo <- max(org + V_lo, 0) / std # if chi-test
      V_up <- (org + V_up) / std
      bnd <- c(V_lo, V_up)
      # truncation test:
      p_val2 <- trunc_test(test = norms$test, bnd = bnd, 
                           type = 'chi', deg = norms$deg, side = 'one', sig2 = res$sig2)
      bounds <- cbind(bounds, bnd)
      # oracle test use true model:
      p_val1 <- oracle_test(model.data)
      # naive test without truncation:
      p_val3 <- naive_test(norms$test, deg = norms$deg, sig2 = res$sig2)
      # fool test without model selection:
      p_val4 <- fool_test(model.data)
      # store the results:
      chi_s <- cbind(chi_s, c(p_val1, p_val2, p_val3, p_val4))
      correct <- correct + 1
    }
    cat('i=',i,'in ip=',ip,'done!\n')
  }
  RES[[ip]] = list(p_values = chi_s, bounds = bounds, corrects = correct)
}


# save.image('type1err.RData')






