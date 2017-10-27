this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
#setwd('/Users/gzchen/Documents/GitHub/thesis')
source('MyFuns.R')
##########################
######generate tree#######

set.seed(1012)
p <- 50
phy <- rtree(p)
phy$tip.label <- 1:p

######Define clusters & plot tree########

#define the associated clusters by the highest internal node
CLS <- c(56,75)

asso_tips <- plot_tree(phy, cls = CLS)

############################################
################generate data###############

ln_par <- ln_param(phy)
# generate parameters used for generating OTUs


gma <- gen_gma(ln_par = ln_par, p1 = 0.02, p2 = 0.05, tree = phy, cls = CLS)
# generate coefficients with given effect size


beta_true <- true_beta(phy, CLS, gma)
# recover full vector



#######generate penalty matrix & Data prep#######

# generate penalty matrix
DW <- gen_D(phy, m = 2, weight = 'max', type = 'wang1')
# data prep. for pen. reg.
ref <- 1

################check selective type I error#####

NN <- 300
correct <- 0
# P_val.norm <- matrix(nrow = NN, ncol = ncol(X_new1))
P_val.chi1 <- NULL
P_val.chi2 <- NULL
P_val.norm.zero1 <- NULL
P_val.norm.zero2 <- NULL
P_val.norm.nonzero1 <- NULL
P_val.norm.nonzero2 <- NULL


for (i in 1:NN) {
  Data <- gen_dat(n = 500, ln_par = ln_par, gma = gma, tree = phy, cls = CLS, sig = 1)
  model.data <- data_prep(Data, DW, ref, alp = 0.26, normlz = F)
  
  G_cen <- model.data$G_cen
  y_cen <- model.data$y_cen
  X_cen1 <- model.data$X_cen1
  D1 <- model.data$D1
  # model
  res <- gen_select(y_cen, X_cen1, D1, btol = 1e-6)
  # model result and assessment
  beta_esti <- esti_beta(res$beta[,res$stop.index], ref) 
  plot_beta_bic(beta_true, beta_esti, res$bic)
  fuse_ass <- assess_fuse(phy, beta_esti, beta_true)
  sparse_ass <- assess_sparse(beta_esti, beta_true)
  
  if (!fuse_ass$nFPR & !sparse_ass$FPR) {
    rdig <- 6
    approx_beta <- round(beta_esti, rdig)
    # different values of beta exist
    level_beta <- unique(approx_beta)
    # how many OTUs correspond to the values in level_beta
    # num_ele_level <- sapply(level_beta, FUN = function(x) sum(approx_beta == x))
    ref_new <- which(level_beta == approx_beta[ref])
    # form new covariates
    X_new <- sapply(level_beta, FUN = function(x) rowSums(model.data$X_cen[,approx_beta == x,drop = F]))
    X_new1 <- cbind(G_cen, X_new[,-ref_new])
    # assumed that ref_new == 1
    # H_0 is false
    beta_nonzero <- beta_true != 0
    beta_nonzero_level <- unique(approx_beta[beta_nonzero])
    index_nonzero <- sapply(beta_nonzero_level, function(x) which(level_beta == x))
    # H_0 is true
    beta_zero <- !beta_nonzero
    beta_zero_level <- unique(approx_beta[beta_zero])
    index_zero <- setdiff(sapply(beta_zero_level, function(x) which(level_beta == x)), ref_new)
    
    for (j in 2:ncol(X_new1)){
      eta <- ginv(X_new1)[j,]
      if (j %in% index_zero) {
        P_val.norm.zero1 <- c(P_val.norm.zero1, test_norm(y_cen, res$gama, res$dd, eta, res$sig2))
        P_val.norm.zero2 <- c(P_val.norm.zero2, test_norm2(y_cen, eta, res$sig2))
      } else {
        P_val.norm.nonzero1 <- c(P_val.norm.nonzero1, test_norm(y_cen, res$gama, res$dd, eta, res$sig2))
        P_val.norm.nonzero2 <- c(P_val.norm.nonzero2, test_norm2(y_cen, eta, res$sig2))
        # if (length(P_val.norm.nonzero1) != length(P_val.norm.nonzero2)) break
      }
    }
    
    # interaction terms
    Intact <- G_cen * X_new[,-ref_new]
    Intact_cen <- t(t(Intact) - colMeans(Intact))
    P <- proj.mat(Intact_cen) %*% (diag(length(y_cen)) - proj.mat(X_new1))

    P_val.chi1 <- c(P_val.chi1, test_chi(y_cen, res$gama, res$dd, P, res$sig2))
    P_val.chi2 <- c(P_val.chi2, test_chi2(y_cen, P, res$sig2))
  }
  print(paste0('i=',i,'done'))
  if (i == 10) break
}


# save.image('type1err.RData')






