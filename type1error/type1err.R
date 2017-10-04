setwd('/Users/gzchen/Documents/GitHub/thesis/type1error')
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
DW <- gen_D(phy, m = 2, weight = 'max', type = 'myown')
# data prep. for pen. reg.
ref <- 1

################check selective type I error#####

NN <- 300
correct <- 0
# P_val.norm <- matrix(nrow = NN, ncol = ncol(X_new1))
P_val.chi <- numeric(NN)

for (i in 1:NN) {
  Data <- gen_dat(n = 300, ln_par = ln_par, gma = gma, tree = phy, cls = CLS, sig = 1)
  model.data <- data_prep(Data, DW, ref, alp = 0.3, normlz = F)
  
  G_cen <- model.data$G_cen
  y_cen <- model.data$y_cen
  X_cen1 <- model.data$X_cen1
  D1 <- model.data$D1
  # model
  res <- gen_select(y_cen, X_cen1, D1, btol = 1e-6)
  # model result and assessment
  beta_esti <- esti_beta(res$beta[,res$stop.index], ref)
  fuse_ass <- assess_fuse(phy, beta_esti, beta_true)
  sparse_ass <- assess_sparse(beta_esti, beta_true)
  # form new covariates
  if (!fuse_ass$nFPR & !sparse_ass$FPR) {
    rdig <- 6
    approx_beta <- round(beta_esti, rdig)
    level_beta <- unique(approx_beta)
    num_ele_level <- sapply(level_beta, FUN = function(x) sum(approx_beta == x))
    ref_new <- which(level_beta == approx_beta[ref])
    
    X_new <- sapply(level_beta, FUN = function(x) rowSums(model.data$X_cen[,approx_beta == x,drop = F]))
    X_new1 <- cbind(G_cen, X_new[,-ref_new])
    # interaction terms
    Intact <- G_cen * X_new[,-ref_new]
    Intact_cen <- t(t(Intact) - colMeans(Intact))
    P <- proj.mat(Intact_cen) %*% (diag(length(y_cen)) - proj.mat(X_new1))
    
    correct <- correct + 1
    # for (j in 1:ncol(X_new1)){
    #   eta <- ginv(X_new1)[j,]
    #   P_val.norm[i,j] <- test_norm(y_cen, res$gama, res$dd, eta, res$sig2)
    # }
    P_val.chi[i] <- test_chi(y_cen, res$gama, res$dd, P, res$sig2)
  } else {
    P_val.chi[i] <- NA
    # P_val.norm[i,] <- NA
  }
  if (i == 2) break
}
# plot_beta_bic(beta_true, beta_esti, res$bic)









