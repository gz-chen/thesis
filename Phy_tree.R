setwd('/Users/gzchen/Documents/GitHub/thesis')
source('MyFuns.R')

##########################
######generate tree#######

set.seed(1012)
p <- 50
phy <- rtree(p)
phy$tip.label <- 1:p



#########################################
######Define clusters & plot tree########

#define the associated clusters by the highest internal node
CLS <- c(56,75)

asso_tips <- plot_tree(phy, cls = CLS)

############################################
################generate data###############

ln_par <- ln_param(phy)
# generate parameters used for generating OTUs


gma <- gen_gma(ln_par = ln_par, p1 = 0.02, p2 = 0.05, tree = phy, cls = CLS, p3 = 0.03)
# generate coefficients with given effect size


beta_true <- true_beta(phy, CLS, gma)
# recover full vector

####################################################
## some checks for fun
Data <- gen_dat(n = 500, ln_par = ln_par, gma = gma, tree = phy, cls = CLS, sig = 1)
lm1 <- lm(Data$y ~ Data$G + Data$Z)
# gma #true parameters
summary(lm1) #run lm to see whether it yields the coef and std. error and R.squared desired
# 
# lm2 <- lm(Data$y ~ Data$G + Data$X[,-1])
# plot(lm2$coefficients) # no sparsity or fusion pattern shown
# summary(lm2)
# 
# lm3 <- lm(Data$y ~ Data$G + Data$X[,c(4:11,24:35)])
# summary(lm3)
# plot(lm3$coefficients)



#################################################
#######generate penalty matrix & Data prep#######


# generate penalty matrix
DW <- gen_D(phy, m = 2, weight = 'max', type = 'myown')
# data prep. for pen. reg.
ref <- 1

##

#generate data

# alp.values <- seq(0, 1,by = 0.05)
# NN <- 10
# STOR <- list()
# 
# for (ii in 1:length(alp.values)) {
#   temp <- matrix(nrow = 4, ncol = NN)
#   for (i in 1:NN){
    Data <- gen_dat(n = 500, ln_par = ln_par, gma = gma, tree = phy, cls = CLS, sig = 1)
    model.data <- data_prep(Data, DW, ref, alp = 0.26, normlz = F)
    
    G_cen <- model.data$G_cen
    y_cen <- c(model.data$y_cen)
    X_cen1 <- model.data$X_cen1
    D1 <- model.data$D1
    
    ##############genlasso#################
    # ## package: genlasso
    # res1 <- select_OTU(y_cen, X_cen1, D1, svd = T)
    # # res1$g_lasso$hit[1:10]
    # res1$stop_index
    # 
    # plot_beta_bic(beta_true, esti_beta(res1$g_lasso$beta[,res1$stop_index], ref), res1$bic)
    
    # plot_beta_bic(beta_true, esti_beta(res1$g_lasso$beta[,114], ref), res1$bic)
    #######################################
    
    ## my own function
    
    res2 <- gen_select2(y_cen, X_cen1, D1, btol = 1e-6,maxsteps = 300)
    # res2$stop.index
    
    # summary of genlasso result
    beta_esti <- esti_beta(res2$beta[,res2$stop.index], ref)
    # plot_beta_bic(beta_true, esti_beta(res2$beta[,20], ref), res2$bic)
    plot_beta_bic(beta_true, beta_esti, res2$bic_n)
    
    fuse_ass <- assess_fuse(phy, beta_esti, beta_true)
    sparse_ass <- assess_sparse(beta_esti, beta_true)
    temp[,i] <- c(fuse_ass$nFPR, fuse_ass$nFNR, sparse_ass$FPR, sparse_ass$FNR)
    print(paste0('i=',i,' done'))
#   }
#   STOR[[ii]] <- temp
#   print(paste0('ii=',ii,' done'))
# }

    
    
    
    

##########################################################################
{
  rdig <- 6
  approx_beta <- round(beta_esti, rdig)
  level_beta <- unique(approx_beta)
  num_ele_level <- sapply(level_beta, FUN = function(x) sum(approx_beta == x))
  ref_new <- which(level_beta == approx_beta[ref])
  
  X_new <- sapply(level_beta, FUN = function(x) rowSums(model.data$X_cen[,approx_beta == x,drop = F]))
  X_new1 <- cbind(G_cen, X_new[,-ref_new])
  
  # lm4 <- lm(y_cen ~ G_cen + X_new[,-ref_new] + G_cen * X_new[,-ref_new])
  # anova(lm4)
  # summary(lm4)
  # 
  # lm5 <- lm(y_cen ~ G_cen + X_new[,-ref_new])
  # summary(lm5)
  
  p_values <- numeric(ncol(X_new1))
  # CIs <- matrix(nrow = 2, ncol = ncol(X_new1))
  for (j in 1: ncol(X_new1)){
    eta <- ginv(X_new1)[j,]
    p_values[j] <- test_norm(y_cen, res2$gama, res2$dd, eta, res2$sig2)
    # CIs[,j] <- test_norm(y_cen, res2$Gama, res2$dd, eta, res2$sig2, type = 'CI')
  }
  # level_beta
  # p_values
  # CIs
}

y <- y_cen; Gama <- res2$gama; dd <- res2$dd; sig2 <- res2$sig2; btol = 1e-7

Intact <- G_cen * X_new[,-ref_new]
Intact_cen <- t(t(Intact) - colMeans(Intact))
P <- proj.mat(Intact_cen) %*% (diag(length(y_cen)) - proj.mat(X_new1))

test_chi(y_cen, res2$gama, res2$dd, P, res2$sig2)


#########################check selective type I error#####

NN <- 300
correct <- 0
P_val.norm <- matrix(nrow = NN, ncol = ncol(X_new1))
P_val.chi <- numeric(NN)

for (i in 1:NN) {
  Data <- gen_dat(n = 300, ln_par = ln_par, gma = gma, tree = phy, cls = CLS, sig = 1)
  model.data <- data_prep(Data, DW, ref, alp = 0.3, normlz = F)
  
  G_cen <- model.data$G_cen
  y_cen <- model.data$y_cen
  X_cen1 <- model.data$X_cen1
  D1 <- model.data$D1
  
  res <- gen_select(y_cen, X_cen1, D1, btol = 1e-6)
  
  beta_esti <- esti_beta(res$beta[,res$stop.index], ref)
  fuse_ass <- assess_fuse(phy, beta_esti, beta_true)
  sparse_ass <- assess_sparse(beta_esti, beta_true)
  
  if (!fuse_ass$nFPR & !sparse_ass$FPR & !fuse_ass$nFPR & !sparse_ass$FPR) {
    correct <- correct + 1
    for (j in 1:ncol(X_new1)){
      eta <- ginv(X_new1)[j,]
      P_val.norm[i,j] <- test_norm(y_cen, res$gama, res$dd, eta, res$sig2)
    }
    
    P_val.chi[i] <- test_chi(y_cen, res$gama, res$dd, P, res$sig2)
  } else {
    P_val.chi[i] <- NA
    P_val.norm[i,] <- NA
  }
}




# plot_beta_bic(beta_true, esti_beta(res2$beta[,100], ref), res2$bic)

############################################
######### model assessment #################

assess_fuse(phy, beta_esti, beta_true)

assess_sparse(beta_esti, beta_true)


###########################################
########## model reference step ###########

# reform new design matrix and reference level








