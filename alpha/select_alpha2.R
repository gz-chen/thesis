this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
#setwd('/Users/gzchen/Documents/GitHub/thesis/alpha')
source('../MyFuns.R')
##########################
######generate tree#######

set.seed(1012)
p <- 50
phy <- rtree(p)
phy$tip.label <- 1:p

#define the associated clusters by the highest internal node
CLS <- c(56,75)

asso_tips <- plot_tree(phy, cls = CLS)

# parameters
ln_par <- ln_param(phy)
# generate parameters used for generating OTUs


gma <- gen_gma(ln_par = ln_par, p1 = 0.02, p2 = 0.05, tree = phy, cls = CLS)
# generate coefficients with given effect size


beta_true <- true_beta(phy, CLS, gma)
# recover full vector


# generate penalty matrix
# DW <- gen_D(phy, m = 2, weight = 'max', type = 'myown')
# data prep. for pen. reg.
ref <- 1



alp.values <- seq(0, 1, by = 0.05)
NN <- 5
STOR <- list()
pp <- length(alp.values)
mtd <- c('myown','wang1')

for (i in 1:NN){
  temp <- matrix(nrow = 5, ncol = pp)
  Data <- gen_dat(n = 300, ln_par = ln_par, gma = gma, tree = phy, cls = CLS, sig = 1)
  STOR[[i]] <- list()
  for (j in 1:2){
    DW <- gen_D(phy, m = 2, weight = 'max', type = mtd[j])
    for (ii in 1:pp) {
      model.data <- data_prep(Data, DW, ref, alp = alp.values[ii], normlz = F)
      G_cen <- model.data$G_cen
      y_cen <- model.data$y_cen
      X_cen1 <- model.data$X_cen1
      D1 <- model.data$D1
      
      res2 <- gen_select2(y_cen, X_cen1, D1, btol = 1e-6)
      beta_esti <- esti_beta(res2$beta[,res2$stop.index], ref)
      plot_beta_bic(beta_true, beta_esti, res2$bic_n)
      # plot_beta_bic(beta_true, esti_beta(res2$beta[,160], ref), res2$bic_n)
      
      fuse_ass <- assess_fuse(phy, beta_esti, beta_true)
      sparse_ass <- assess_sparse(beta_esti, beta_true)
      temp[,ii] <- c(fuse_ass$nFPR, fuse_ass$nFNR, sparse_ass$FPR, sparse_ass$FNR, res2$bic_n[res2$stop.index])
      print(paste0("Stop.index = ", res2$stop.index))
      print(paste0("BIC = ", res2$bic_n[res2$stop.index]))
      print(paste0('i=',i,'; type=',j,'; ii=',ii,' done!'))
    }
    STOR[[i]][[j]] <- temp
  }
}



####################################################
## some checks for fun
# Data <- gen_dat(n = 5000, ln_par = ln_par, gma = gma, tree = phy, cls = CLS, sig = 1)
# lm1 <- lm(Data$y ~ Data$G + Data$Z)
# gma #true parameters
# summary(lm1) #run lm to see whether it yields the coef and std. error and R.squared desired
# 
# lm2 <- lm(Data$y ~ Data$G + Data$X[,-1])
# plot(lm2$coefficients) # no sparsity or fusion pattern shown
# summary(lm2)
# 
# lm3 <- lm(Data$y ~ Data$G + Data$X[,c(4:11,24:35)])
# summary(lm3)
# plot(lm3$coefficients)
