this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
#setwd('/Users/gzchen/Documents/GitHub/thesis')
print(this.dir)
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


gma <- gen_gma(ln_par = ln_par, p1 = 0.02, p2 = 0.05, tree = phy, cls = CLS)
# generate coefficients with given effect size


beta_true <- true_beta(phy, CLS, gma)
# recover full vector

####################################################
## some checks for fun
Data <- gen_dat(n = 5000, ln_par = ln_par, gma = gma, tree = phy, cls = CLS, sig = 1)
lm1 <- lm(Data$y ~ Data$G + Data$Z)
gma #true parameters
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

alp.values <- seq(0.2, 0.5, by = 0.02)
NN <- 100
STOR <- list()

for (ii in 1:length(alp.values)) {
  temp <- matrix(nrow = 5, ncol = NN)
  for (i in 1:NN){
    Data <- gen_dat(n = 300, ln_par = ln_par, gma = gma, tree = phy, cls = CLS, sig = 1)
    model.data <- data_prep(Data, DW, ref, alp = alp.values[ii], normlz = F)
    
    G_cen <- model.data$G_cen
    y_cen <- model.data$y_cen
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
    
    res2 <- gen_select2(y_cen, X_cen1, D1, btol = 1e-6)
    # res2$stop.index
    
    # summary of genlasso result
    beta_esti <- esti_beta(res2$beta[,res2$stop.index], ref)
    # plot_beta_bic(beta_true, esti_beta(res2$beta[,160], ref), res2$bic)
    # plot_beta_bic(beta_true, beta_esti, res2$bic)
    fuse_ass <- assess_fuse(phy, beta_esti, beta_true)
    sparse_ass <- assess_sparse(beta_esti, beta_true)
    temp[,i] <- c(fuse_ass$nFPR, fuse_ass$nFNR, sparse_ass$FPR, sparse_ass$FNR, res2$bic[res2$stop.index])
    print(paste0('i=',i,' done'))
  }
  STOR[[ii]] <- temp
  print(paste0('ii=',ii,' done'))
}
