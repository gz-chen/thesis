setwd('/Users/gzchen/Documents/GitHub/thesis/10_2')
rm(list = ls())
load('Simu10_2')

pdf('new.pdf')
par(mfrow = c(3,3))
for (i in 1:21) {
  beta_esti <- STOR[[i]]$beta_esti
  matplot(cbind(beta_esti,beta_true), col = c('black','red'), pch = c(1,2),
          ylab = 'beta',xlab = paste0('alpha=',(i-1)*0.05), main = 'The coefficients')
  legend(30, -0.1, c('Est.','True'), pch = c(1,2), col = c('black','red'))
}
dev.off()
