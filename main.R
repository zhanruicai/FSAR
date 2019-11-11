rm(list = ls())
#setwd("/Users/zhanruicai/Library/Mobile Documents/com~apple~CloudDocs/Research/FDA_Network_code/s0")
library(Matrix)
library(poweRlaw)
set.seed(1234)
source("FSAR_functions.R")

N = 100
Num.Observed = 50

num.points = 25
Nrep = 100
rhobeta = rep(list(matrix(0, nrow = num.points, ncol = 2)), length(Nrep))
CST = rep(list(matrix(0, num.points, num.points), length(Nrep)))

for (i in 1:Nrep) {
  cat(i, "\r")
  ## Generate data
  dat = get.func.network(n = N, min.points = Num.Observed - 10, max.points = Num.Observed, net.type = "SBM")
  ## Estimate the coefficients on the entire timeline
  ans = est_entire_timeline(Y = dat$Y, X = dat$v1,Tpoints = dat$Tpoints, W = dat$W, num.points)
  rhobeta[[i]] = ans$est_coef
  cov1 = est_covariance(ans$est_coef, ans$timeline, Y = dat$Y, X = dat$v1,Tpoints = dat$Tpoints, W = dat$W)$cov
  CST[[i]] = cov1
}

## True parameters setting
true_coef = cbind((15-ans$timeline)/40, (ans$timeline-5)^2/25, ans$timeline^0.5/4)

esti_coef = apply(simplify2array(rhobeta), 1:2, median)

MSE_coef =colMeans((true_coef-esti_coef)^2)
print(c('rho', 'beta1','beta2'))
r1 = sqrt(MSE_coef)
print(specify_decimal(r1, 4))


rhobeta$esti_coef = esti_coef
rhobeta$true_coef = true_coef
rhobeta$MSE_coef = MSE_coef
rhobeta$RMSE_coef = sqrt(MSE_coef)

### Estimate the covariance matrix
true_cov = getCovCST(ans$timeline)
esti_cov = apply(simplify2array(CST), 1:2, median)
EUC_dis = norm(true_cov-esti_cov, "F")
print('Euc Norm')
print(EUC_dis)
CST$true_cov = true_cov
CST$esti_cov = esti_cov
CST$EUC_dis = EUC_dis

all_result <- list(CST= CST, rhobeta = rhobeta)
filename = paste("Result_SBM", N, "-", Num.Observed, ".rda", sep = "")
save(all_result, file = filename)





