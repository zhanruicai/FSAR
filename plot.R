rm(list = ls())
library(latex2exp)

N = 100
Num.Observed = 50
filename = paste("Result_SBM", N, "-", Num.Observed, ".rda", sep = "")
load(filename)
rhobeta <- all_result$rhobeta
CST <- all_result$CST
all_coef <- rhobeta[c(1:100)]
esti_mean <- apply(simplify2array(all_coef), 1:2, mean)
sd_coef = apply(simplify2array(all_coef), 1:2, sd)
up_line = esti_mean+1.959964*sd_coef
lo_line = esti_mean-1.959964*sd_coef
timeline = seq(from = 1.5, to = 8.5, length.out = 25)

par(mfrow = c(1,3))
plot(timeline, esti_mean[,1], type = "l", ylim = c(0,0.5), main = TeX('$\\widehat{\\rho(t)}$'), ylab = "", xlab = "Time")
lines(up_line[,1], type="l", lty=2)
lines(lo_line[,1], type="l", lty=2)
plot(timeline, esti_mean[,2], type = "l", ylim = c(-0.3,0.8), main = TeX('$\\widehat{\\beta_1(t)}$'), ylab = "", xlab = "Time")
lines(timeline, up_line[,2], type="l", lty=2)
lines(timeline, lo_line[,2], type="l", lty=2)
plot(timeline, esti_mean[,3], type = "l", ylim = c(0.2,0.9), main = TeX('$\\widehat{\\beta_2(t)}$'), ylab = "", xlab = "Time")
lines(timeline, up_line[,3], type="l", lty=2)
lines(timeline, lo_line[,3], type="l", lty=2)
