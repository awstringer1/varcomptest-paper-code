### Analyze the oxide bootstrap results ###

library(tidyverse)
basepath <- getwd()
resultspath <- basepath
figurepath <- basepath

results <- readr::read_csv(file.path(resultspath, "oxideanalysis.csv"))

## Histograms of the MLE ##

# Values copied from Julia
tau1mle <- 10.3352
tau2mle <- 2.8534
# LRT
# julia> -vmod.opt.val + vmod.optcond.val
lrt <- 2.8877796463152805

HISTBINS <- 30
PLOTTEXTSIZE <- 2
PLOTLINEWIDTH <- 2
PLOTMARGIN <- 5
PLOTSIZE <- 7

## Tau 1 ##
pdf(file = file.path(figurepath, "oxide-tau1.pdf"), height = PLOTSIZE, width = PLOTSIZE)
par(mar = rep(PLOTMARGIN, 4))
hist(
  results$Lot, 
  freq = FALSE, 
  breaks = HISTBINS,
  main = expression("Bootstrap samples of"~widehat(tau)[1]~"under"~H[0]),
  xlab = expression(widehat(tau)[1]),
  ylab = "Bootstrap density",
  cex.main = PLOTTEXTSIZE,
  cex.axis = PLOTTEXTSIZE,
  cex.lab = PLOTTEXTSIZE
)
abline(v = tau1mle, lty = "dashed", lwd = PLOTLINEWIDTH)
dev.off()

## Tau 2 ##
pdf(file = file.path(figurepath, "oxide-tau2.pdf"), height = PLOTSIZE, width = PLOTSIZE)
par(mar = rep(PLOTMARGIN, 4))
hist(
  results$WaferLot, 
  freq = FALSE, 
  breaks = HISTBINS,
  main = expression("Bootstrap samples of"~widehat(tau)[2]~"under"~H[0]),
  xlab = expression(widehat(tau)[2]),
  ylab = "Bootstrap density",
  cex.main = PLOTTEXTSIZE,
  cex.axis = PLOTTEXTSIZE,
  cex.lab = PLOTTEXTSIZE
)
abline(v = tau2mle, lty = "dashed", lwd = PLOTLINEWIDTH)
dev.off()

## Tau1 - Tau2 ##

pdf(file = file.path(figurepath, "oxide-taudiff.pdf"), height = PLOTSIZE, width = PLOTSIZE)
par(mar = rep(PLOTMARGIN, 4))
hist(
  results$Lot - results$WaferLot, 
  freq = FALSE, 
  breaks = HISTBINS,
  main = expression("Bootstrap samples of"~widehat(tau)[1] - widehat(tau)[2]~"under"~H[0]),
  xlab = expression(widehat(tau)[1] - widehat(tau)[2]),
  ylab = "Bootstrap density",
  cex.main = PLOTTEXTSIZE,
  cex.axis = PLOTTEXTSIZE,
  cex.lab = PLOTTEXTSIZE
)
abline(v = tau1mle - tau2mle, lty = "dashed", lwd = PLOTLINEWIDTH)
dev.off()

## LRT ##

pdf(file = file.path(figurepath, "oxide-lrt.pdf"), height = PLOTSIZE, width = PLOTSIZE)
par(mar = rep(PLOTMARGIN, 4))
hist(
  results$lrt, 
  freq = FALSE, 
  breaks = HISTBINS,
  main = expression("Bootstrap samples of LRT under"~H[0]),
  xlab = "LRT",
  ylab = "Bootstrap density",
  cex.main = PLOTTEXTSIZE,
  cex.axis = PLOTTEXTSIZE,
  cex.lab = PLOTTEXTSIZE
)
abline(v = lrt, lty = "dashed", lwd = PLOTLINEWIDTH)
dev.off()

