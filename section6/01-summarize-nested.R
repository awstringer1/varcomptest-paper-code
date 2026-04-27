### Summarize simulation data, make plots ###



PLOTWIDTH <- PLOTHEIGHT <- 7
GGWIDTH <- PLOTWIDTH
GGHEIGHT <- PLOTHEIGHT
PLOTTEXTSIZE <- 1.5
GGTEXTSIZE <- 12
PLOTMARGIN <- 4

args <- commandArgs(TRUE) # Returns character(0) if interactive
if (length(args) > 0) {
  JOBID <- args[1]
} else {
  JOBID <- "1323158"
}

## Packages ##
pkgs <- c(
  "dplyr",
  "readr",
  "ggplot2",
  "dplyr"
)
for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste0("Could not find package ", pkg, ", installing from CRAN.\n"))
    install.packages(pkg)
    require(pkg, character.only = TRUE, quietly = TRUE)
  }
}

## Set paths ##
homedir <- getwd()
basepath <- homedir
resultspath <- file.path(basepath, "results")
if(!dir.exists(resultspath)) dir.create(resultspath)
figurespath <- file.path(basepath, "figures")
if(!dir.exists(figurespath)) dir.create(figurespath)

# Get the file(s) with that job id
files <- list.files(resultspath)
files <- files[grep(JOBID, files)]
if (length(files) > 1) stop(paste0("Directory contains ", length(files), " files with JOBID = ", JOBID))
if (length(files) == 0) stop(paste0("No files found for JOBID = ", JOBID))

cat("Loading simulation results from ", files, "\n", sep = "")

## load p values
pvals <- read_csv(file.path(resultspath, files), progress = FALSE, show_col_types = FALSE)

simname <- sub("^simulationresults-", "", files)
simname <- sub(".csv$", "", simname)

## QQ plot
qqpvals <- function(var, res = pvals, saveplot = FALSE) {
  
  H0 = ""
  H1 = ""
  if (var == "pvalboot") {
    H0 = expression("H"[0]*":"~tau[1]~"="~tau[2]*", H"[1]*":"~tau[1] != tau[2])
  } else if (var == "pvaloneside") {
    H0 = expression("H"[0]*":"~tau[1]~"="~tau[2]*", H"[1]*":"~tau[1] > tau[2])
  } else if (var == "pvalbootzero") {
    H0 = "tau = 0"
    H1 = "tau != 0"
  } else {
    stop("Unknown variable.")
  }


  
  res <- filter(res, tau1 == tau2) %>% 
    mutate(tau = factor(tau1), taudiffest = tau1est - tau2est, tauvalest = (tau1est + tau2est) / 2) %>%
    group_by(m, n, r, tau) %>%
    arrange(.data[[var]], .by_group = TRUE) %>%
    mutate(
      pp = ppoints(n()),
      qtheoretical = qunif(pp),
      qempirical = quantile(.data[[var]], probs = pp)
    )

  if (var == "pvalboot") grid <- seq(0, 1, length.out = 1e04)
  if (var == "pvaloneside") grid <- seq(0, 0.5, length.out = 1e04)
  kstests <- res %>%
    summarize(
      ks_stat = max(abs(ecdf(.data[[var]])(grid) - grid)),
      ks_stat_p05 = max(abs(ecdf(.data[[var]])(grid) - grid)[grid <= .05]),
      ks_idx = which.max(abs(ecdf(.data[[var]])(grid) - grid)),
      ks_val = grid[which.max(abs(ecdf(.data[[var]])(grid) - grid))]
    )
  
  res <- res %>%
    left_join(kstests, by = c("m", "n", "r", "tau"))
  
  qqplt <- ggplot(res, aes(x = qtheoretical, y = qempirical)) +
    theme_bw() +
    facet_grid(m + n ~ r, labeller = label_both) +
    geom_point(aes(colour = tau), alpha = .05, pch = ".") +
    geom_abline(slope = 1, intercept = 0, linewidth = .2) +
    geom_vline(aes(colour = tau, xintercept = ks_val), linetype = "dotted") +
    labs(
      title = "Uniform QQ-plot of bootstrap p-values",
      subtitle = H0,
      x = "Unif(0,1) quantiles", y = "Empirical quantiles",
      colour = "tau1 = tau2"
    ) +
    scale_colour_brewer(type = "seq", palette = "Greys") +
    theme(text = element_text(size = GGTEXTSIZE), legend.position = "bottom") +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))

  ksplt <- ggplot(kstests, aes(x = tau)) +
    theme_bw() +
    facet_grid(m + n ~ r, labeller = label_both) +
    geom_point(aes(y = ks_stat), pch = 16) +
    geom_point(aes(y = ks_stat_p05), pch = 15) +
    labs(
      title = "KS statistics for uniformity of bootstrap p-values",
      subtitle = H0,
      x = expression(tau[1]~"="~tau[2]), y = "Empirical quantiles"
    ) +
    scale_colour_brewer(type = "seq", palette = "Greys") +
    theme(text = element_text(size = GGTEXTSIZE), legend.position = "bottom") +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    coord_cartesian(ylim = c(0, .15))
  
  tauplt <- ggplot(res) +
    theme_bw() + 
    theme(text = element_text(size = GGTEXTSIZE), legend.position = "bottom", axis.text.x = element_text(angle = 45)) +
    facet_grid(m + n ~ r, labeller = label_both) +
    geom_boxplot(aes(x = tau, y = tauvalest), outlier.size=.3) +
    geom_hline(aes(yintercept = tau1), linetype = "dotted") +
    labs(
      title = expression("Estimated common value"~(hat(tau)[1]+hat(tau)[2])/2),
      subtitle = H0,
      x = expression(tau), y = expression((hat(tau)[1]+hat(tau)[2])/2)
    )


  plotnameqq <- paste0(simname, "-qq-", var, ".pdf")
  plotnameks <- paste0(simname, "-ks-", var, ".pdf")
  plotnametau <- paste0(simname, "-tau-", var, ".pdf")
  
  plotpathqq <- file.path(figurespath, plotnameqq)
  plotpathks <- file.path(figurespath, plotnameks)
  plotpathtau <- file.path(figurespath, plotnametau)
  
  if (saveplot) {
    ggsave(file = plotpathqq, plot = qqplt, width = 7, height = 7)
    ggsave(file = plotpathks, plot = ksplt, width = 7, height = 7)
    ggsave(file = plotpathtau, plot = tauplt, width = 7, height = 7)
  } else {
    return(list(qqplt, ksplt, tauplt))
  }
}

vars <- c("pvaloneside", "pvalboot")
for (v in vars) {
  qqpvals(v, saveplot = TRUE)
}

## power plot

powerplot <- function(res = pvals, saveplot = FALSE) {
  alpha <- .05
  # One sided alternative
  res$taudiff <- res$tau1 - res$tau2
  # Two-sided alternative
  res$abstaudiff <- abs(res$tau1 - res$tau2)
  powertwoside <- filter(res, tau1 >= tau2) %>%
    group_by(m, n, r, tau1, tau2, abstaudiff) %>%
    summarize(
      powermean = mean(pvalboot <= alpha),
      powerse = sqrt(powermean * (1 - powermean) / n())
    )
  poweroneside <- filter(res, tau1 >= tau2) %>%
    group_by(m, n, r, tau1, tau2, taudiff) %>%
    summarize(
      powermean = mean(pvalboot <= alpha),
      powerse = sqrt(powermean * (1 - powermean) / n())
    )
  facetlabeller <- labeller(
    m = function(x) paste0("N. Groups = ", x),
    n = function(x) paste0("Group Size = ", x)
  )

  plottwoside <- ggplot(powertwoside, aes(x = abstaudiff)) +
    theme_bw() + 
    theme(text = element_text(size = GGTEXTSIZE), legend.position = "bottom", axis.text.x = element_text(angle = 45)) +
    facet_grid(m + n ~ r, labeller = label_both) +
    geom_line(aes(y = powermean, linetype = factor(tau1))) +
    geom_line(aes(y = powermean + 2 * powerse, group = factor(tau1)), linetype = "solid", alpha = .2) +
    geom_line(aes(y = powermean - 2 * powerse, group = factor(tau1)), linetype = "solid", alpha = .2) +
    labs(x = expression("|"*tau[1]-tau[2]*"|"), y = "Empirical Power", linetype = "Reps") +
    scale_y_continuous(breaks = seq(0, 1, by = .25), labels = scales::percent_format())
    
  plotoneside <- ggplot(poweroneside, aes(x = taudiff)) +
    theme_bw() + 
    theme(text = element_text(size = GGTEXTSIZE), legend.position = "bottom", axis.text.x = element_text(angle = 45)) +
    facet_grid(m + n ~ r, labeller = label_both) +
    geom_line(aes(y = powermean, linetype = factor(tau1))) +
    geom_line(aes(y = powermean + 2 * powerse, group = factor(tau1)), linetype = "solid", alpha = .2) +
    geom_line(aes(y = powermean - 2 * powerse, group = factor(tau1)), linetype = "solid", alpha = .2) +
    labs(x = expression(tau[1]-tau[2]), y = "Empirical Power", linetype = "Reps") +
    scale_y_continuous(breaks = seq(0, 1, by = .25), labels = scales::percent_format())
    

  if(saveplot) {
    plotnametwoside <- paste0(simname, "-power-twoside.pdf")
    plotpathtwoside <- file.path(figurespath, plotnametwoside)
    ggsave(
      file = plotpathtwoside,
      plot = plottwoside,
      width = GGWIDTH, height = GGHEIGHT
    )
    plotnameoneside <- paste0(simname, "-power-oneside.pdf")
    plotpathoneside <- file.path(figurespath, plotnameoneside)
    ggsave(
      file = plotpathoneside,
      plot = plotoneside,
      width = GGWIDTH, height = GGHEIGHT
    )
  } else {
    return(list(twoside = plottwoside, oneside = plotoneside))
  }
}

powerplot(saveplot = TRUE)
