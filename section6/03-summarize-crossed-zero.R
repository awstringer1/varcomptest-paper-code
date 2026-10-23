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
  JOBID <- "2358692"
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
homedir <- path.expand("~")
basepath <- file.path(homedir, "work/projects/testing-random-effects/julia")
resultspath <- file.path(basepath, "results")
stopifnot(dir.exists(resultspath))
figurespath <- file.path(basepath, "figures")
stopifnot(dir.exists(figurespath))

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

  
  res <- filter(res, tau1 == 0, tau2 == 0) %>% 
    group_by(m, n, r) %>%
    arrange(.data[[var]], .by_group = TRUE) %>%
    mutate(
      pp = ppoints(n()),
      qtheoretical = qunif(pp),
      qempirical = quantile(.data[[var]], probs = pp)
    )

  grid <- seq(0, 1, length.out = 1e04)
  if (var == "pvaloneside") grid <- seq(0, 0.5, length.out = 1e04)
  kstests <- res %>%
    summarize(
      ks_stat = max(abs(ecdf(.data[[var]])(grid) - grid)),
      ks_stat_p05 = max(abs(ecdf(.data[[var]])(grid) - grid)[grid <= .05]),
      ks_idx = which.max(abs(ecdf(.data[[var]])(grid) - grid)),
      ks_val = grid[which.max(abs(ecdf(.data[[var]])(grid) - grid))],
      .groups = "drop_last"
    )
  
  res <- res %>%
    left_join(kstests, by = c("m", "n", "r"))
  
  qqplt <- ggplot(res, aes(x = qtheoretical, y = qempirical)) +
    theme_bw() +
    facet_grid(
      m + n ~ r,
      labeller = labeller(
        m = label_both,
        n = label_both,
        r = function(x) {
          ifelse(
            x == -1,
            "Balanced",
            paste0("Unbalanced: rho = ", x)
          )
        }
      )
    ) +
    geom_point(alpha = 1, pch = ".") +
    geom_abline(slope = 1, intercept = 0, linewidth = .2) +
    labs(
      title = "Uniform QQ-plot of bootstrap p-values",
      subtitle = H0,
      x = "Unif(0,1) quantiles", y = "Empirical quantiles"
    ) +
    theme(text = element_text(size = GGTEXTSIZE), legend.position = "bottom") +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))

  plotnameqq <- paste0(simname, "-qq-", var, ".pdf")
  
  plotpathqq <- file.path(figurespath, plotnameqq)
  
  if (saveplot) {
    ggsave(file = plotpathqq, plot = qqplt, width = 7, height = 7)
  } else {
  return(list(qqplt))
  }
}

vars <- c("pvalfull", "pvaltau1", "pvaltau2")
for (v in vars) {
  qqpvals(v, saveplot = TRUE)
}

## power plot


powerplot <- function(var, res = pvals, saveplot = FALSE) {
  alpha <- 0.05

  if (!var %in% names(res)) {
    stop("Column '", var, "' does not exist in res.")
  }

  power <- res %>%
    group_by(m, n, r, tau1, tau2) %>%
    summarize(
      nsim = sum(!is.na(.data[[var]])),
      powermean = mean(.data[[var]] <= alpha, na.rm = TRUE),
      powerse = sqrt(powermean * (1 - powermean) / nsim),
      .groups = "drop"
    ) %>%
    mutate(grp = paste0(n, ",", r))

  p <- ggplot(power, aes(x = m)) +
    theme_bw() +
    theme(
      text = element_text(size = GGTEXTSIZE),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    facet_grid(tau1 ~ tau2, labeller = label_both) +
    geom_line(
      aes(
        y = powermean,
        linetype = factor(grp),
        group = factor(grp)
      )
    ) +
    geom_line(
      aes(
        y = powermean + 2 * powerse,
        group = factor(grp)
      ),
      linetype = "solid",
      alpha = 0.2
    ) +
    geom_line(
      aes(
        y = powermean - 2 * powerse,
        group = factor(grp)
      ),
      linetype = "solid",
      alpha = 0.2
    ) +
    labs(
      x = "m",
      y = "Empirical Power",
      linetype = "(n, Balance)"
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, by = 0.25),
      labels = scales::percent_format()
    )

  if (saveplot) {
    plotname <- paste0(simname, "-power-", var, ".pdf")
    plotpath <- file.path(figurespath, plotname)

    ggsave(
      filename = plotpath,
      plot = p,
      width = GGWIDTH,
      height = GGHEIGHT
    )
  }

  invisible(p)
}

vars <- c("pvalfull", "pvaltau1", "pvaltau2")
for (v in vars) {
  powerplot(v, saveplot = TRUE)
}
