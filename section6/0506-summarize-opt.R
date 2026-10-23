### Summarize simulation data, both crossed and nested ###

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
  JOBID <- "2396257"
}

## Packages ##
pkgs <- c(
  "dplyr",
  "readr",
  "ggplot2",
  "dplyr",
  "stringr",
  "knitr",
  "scales",
  "ggplot2"
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
pattern <- paste0("sim.*opt.*", JOBID, ".*\\.csv")
files <- files[grep(pattern, files)]
if (length(files) > 1) stop(paste0("Directory contains ", length(files), " files with JOBID = ", JOBID))
if (length(files) == 0) stop(paste0("No files found for JOBID = ", JOBID))

cat("Loading simulation results from ", files, "\n", sep = "")

## load optimization results
simresults <- read_csv(file.path(resultspath, files), progress = FALSE, show_col_types = FALSE)
simresults$taudiff <- round(simresults$tau2 - simresults$tau1, 4)
eps <- 1e-06
simresults <- simresults %>%
  mutate(
    taudiff = round(tau2 - tau1, 4),
    failed = !is.finite(gradnorm) | gradnorm > eps
  ) 

simname <- sub("^simulationresults-", "", files)
simname <- sub(".csv$", "", simname)


## Summaries ##

get_failure_table <- function(frame) {
  frame %>%
    mutate(
      Failed = sprintf(
        "%s/%s (%.4f%%)",
        formatC(fails, format = "d", big.mark = ","),
        formatC(total, format = "d", big.mark = ","),
        100 * fail_rate
      )
    ) %>%
    arrange(desc(fail_rate)) %>%
    dplyr::select(-fails, -total, -fail_rate) %>%
    tidyr::pivot_wider(names_from = type, values_from = Failed) %>%
    knitr::kable(
      format = "latex",
      digits = 3,
      booktabs = TRUE,
      linesep = "",
      format.args = list(big.mark = ",", scientific = FALSE)
    ) |>
    gsub(
      pattern = "\\\\(toprule|midrule|bottomrule)",
      replacement = ""
    )
}

# Overall failure rate
failure_overall <- simresults %>%
  group_by(type) %>%
  summarize(
    fails = sum(failed),
    total = n(),
    fail_rate = mean(failed),
    .groups = "drop"
  )

# Failure rate by design
failure_design <- simresults %>%
  group_by(type, m, n, r) %>%
  summarize(
    fails = sum(failed),
    total = n(),
    fail_rate = mean(failed),
    .groups = "drop"
  )

# Failure rate by starting value
failure_starting <- simresults %>%
  group_by(type, startingvalue) %>%
  summarize(
    fails = sum(failed),
    total = n(),
    fail_rate = mean(failed),
    .groups = "drop"
  )

# Failure rate by variance components
failure_varcomp <- simresults %>%
  group_by(type, tau1) %>%
  summarize(
    fails = sum(failed),
    total = n(),
    fail_rate = mean(failed),
    .groups = "drop"
  )

failure_varcompdiff <- simresults %>%
  group_by(type, taudiff) %>%
  summarize(
    fails = sum(failed),
    total = n(),
    fail_rate = mean(failed),
    .groups = "drop"
  )

failure_varcompdiffstarting <- simresults %>%
  group_by(type, startingvalue, taudiff) %>%
  summarize(
    fails = sum(failed),
    total = n(),
    fail_rate = mean(failed),
    .groups = "drop"
  )


get_failure_table(failure_overall)
get_failure_table(failure_design)
get_failure_table(failure_starting)
get_failure_table(failure_varcomp)
get_failure_table(failure_varcompdiff)
print(get_failure_table(failure_varcompdiffstarting), n = Inf)
