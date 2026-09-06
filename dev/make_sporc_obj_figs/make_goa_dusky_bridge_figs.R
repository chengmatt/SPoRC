# Purpose: Bridge the 2024 GOA dusky rockfish assessment to SPoRC and render its figures
# Creator: Matthew LH. Cheng
# Date Created: 8/10/26
#
# sgl_rg_dusky_data has no ADMB parameter vector, so there is no stage that seeds SPoRC at the
# assessment's estimate; the comparison is between the two optimized fits

library(here)
library(dplyr)
source(here("tests", "testthat", "helper-bridge_goa_dusky.R"))
source(here("dev", "make_sporc_obj_figs", "helper-bridge_figs.R"))

label <- "2024 GOA Dusky Rockfish Assessment"

dat <- SPoRC::sgl_rg_dusky_data
yrs <- dat$years
n_yrs <- length(yrs)

# Read the assessment's report file ---------------------------------------------
# The report file is a flat label-then-numbers dump, so each quantity is pulled
# by its label rather than by line number.
read_dusky_rep <- function(path) {

  ln <- readLines(path, warn = FALSE)

  row_after <- function(tag) {
    i <- grep(tag, ln, fixed = TRUE)[1]
    txt <- trimws(sub(tag, "", ln[i], fixed = TRUE))
    as.numeric(strsplit(txt, "[[:space:]]+")[[1]])
  }

  # The numbers at age block is a header of ages followed by one row per year,
  # each row led by its year.
  i_naa <- grep("^Numbers[[:space:]]", ln)[1]
  ages <- as.numeric(strsplit(trimws(sub("^Numbers", "", ln[i_naa])), "[[:space:]]+")[[1]])
  naa_rows <- ln[(i_naa + 1):(i_naa + length(yrs))]
  naa <- t(sapply(naa_rows, function(x) {
    as.numeric(strsplit(trimws(x), "[[:space:]]+")[[1]])[-1]
  }, USE.NAMES = FALSE))

  list(
    ages = ages,
    SSB = row_after("SpBiom "),
    Tot_biom = row_after("Tot_biom "),
    Fmort = row_after("Fully_selected_F "),
    Rec = naa[, 1],
    sel_fsh = row_after("Fishery_Selectivity "),
    sel_srv = row_after("TWL Survey_Selectivity ")
  )

} # end read_dusky_rep

admb <- read_dusky_rep(here("dev", "dev_data", "dusky.rep"))

stopifnot(length(admb$SSB) == n_yrs,
          length(admb$Rec) == n_yrs,
          identical(as.numeric(admb$ages), as.numeric(dat$mod_ages)))

# Optimize ----------------------------------------------------------------------
input_list <- build_goa_dusky_input(dat)

devtools::load_all(here('R'))

est <- fit_model(
  input_list$data,
  input_list$par,
  input_list$map,
  random = NULL,
  newton_loops = 3,
  silent = T,
  do_optim = T
)

est$sdrep <- RTMB::sdreport(est)

cat("=== Optimized ===\n")
cat("free parameters:", length(est$optim$par), "\n")
cat("final jnLL:", est$optim$objective, "  max |gradient|:", max(abs(est$gr(est$optim$par))), "\n")
cat("pdHess:", est$sdrep$pdHess, "\n")

# Compare -----------------------------------------------------------------------
sdr <- est$sdrep
rep <- est$rep

ssb <- as.vector(rep$SSB)[1:n_yrs]
rec <- as.vector(rep$Rec)[1:n_yrs]

ggplot2::ggsave(
  here("vignettes", "figures", "x_goa_dusky_ts_comparison.png"),
  bridge_ts_figure(yrs, ssb, rec, admb$SSB, admb$Rec, label,
                                 ssb_se = bridge_se(sdr, "log_SSB", ssb),
                                 rec_se = bridge_se(sdr, "log_Rec", rec)),
  width = 17,
  height = 9,
  dpi = 150
)

# Selectivity is time invariant, so a single curve per gear has everything.
sel_df <- bind_rows(
  bridge_sel_rows(dat$mod_ages, rep$fish_sel[1, 1, 1, 1, , 1, 1], admb$sel_fsh, "Fishery", label),
  bridge_sel_rows(dat$mod_ages, rep$srv_sel[1, 1, 1, 1, , 1, 1], admb$sel_srv, "Survey", label)
)

ggplot2::ggsave(
  here("vignettes", "figures", "x_goa_dusky_sel_comparison.png"),
  bridge_sel_figure(sel_df),
  width = 12,
  height = 7,
  dpi = 150
)

cat("\n=== Optimized SPoRC against the assessment ===\n")
print(rbind(bridge_cmp("SSB", ssb, admb$SSB),
            bridge_cmp("Recruitment", rec, admb$Rec),
            bridge_cmp("Fmort", as.vector(rep$Fmort), admb$Fmort),
            bridge_cmp("Fishery selectivity", as.vector(rep$fish_sel[1, 1, 1, 1, , 1, 1]), admb$sel_fsh),
            bridge_cmp("Survey selectivity", as.vector(rep$srv_sel[1, 1, 1, 1, , 1, 1]), admb$sel_srv)),
      row.names = FALSE, digits = 4)
