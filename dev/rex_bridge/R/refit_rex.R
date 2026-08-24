# Refit the rex sole SPoRC model from the SS3 estimate and save the fit with its
# standard errors, for the bridge figures.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source("tests/testthat/helper-bridge_goa_rex.R")
data("mlt_rg_goa_rex_data"); dat <- mlt_rg_goa_rex_data
input <- seed_goa_rex_mle(suppressWarnings(suppressMessages(build_goa_rex_input(dat))), dat)
t0 <- Sys.time()
seed <- fit_model(input$data, input$par, input$map, do_optim = FALSE, silent = TRUE)
cat("objective at the SS3 estimate:", seed$fn(seed$par), "\n")
est <- fit_model(input$data, input$par, input$map, do_optim = TRUE, newton_loops = 3, silent = TRUE)
cat("objective after refit:", est$optim$objective, " max|grad|:", max(abs(est$gr(est$optim$par))), " in", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
est$sdrep <- RTMB::sdreport(est, hessian.fixed = est$he(est$optim$par))
cat("pdHess:", est$sdrep$pdHess, "\n")
saveRDS(list(input = input, seed_rep = seed$rep, seed_obj = seed$fn(seed$par), est_rep = est$rep, optim = est$optim,
             sdrep_summary = summary(est$sdrep), par_names = names(est$optim$par)),
        "dev/rex_bridge/output/rex_sporc_refit.rds")
cat("saved\n")
