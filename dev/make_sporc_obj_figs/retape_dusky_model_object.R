# Purpose: Re-tape the packaged dusky model object after a model function's arguments change
# Creator: Matthew LH. Cheng
# Date Created: 11/4/26
#
# data/dusky_rtmb_model.rda keeps its objective by value, so a change to a model function's arguments
# breaks it even when the model is the same. This rebuilds the objective from the object's own data,
# parameters and mapping, and checks the stored optimum is reproduced before saving over it.

library(here)
devtools::load_all(here())

load(here("data", "dusky_rtmb_model.rda"))
old <- dusky_rtmb_model

# fields the current model no longer reads
data <- old$data
data$t_fish <- NULL

new <- fit_model(data, old$parameters, old$mapping, random = NULL, do_optim = FALSE, silent = TRUE)

# the stored optimum predates named parameter vectors, and the reference point uncertainty reads the names
par <- old$optim$par
stopifnot(length(par) == length(new$par))
names(par) <- names(new$par)

fn_new <- new$fn(par)
grad <- max(abs(new$gr(par)))
cat("objective, stored vs re-taped:", old$optim$objective, fn_new, "\n")
cat("max |gradient| at the stored optimum:", grad, "\n")
stopifnot(abs(fn_new - old$optim$objective) < 1e-8)
stopifnot(grad < 1e-6)

new$rep <- new$report(par)
for(q in c("SSB", "Rec", "NAA", "CAA", "PredSrvIdx")) {
  cat(q, "max relative difference:", max(abs(new$rep[[q]] - old$rep[[q]]) / pmax(abs(old$rep[[q]]), 1e-12)), "\n")
}

new$optim <- old$optim
new$optim$par <- par
new$sdrep <- RTMB::sdreport(new, par.fixed = par)
new$provenance <- old$provenance

dusky_rtmb_model <- new
usethis::use_data(dusky_rtmb_model, overwrite = TRUE)
