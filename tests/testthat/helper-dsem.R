# Reference density for the dsem tests: the covariance of the grid built densely from the arrows,
# Sigma = (I - B)^-1 t(Gamma) Gamma (I - B)^-T, evaluated as a multivariate normal. No precision, no dgmrf.

dense_dsem_sigma <- function(dsem_beta, ln_dsem_sd, x_grid, dsem_model) {

  n_t <- nrow(x_grid)
  n_cells <- length(x_grid)
  arrow_value <- get_dsem_arrow_values(dsem_beta, ln_dsem_sd, dsem_model)
  arrows <- dsem_model$arrows

  # fill B and Gamma cell by cell
  B <- Gamma <- matrix(0, n_cells, n_cells)
  for(i in seq_len(nrow(arrows))) {
    for(t in seq_len(n_t)) {
      if(t - arrows$lag[i] < 1) next # no earlier value to read
      row <- (arrows$to_idx[i] - 1) * n_t + t
      col <- (arrows$from_idx[i] - 1) * n_t + t - arrows$lag[i]

      # a moderated arrow takes its moderating series' value in the year it points to
      value <- arrow_value[i]
      if(arrows$mod_idx[i] > 0) {
        value <- x_grid[t, arrows$mod_idx[i]]
        if(arrows$type[i] != "path" && isTRUE(dsem_model$mod_var_logscale)) value <- exp(value)
      }

      if(arrows$type[i] == "path") B[row, col] <- value
      else Gamma[row, col] <- value
    } # end t loop
  } # end i loop

  A <- solve(diag(n_cells) - B) # response of every cell to every innovation
  return(A %*% t(Gamma) %*% Gamma %*% t(A))
}

dense_dsem_nLL <- function(dsem_beta, ln_dsem_sd, x_grid, mu_grid, dsem_model) {
  Sigma <- dense_dsem_sigma(dsem_beta, ln_dsem_sd, x_grid, dsem_model)
  L <- t(chol(Sigma))
  z <- forwardsolve(L, as.vector(x_grid - mu_grid))
  return(sum(log(diag(L))) + 0.5 * sum(z^2) + 0.5 * length(x_grid) * log(2 * pi))
}

# the variance of each unknown cell given the known ones by the Schur complement of the dense covariance,
# a second route to what get_dsem_margvar takes from the precision's unknown block
dense_dsem_margvar <- function(dsem_beta, ln_dsem_sd, x_grid, dsem_model, known) {
  Sigma <- dense_dsem_sigma(dsem_beta, ln_dsem_sd, x_grid, dsem_model)
  u <- which(!known)
  k <- which(known)
  v <- diag(Sigma)[u]
  if(length(k) > 0) v <- v - rowSums((Sigma[u, k, drop = FALSE] %*% solve(Sigma[k, k, drop = FALSE])) * Sigma[u, k, drop = FALSE])
  out <- rep(0, length(known))
  out[u] <- v
  return(matrix(out, nrow(x_grid), ncol(x_grid)))
}

# small models covering each arrow type, with values for their parameters
dsem_test_cases <- list(
  two_series = list(
    arrows = "env -> rec, 0, b0\nenv -> rec, 1, b1\nenv -> env, 1, rho\nenv <-> env, 0, sd_env\nrec <-> rec, 0, sd_rec",
    variables = c("env", "rec"),
    values = c(b0 = 0.4, b1 = 0.25, rho = 0.55, sd_env = 0.9, sd_rec = 0.65)),
  chain_shared_fixed = list(
    arrows = "a -> b, 0, p1\nb -> c, 0, p2\na -> c, 2, p3\na -> a, 1, r\nb -> b, 1, r\na <-> a, 0, NA, 1\nb <-> b, 0, sd_b\nc <-> c, 0, sd_c",
    variables = c("a", "b", "c"),
    values = c(p1 = 0.7, p2 = -0.4, p3 = 0.2, r = 0.3, sd_b = 0.5, sd_c = 0.3)),
  same_year_loop = list(
    arrows = "a -> b, 0, p1\nb -> a, 0, p2\na -> a, 1, r\na <-> a, 0, sd_a\nb <-> b, 0, sd_b",
    variables = c("a", "b"),
    values = c(p1 = 0.5, p2 = -0.3, r = 0.4, sd_a = 0.8, sd_b = 0.6)),
  covariance = list(
    arrows = "a -> b, 1, p\na -> a, 1, r\na <-> a, 0, sd_a\nb <-> b, 0, sd_b\na <-> b, 0, c_ab",
    variables = c("a", "b"),
    values = c(p = 0.5, r = 0.4, sd_a = 0.8, sd_b = 0.6, c_ab = 0.3)),
  random_walk = list(
    arrows = "rec -> rec, 1, NA, 1\nrec <-> rec, 0, sd_rec",
    variables = "rec",
    values = c(sd_rec = 0.7))
)

# models with an arrow named after a series, so its value changes year to year. mod_var_logscale
# says whether a moderated sd is the exponential of its series, and cpp_pars is dsem's beta_z order
dsem_moderated_cases <- list(
  moderated_path = list(
    arrows = "env -> env, 1, rho\nb -> b, 1, rho_b\nenv -> rec, 0, b\nenv <-> env, 0, sd_env\nb <-> b, 0, sd_b\nrec <-> rec, 0, sd_rec",
    variables = c("env", "b", "rec"),
    values = c(rho = 0.5, rho_b = 0.3, sd_env = 0.9, sd_b = 0.4, sd_rec = 0.7),
    mod_var_logscale = FALSE,
    cpp_pars = c("rho", "rho_b", "sd_env", "sd_b", "sd_rec")),
  moderated_variance = list(
    arrows = "regime -> regime, 1, rho_r\nregime <-> regime, 0, sd_regime\nrec <-> rec, 0, regime",
    variables = c("regime", "rec"),
    values = c(rho_r = 0.4, sd_regime = 0.6),
    mod_var_logscale = TRUE,
    cpp_pars = c("rho_r", "sd_regime")),
  moderated_both = list(
    arrows = "env -> env, 1, rho\nb -> b, 1, rho_b\nenv -> rec, 0, b\nenv <-> env, 0, sd_env\nb <-> b, 0, sd_b\nrec <-> rec, 0, b",
    variables = c("env", "b", "rec"),
    values = c(rho = 0.5, rho_b = 0.3, sd_env = 0.9, sd_b = 0.4),
    mod_var_logscale = TRUE,
    cpp_pars = c("rho", "rho_b", "sd_env", "sd_b"))
)

# parameter vectors of a test case, in the order the reader numbers them
dsem_case_pars <- function(dsem_model, values) {
  list(dsem_beta = unname(values[dsem_model$beta_names]),
       ln_dsem_sd = unname(log(values[dsem_model$ln_sd_names])))
}

# Hold named arrow parameters at a value the way the interface does it: the arrow's name becomes NA and
# the value goes in the fourth field. Arrows come in as a vector or one string, one arrow per line.
hold_arrows <- function(arrows, values) {
  if(is.null(values)) return(arrows)
  lines <- trimws(sub("#.*$", "", unlist(strsplit(paste(arrows, collapse = "\n"), "\n"))))
  fields <- lapply(strsplit(lines[nzchar(lines)], ","), trimws)
  vapply(fields, function(f) {
    if(length(f) >= 3 && f[3] %in% names(values)) paste(c(f[1:2], "NA", format(values[[f[3]]], digits = 15)), collapse = ", ")
    else paste(f, collapse = ", ")
  }, "")
}
