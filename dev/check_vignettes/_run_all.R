# Purpose: run the code in every vignette so a broken chunk shows up as a failure
# Creator: Matthew LH. Cheng
# Date Created: 9/8/26
#
# Usage:
#   Rscript dev/check_vignettes/_run_all.R fast   # reference vignettes, no model fits
#   Rscript dev/check_vignettes/_run_all.R fit    # + case studies and bridges
#   Rscript dev/check_vignettes/_run_all.R all    # + spatial, simulation and closed loop
#   Rscript dev/check_vignettes/_run_all.R slow   # only the vignettes that run for hours
#
#   fast, fit and all are cumulative and default to "fit". slow stands on its own: those
#   vignettes take hours each, so "all" leaves them out and they are asked for by name.
#
#   Each vignette runs in its own Rscript process against devtools::load_all, so a failure
#   stays contained. Per vignette stdout and stderr land in dev/scratch/vignette_logs,
#   which is gitignored. Chunks headed purl = FALSE are skipped, and their number is
#   reported, so fragments cannot pile up unseen.

# Tiers ----------------------------------------------------------------------
# fast: reference and option vignettes, seconds each
tier_fast <- c(
  "a_model_dimensions.Rmd",
  "b_model_parameters.Rmd",
  "c_model_equations.Rmd",
  "d_model_report.Rmd",
  "m_simulation_dimensions.Rmd",
  "aj_development_roadmap.Rmd",
  "contributing.Rmd",
  "architecture.Rmd",
  "ta_option_reference.Rmd",
  "t_model_options.Rmd",
  "v_structuring_your_data.Rmd",
  "k_defining_priors.Rmd",
  "j_starting_mapping.Rmd",
  "q_movement_param.Rmd"
)

# fit: one or two model fits each, writes nothing outside the session
tier_fit <- c(
  "e_single_region_sablefish_case_study.Rmd",
  "f_single_region_ebs_pollock_case_study.Rmd",
  "w_goa_northern_rockfish_case_study.Rmd",
  "x_goa_dusky_rockfish_case_study.Rmd",
  "y_bsai_northern_rockfish_case_study.Rmd",
  "z_bsai_rougheye_rockfish_case_study.Rmd",
  "aa_bsai_pacific_ocean_perch_case_study.Rmd",
  "ab_bsai_atka_mackerel_case_study.Rmd",
  "ac_bsai_northern_rock_sole_case_study.Rmd",
  "ac_wc_sablefish_case_study.Rmd",
  "ad_goa_rex_sole_case_study.Rmd",
  "ae_ebs_pacific_cod_case_study.Rmd",
  "ah_north_sea_sandeel_case_study.Rmd",
  "af_growth_options.Rmd",
  "i_reference_points.Rmd",
  "s_discard_retention.Rmd",
  "u_osa_residuals.Rmd"
)

# all: spatial fits, simulation studies and closed loop, tens of minutes each
tier_all <- c(
  "g_spatial_sablefish_case_study.Rmd",
  "r_natal-homing-pop-lrgr-rg.Rmd",
  "ak_state_space_naa_spatial_sim.Rmd",
  "l_simulation_testing.Rmd",
  "o_get_started.Rmd",
  "h_closed_loop_simulations.Rmd",
  "p_single_region_dusky_alt_mp_testing.Rmd"
)

# vignettes needing data that does not ship with the package. none at present
tier_external <- character(0)

# slow: hours per vignette. six Laplace fits totalling six and a half hours, of which
# the two semi-parametric surfaces are two hours forty and two hours ten on their own
tier_slow <- c(
  "n_single_region_ebs_pollock_randomeff_case_study.Rmd"
)

tiers <- list(
  fast = tier_fast,
  fit = c(tier_fast, tier_fit),
  all = c(tier_fast, tier_fit, tier_all),
  slow = tier_slow
)

# time cap per vignette in seconds, so one stalled fit does not hold the run
caps <- c(fast = 600, fit = 3600, all = 14400, slow = 43200)

args <- commandArgs(trailingOnly = TRUE)
tier <- if(length(args) > 0) args[1] else "fit"
if(!tier %in% names(tiers)) stop("Unknown tier '", tier, "'. Use one of: ", paste(names(tiers), collapse = ", "))

vignettes <- tiers[[tier]]
root <- here::here()
vignette_dir <- file.path(root, "vignettes")
log_dir <- file.path(root, "dev", "scratch", "vignette_logs")
result_dir <- file.path(log_dir, "results")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# a vignette in no tier never runs, so it rots unnoticed. name the orphans here rather
# than letting the directory and the tiers fall out of step in silence
listed <- c(unlist(tiers), tier_external)
orphans <- setdiff(list.files(vignette_dir, pattern = "[.]Rmd$"), listed)
if(length(orphans) > 0) message("In no tier, never run: ", paste(orphans, collapse = ", "))
if(length(tier_external) > 0) message("Needs data outside the package, skipped: ", paste(tier_external, collapse = ", "))

# how much of the tier's code this run covers, so a growing pile of skipped chunks
# does not pass for a clean run
source(file.path(root, "dev", "check_vignettes", "helper-chunks.R"))
counts <- vapply(file.path(vignette_dir, vignettes), function(path) {
  chunks <- read_vignette_chunks(path, keep_fragments = TRUE)
  c(length(chunks), sum(vapply(chunks, function(x) x$fragment, logical(1))))
}, numeric(2))
message(sum(counts[1, ]), " chunks in this tier, ", sum(counts[2, ]), " skipped as fragments")

# Run ------------------------------------------------------------------------
# a vignette that fails does not stop the run. the failures are collected and reported
# at the end, so one broken case study does not hide the state of the rest
results <- data.frame(vignette = vignettes, status = NA_character_, minutes = NA_real_,
                      chunks = NA_character_, line = NA_character_, error = NA_character_,
                      warnings = NA_character_)

message("Running tier '", tier, "': ", length(vignettes), " vignettes")

for(i in seq_along(vignettes)) {

  vignette <- vignettes[i]
  path <- file.path(vignette_dir, vignette)
  stem <- tools::file_path_sans_ext(vignette)
  log_path <- file.path(log_dir, paste0(stem, ".log"))
  result_path <- file.path(result_dir, paste0(stem, ".tsv"))
  unlink(result_path)

  if(!file.exists(path)) {
    message("  [", i, "/", length(vignettes), "] ", vignette, " MISSING")
    results$status[i] <- "missing"
    next
  } # end missing vignette

  message("  [", i, "/", length(vignettes), "] ", vignette)
  exit_code <- system2(
    file.path(R.home("bin"), "Rscript"),
    args = c("--vanilla", shQuote(file.path(root, "dev", "check_vignettes", "run_one.R")),
             shQuote(path), caps[[tier]], shQuote(result_path)),
    stdout = log_path,
    stderr = log_path,
    timeout = caps[[tier]] + 60
  )

  # the worker writes its one row last, so a missing row means it never got there. setTimeLimit
  # cannot interrupt a long call inside compiled code, which is what the wall clock cap is for
  if(!file.exists(result_path)) {
    results$status[i] <- if(exit_code == 124) "over the cap" else "killed"
    message("      ", toupper(results$status[i]), ", see ", log_path)
    next
  } # end no row written

  row <- strsplit(readLines(result_path, warn = FALSE)[1], "\t")[[1]]
  length(row) <- 8
  results$status[i] <- row[2]
  results$minutes[i] <- as.numeric(row[3])
  results$chunks[i] <- row[4]
  results$line[i] <- row[5]
  results$error[i] <- row[7]
  results$warnings[i] <- row[8]

  message("      ", results$status[i], " in ", round(results$minutes[i], 1), " min, ",
          results$chunks[i], " chunks",
          if(!is.na(results$warnings[i]) && results$warnings[i] != "0")
            paste0(", ", results$warnings[i], " distinct warnings") else "")
  if(results$status[i] != "ok") message("      line ", results$line[i], ": ", results$error[i])

} # end i loop

# Report ---------------------------------------------------------------------
write.csv(results, file.path(log_dir, paste0("results_", tier, ".csv")), row.names = FALSE)
failed <- results[!is.na(results$status) & results$status != "ok", ]

message("\nTier '", tier, "' finished in ", round(sum(results$minutes, na.rm = TRUE), 1), " min")
if(nrow(failed) > 0) {
  message(nrow(failed), " vignette(s) failed. Logs in ", log_dir)
  for(j in seq_len(nrow(failed))) {
    message("  ", failed$vignette[j], " (", failed$status[j], ") line ", failed$line[j], ": ", failed$error[j])
  } # end j loop
  quit(status = 1)
} # end failure report

message("All vignettes ran. Logs in ", log_dir)
