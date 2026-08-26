# Purpose: Rebuild the exported data objects, bridge figures, and case study
#          objects in dev/make_sporc_obj_figs. Scripts are grouped into tiers by
#          cost and run in dependency order: the data objects write data/*.rda
#          with usethis::use_data, and everything downstream reads those back
#          through SPoRC::<object>, so the data tier always runs first.
#
#          Each script runs in its own Rscript process. That keeps a failure
#          contained to one script, and guarantees the devtools::load_all at the
#          top of a downstream script sees the data/*.rda a previous tier just
#          wrote.
#
# Usage:
#   Rscript dev/make_sporc_obj_figs/_run_all.R data     # data objects only, fast
#   Rscript dev/make_sporc_obj_figs/_run_all.R bridge   # data + bridge figures
#   Rscript dev/make_sporc_obj_figs/_run_all.R objects  # + spatial/case study objects
#   Rscript dev/make_sporc_obj_figs/_run_all.R all      # + long examples and sims
#
#   Tiers are cumulative and default to "bridge". Per-script stdout and stderr
#   land in dev/scratch/make_logs, which is gitignored.
#
# Creator: Matthew LH. Cheng

# Tiers ----------------------------------------------------------------------
# data: no model fits, reads dev/dev_data and writes data/*.rda
tier_data <- c(
  "make_sgl_rg_sabie_data_object.R",
  "make_sgl_dusky_data_object.R",
  "make_goa_northern_data_object.R",
  "make_bsai_nork_data_object.R",
  "make_bsai_pop_data_object.R",
  "make_bsai_rougheye_data_object.R",
  "make_ebs_pollock_data_object.R",
  "make_bsai_atka_data_object.R",
  "make_bsai_nrs_data_object.R",
  "make_three_rg_sablefish_data_spt_comparison.R",
  "make_five_rg_sablefish_data_spt_comparison.R"
)

# bridge: single region fits against an ADMB target, writes vignettes/figures
tier_bridge <- c(
  "make_sgl_rg_sablefish_bridge_figs.R",
  "make_goa_dusky_bridge_figs.R",
  "make_goa_northern_bridge_figs.R",
  "make_bsai_northern_bridge_figs.R",
  "make_bsai_pop_bridge_figs.R",
  "make_bsai_rougheye_bridge_figs.R",
  "make_ebs_pollock_bridge_figs.R",
  "make_bsai_atka_bridge_figs.R",
  "make_bsai_nrs_bridge_figs.R",
  "make_north_sea_sandeel_bridge_figs.R"
)

# objects: multi region and spatial fits, writes data/*.rda and dev/dev_output
tier_objects <- c(
  "make_sgl_rg_sablefish_bridge_objects.R",
  "make_sgl_rg_sablefish_data_spt_comparison.R",
  "make_three_rg_sablefish_spt_comparison.R",
  "make_five_rg_sablefish_spt_comparison.R",
  "make_spatial_sablefish_figs.R",
  "make_dusky_2024_model_object.R"
)

# all: examples and simulation studies, hours of runtime
tier_all <- c(
  "make_full_dusky_asmt_example_plots.R",
  "make_osa_residuals_example.R",
  "make_discarding_example.R",
  "make_reference_pts_proj_plots.R",
  "make_simulation_testing_plots.R",
  "make_sim_test_pop_lrgr_rg.R",
  "make_dusky_closed_loop_simulation_plots.R"
)

tiers <- list(
  data = tier_data,
  bridge = c(tier_data, tier_bridge),
  objects = c(tier_data, tier_bridge, tier_objects),
  all = c(tier_data, tier_bridge, tier_objects, tier_all)
)

args <- commandArgs(trailingOnly = TRUE)
tier <- if(length(args) > 0) args[1] else "bridge"
if(!tier %in% names(tiers)) stop("Unknown tier '", tier, "'. Use one of: ", paste(names(tiers), collapse = ", "))

scripts <- tiers[[tier]]
root <- here::here()
script_dir <- file.path(root, "dev", "make_sporc_obj_figs")
log_dir <- file.path(root, "dev", "scratch", "make_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# Run ------------------------------------------------------------------------
# A script that fails does not stop the run. The failures are collected and
# reported at the end, so one broken bridge does not hide the state of the rest.
results <- data.frame(script = scripts, status = NA_integer_, minutes = NA_real_)

message("Running tier '", tier, "': ", length(scripts), " scripts")

for(i in seq_along(scripts)) {
  script <- scripts[i]
  path <- file.path(script_dir, script)
  log_path <- file.path(log_dir, paste0(tools::file_path_sans_ext(script), ".log"))

  if(!file.exists(path)) {
    message("  [", i, "/", length(scripts), "] ", script, " MISSING")
    results$status[i] <- 127L
    next
  } # end missing script

  message("  [", i, "/", length(scripts), "] ", script)
  started <- Sys.time()
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    args = c("--vanilla", shQuote(path)),
    stdout = log_path,
    stderr = log_path
  )
  results$status[i] <- status
  results$minutes[i] <- as.numeric(difftime(Sys.time(), started, units = "mins"))

  message("      ", if(status == 0) "ok" else paste0("FAILED (", status, ")"),
          " in ", round(results$minutes[i], 1), " min")
} # end i loop

# Report ---------------------------------------------------------------------
failed <- results[!is.na(results$status) & results$status != 0, ]

message("\nTier '", tier, "' finished in ", round(sum(results$minutes, na.rm = TRUE), 1), " min")
if(nrow(failed) > 0) {
  message(nrow(failed), " script(s) failed. Logs in ", log_dir)
  for(j in seq_len(nrow(failed))) message("  ", failed$script[j])
  quit(status = 1)
} # end failure report

message("All scripts succeeded. Logs in ", log_dir)
