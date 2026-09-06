# Shared rendering for the bridge scripts in this directory: SSB and recruitment against the assessment
# with percent differences, selectivity against the assessment, and a median/maximum difference table.
#
# Only the comparison layer lives here. Each species script keeps its own specification, likelihood
# crosswalk and diagnostics, which is where the species specific reasoning sits.
#
# Sourced by make_goa_northern, make_goa_dusky, make_bsai_northern, make_bsai_rougheye, make_bsai_pop,
# make_sgl_rg_sablefish and make_ebs_pollock _bridge_figs.R.
# Creator: Matthew LH. Cheng

# standard errors on the natural scale for a log-scale sdreport entry, NA when the quantity was not
# reported at the length asked for. exact = TRUE also returns NA when the sdreport has more entries
bridge_se <- function(sdrep, nm, vals, exact = FALSE) {
  se <- sdrep$sd[names(sdrep$value) == nm]
  if(exact && length(se) != length(vals)) return(rep(NA_real_, length(vals)))
  if(length(se) < length(vals)) return(rep(NA_real_, length(vals)))
  se[seq_along(vals)] * vals
} # end bridge_se

# median and maximum absolute percent difference of a against b. signed = TRUE reports the median
# of the signed difference, which the sablefish and pollock bridges quote
bridge_cmp <- function(lab, a, b, signed = FALSE) {
  dev <- 100 * (a - b) / b
  data.frame(quantity = lab,
             median_pct = if(signed) stats::median(dev) else stats::median(abs(dev)),
             max_abs_pct = max(abs(dev)))
} # end bridge_cmp

# spawning biomass and recruitment against the assessment, beside the percent difference, which
# gets its own panel. mark_year draws a dashed rule; legend_nrow wraps the legend, NULL for one row
bridge_ts_figure <- function(
  yrs,
  ssb,
  rec,
  admb_ssb,
  admb_rec,
  label,
  ssb_se = NA,
  rec_se = NA,
  mark_year = NULL,
  base_size = 20,
  legend_nrow = 2,
  ref_name = "ADMB"
) {

  ts_df <- dplyr::bind_rows(
    data.frame(Year = yrs, value = ssb, se = ssb_se, Par = "Spawning Biomass", type = "SPoRC"),
    data.frame(Year = yrs, value = rec, se = rec_se, Par = "Recruitment", type = "SPoRC"),
    data.frame(Year = yrs, value = admb_ssb, se = NA, Par = "Spawning Biomass", type = label),
    data.frame(Year = yrs, value = admb_rec, se = NA, Par = "Recruitment", type = label)
  )

  p <- ggplot2::ggplot(ts_df, ggplot2::aes(
    x = Year,
    y = value,
    ymin = value - 1.96 * se,
    ymax = value + 1.96 * se,
    color = type,
    fill = type
  )) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~Par, scales = "free", ncol = 1) +
    ggplot2::geom_ribbon(alpha = 0.3, color = NA) +
    ggthemes::scale_color_colorblind() +
    ggthemes::scale_fill_colorblind() +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    ggplot2::theme(legend.position = "top") +
    ggplot2::labs(x = "Year", y = "Value", color = "Type", fill = "Type")

  if(!is.null(legend_nrow)) {
    p <- p + ggplot2::guides(color = ggplot2::guide_legend(nrow = legend_nrow),
                             fill = ggplot2::guide_legend(nrow = legend_nrow))
  } # end if

  pd_df <- dplyr::bind_rows(
    data.frame(Year = yrs, value = 100 * (ssb - admb_ssb) / admb_ssb, Par = "Spawning Biomass"),
    data.frame(Year = yrs, value = 100 * (rec - admb_rec) / admb_rec, Par = "Recruitment")
  )

  p_pd <- ggplot2::ggplot(pd_df, ggplot2::aes(x = Year, y = value)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4) +
    ggplot2::geom_line(color = "#E69F00") +
    ggplot2::geom_point(size = 2, color = "#E69F00") +
    ggplot2::facet_wrap(~Par, scales = "free_y", ncol = 1) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::labs(x = "Year", y = paste0("SPoRC vs ", ref_name, " (%)"))

  if(!is.null(mark_year)) {
    p_pd <- p_pd + ggplot2::geom_vline(xintercept = mark_year, linetype = 2, linewidth = 0.4)
  } # end if

  patchwork::wrap_plots(p, p_pd, ncol = 2, widths = c(1.25, 1))

} # end bridge_ts_figure

# one curve per panel, SPoRC against the assessment. sel_df is long with Age, value, type and the
# column named by facet_by, the gear for a time-invariant fleet and the year for a surface
bridge_sel_figure <- function(
  sel_df,
  facet_by = "Gear",
  ylab = "Selectivity",
  base_size = 20,
  legend = "top",
  nrow = NULL
) {

  p <- ggplot2::ggplot(sel_df, ggplot2::aes(x = Age, y = value, color = type, linetype = type)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by)), nrow = nrow) +
    ggthemes::scale_color_colorblind() +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(legend.position = legend) +
    ggplot2::labs(x = "Age", y = ylab, color = "Type", linetype = "Type")

  if(legend != "none") {
    p <- p + ggplot2::guides(color = ggplot2::guide_legend(nrow = 2),
                             linetype = ggplot2::guide_legend(nrow = 2))
  } # end if

  p

} # end bridge_sel_figure

# The two rows a time invariant fleet needs: SPoRC and the assessment on one
# gear. Saves each species script writing the same four data.frame calls.
bridge_sel_rows <- function(ages, sporc, admb, gear, label, facet_by = "Gear") {
  out <- dplyr::bind_rows(
    data.frame(Age = ages, value = as.vector(sporc), gear = gear, type = "SPoRC"),
    data.frame(Age = ages, value = as.vector(admb), gear = gear, type = label)
  )
  names(out)[names(out) == "gear"] <- facet_by
  out
} # end bridge_sel_rows
