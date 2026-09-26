# Stage 3 of 3: post fit
#
# Draws a fitted dsem as a picture of its arrows, one node per deviation series and one edge
# per arrow, so the paths, sd lines and covariances can be read off at a glance.

#' Plot a dsem as a graph of its arrows
#'
#' One node per series, and a further node \code{lag(series, k)} for a source read \code{k}
#' years back. A path is a solid edge, blue when its estimate is positive and red when
#' negative; an sd line is a grey loop on its own node and a covariance a grey dashed line
#' with no head. Nodes are laid out in layers, in base graphics.
#'
#' @param x Output of \code{read_dsem_arrows}, or a fitted model whose data holds one
#'   (\code{Setup_Mod_DSEM} then \code{fit_model}).
#' @param sd_rep Optional \code{RTMB::sdreport} of that fit, for the p values. Read from
#'   \code{x$sd_rep} when the fit has one.
#' @param edge_label What to write on each edge: the estimate with stars for its p value
#'   (\code{***} below 0.001, \code{**} below 0.01, \code{*} below 0.05), the estimate
#'   alone, or the parameter name. A fixed arrow shows its value in every case, and a bare
#'   arrow set has no estimates, so it is labeled with names.
#' @param digits Digits to show for an estimate.
#' @param ... Passed to \code{igraph::plot.igraph}, such as \code{vertex.size}.
#'
#' @return The \code{igraph} object, invisibly. Its edge attributes hold each arrow's lag,
#'   estimate, standard error and p value.
#'
#' @seealso \code{\link{Setup_Mod_DSEM}}
#'
#' @export plot_dsem_dag
plot_dsem_dag <- function(x,
                          sd_rep = NULL,
                          edge_label = c("value_and_stars", "value", "name"),
                          digits = 2,
                          ...) {

  if(!requireNamespace("igraph", quietly = TRUE)) {
    stop("plot_dsem_dag needs the igraph package. Install it with install.packages('igraph').")
  }

  edge_label <- match.arg(edge_label)

  # Arrow Values ------------------------------------------------------------

  # a fitted model has a value for every arrow; a bare arrow set only for the fixed ones
  if(!is.null(x$env) && !is.null(x$data$dsem_model)) {

    dsem_model <- x$data$dsem_model
    pars <- x$env$parList()
    arrow_value <- as.numeric(get_dsem_arrow_values(pars$dsem_beta, pars$ln_dsem_sd, dsem_model))
    if(is.null(sd_rep) && !is.null(x$sd_rep)) sd_rep <- x$sd_rep # the fit's own sdreport, when it has one

  } else if(!is.null(x$arrows)) {

    dsem_model <- x
    arrow_value <- ifelse(x$arrows$par == 0 & x$arrows$mod_idx == 0, x$arrows$start, NA_real_) # only the fixed ones
    edge_label <- "name" # nothing is estimated yet, so names are all there is to write

  } else stop("x should be the output of read_dsem_arrows, or a model fitted with a dsem.")

  arrows <- dsem_model$arrows
  n_arrows <- nrow(arrows)
  estimated <- arrows$par > 0 & arrows$mod_idx == 0 # a free arrow with one value over all years
  arrow_value[arrows$mod_idx > 0] <- NA_real_ # a moderated coefficient changes by year, so it has no one value

  # Standard Errors and p Values --------------------------------------------

  se <- rep(NA_real_, n_arrows)

  if(!is.null(sd_rep)) {

    fixed <- summary(sd_rep, "fixed")
    se_beta <- fixed[rownames(fixed) == "dsem_beta", "Std. Error"] # paths and covariances
    se_ln_sd <- fixed[rownames(fixed) == "ln_dsem_sd", "Std. Error"] # sd lines, on the log scale

    for(i in which(estimated)) {
      if(arrows$type[i] == "sd") se[i] <- arrow_value[i] * se_ln_sd[arrows$par[i]] # delta method back off the log scale
      else se[i] <- se_beta[arrows$par[i]]
    } # end i loop

  } # end if an sdreport was given

  p_value <- ifelse(arrows$type != "sd", 2 * stats::pnorm(-abs(arrow_value / se)), NA_real_) # two sided, none for an sd
  stars <- cut(p_value, breaks = c(0, 0.001, 0.01, 0.05, 1), labels = c("***", "**", "*", ""), include.lowest = TRUE)
  stars <- ifelse(is.na(p_value), "", as.character(stars))

  # Edge Labels and Styling -------------------------------------------------

  label <- switch(edge_label,
                  name = arrows$name,
                  value = as.character(round(arrow_value, digits)),
                  value_and_stars = paste0(round(arrow_value, digits), stars))

  label[!estimated] <- as.character(round(arrow_value[!estimated], digits)) # a fixed arrow shows its value
  label[arrows$mod_idx > 0] <- paste0("x ", arrows$name[arrows$mod_idx > 0]) # a moderated one shows its moderator

  # a path is colored by the sign of its estimate, and everything else is grey
  edge_color <- ifelse(arrow_value < 0, "#B2182B", "#2166AC")
  edge_color[arrows$type != "path" | is.na(arrow_value)] <- "grey55"

  edge_lty <- ifelse(arrows$type == "path", 1, 2) # a covariance or sd line is dashed
  edge_head <- ifelse(arrows$type == "path", ">", "-") # a two headed line is drawn with no head

  # a lagged source is its own node, lag(series, k), which keeps a lagged self path from drawing as a loop
  edge_from <- ifelse(arrows$lag == 0, arrows$from, paste0("lag(", arrows$from, ",", arrows$lag, ")"))

  edges <- data.frame(from = edge_from,
                      to = arrows$to,
                      label = label,
                      color = edge_color,
                      lty = edge_lty,
                      arrow.mode = edge_head,
                      lag = arrows$lag,
                      estimate = arrow_value,
                      se = se,
                      p_value = p_value,
                      stringsAsFactors = FALSE)

  # a moderating series sets an arrow's coefficient, so draw it dotted into the series that arrow points at
  moderated <- which(arrows$mod_idx > 0)

  if(length(moderated) > 0) {
    edges <- rbind(edges, data.frame(from = arrows$name[moderated],
                                     to = arrows$to[moderated],
                                     label = "moderates",
                                     color = "grey30",
                                     lty = 3,
                                     arrow.mode = ">",
                                     lag = 0L,
                                     estimate = NA_real_,
                                     se = NA_real_,
                                     p_value = NA_real_,
                                     stringsAsFactors = FALSE))
  } # end if any arrow is moderated

  # Draw It -----------------------------------------------------------------

  nodes <- data.frame(name = union(dsem_model$variables, edges$from)) # the lag nodes are not series, so add them here
  g <- igraph::graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
  coords <- igraph::layout_with_sugiyama(g)$layout # layers, sources above the series they feed

  dots <- list(...)

  # igraph 2.2.2 indexes the node size by node when it places a loop, so one size for all loses every
  # loop but the first node's, and a loop on an outer node needs the margin to point outward into
  dots$vertex.size <- rep(if(is.null(dots$vertex.size)) 15 else dots$vertex.size, length.out = igraph::vcount(g))
  if(is.null(dots$margin)) dots$margin <- 0.3

  do.call(igraph::plot.igraph, c(list(g, layout = coords), dots))

  return(invisible(g))

} # end function
