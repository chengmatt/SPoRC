# Stage 3 of 3: post fit
#
# Runs test for serial correlation in residuals.

#' Runs Test for Residual Randomness
#'
#' Performs a nonparametric runs test to evaluate whether a sequence of
#' residuals is randomly distributed around a reference mean. The function
#' also computes three–sigma control limits used to identify potential
#' residual outliers.
#'
#' This implementation is adapted from diagnostics used in
#' Stock Synthesis (SS3). The runs test evaluates whether residuals exhibit
#' non-random structure (e.g., positive or negative autocorrelation).
#'
#' @param x Numeric vector of residuals.
#' @param type Character string specifying the assumed mean of the residuals.
#'   If \code{"resid"} (default), the residual mean is assumed to be zero.
#'   Otherwise, the empirical mean of \code{x} is used.
#' @param mixing Character string specifying the alternative hypothesis for
#'   the runs test:
#'   \itemize{
#'     \item \code{"two.sided"} – tests for both positive and negative
#'     autocorrelation (default).
#'     \item \code{"less"} – left-tailed test detecting positive
#'     autocorrelation.
#'   }
#'
#' @returns A list containing:
#' \itemize{
#'   \item \code{sig3lim} – Numeric vector of length two giving the lower and
#'   upper three–sigma control limits for the residuals.
#'   \item \code{p.runs} – P-value from the runs test for randomness.
#' }
#'
#' A small p-value (e.g., \code{< 0.05}) indicates evidence that the residual
#' sequence is not random and may exhibit autocorrelation or other systematic
#' patterns.
#'
#' @export do_runs_test
#' @family Model Diagnostics
#' @import randtests
do_runs_test <- function(x,
                         type = NULL,
                         mixing = "two.sided"
                         ) {

  if(is.null(type)) type="resid"
  if(type=="resid"){
    mu = 0}else{mu = mean(x, na.rm = TRUE)}
  alternative=c("two.sided","left.sided")[which(c("two.sided", "less")%in%mixing)]
  # Average moving range
  mr  <- abs(diff(x - mu))
  amr <- mean(mr, na.rm = TRUE)
  # Upper limit for moving ranges
  ulmr <- 3.267 * amr
  # Remove moving ranges greater than ulmr and recalculate amr, Nelson 1982
  mr  <- mr[mr < ulmr]
  amr <- mean(mr, na.rm = TRUE)
  # Calculate standard deviation, Montgomery, 6.33
  stdev <- amr / 1.128
  # Calculate control limits
  lcl <- mu - 3 * stdev
  ucl <- mu + 3 * stdev
  if(nlevels(factor(sign(x)))>1){
    # Make the runs test non-parametric
    runstest = randtests::runs.test(x,threshold = 0,alternative = alternative)
    if(is.na(runstest$p.value)) p.value =0.001
    pvalue = round(runstest$p.value,3)} else {
      pvalue = 0.001
    }
  return(list(sig3lim=c(lcl,ucl),p.runs= pvalue))
}
