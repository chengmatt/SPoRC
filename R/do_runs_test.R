#' Runs test function taken from SS3 diags.
#'
#' @param x Vector of residuals
#' @param type Whether to use mean 0 assumption of mean of residuals (default = use mean 0)
#' @param mixing Type of test to do, less = left tailed test that detects positive autocorrelation, two.sided = two sided test that tests whether there is positive and/or negative autocorrealtion. The null is that there isn't any, rejecting the null (<0.05) indictes that there is some non-randomness.
#'
#' @returns List object with p value and limits for a three-sigma limit - (potential data outlier, where residual is > 3 standard deviations away from a mean of 0)
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
