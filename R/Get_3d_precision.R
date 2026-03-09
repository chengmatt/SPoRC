#' Construct a 3D sparse precision matrix over ages, years, and cohorts
#'
#' Builds a sparse \eqn{(n_{\text{ages}} \times n_{\text{yrs}}) \times
#' (n_{\text{ages}} \times n_{\text{yrs}})} precision matrix \eqn{Q} for a
#' Gaussian Markov random field (GMRF) with simultaneous autoregressive (SAR)
#' structure across three biological dimensions: age, year, and cohort
#' (age-year diagonal). The matrix is constructed via the path-matrix
#' factorisation \eqn{Q = (I - B)^\top \Omega^{-1} (I - B)}, where \eqn{B}
#' encodes the partial correlations and \eqn{\Omega} is a diagonal variance
#' matrix. Two variance parameterisations are supported: marginal (stationary)
#' and conditional (non-stationary).
#'
#' @param n_ages Integer. Number of age classes.
#' @param n_yrs Integer. Number of years.
#' @param pcorr_age Numeric. Partial correlation along the age dimension
#'   (i.e., between adjacent ages within the same year).
#' @param pcorr_year Numeric. Partial correlation along the year dimension
#'   (i.e., between adjacent years within the same age).
#' @param pcorr_cohort Numeric. Partial correlation along the cohort diagonal
#'   (i.e., between the \eqn{(a-1, y-1)} and \eqn{(a, y)} cell).
#' @param ln_var_value Numeric. Log of the target variance. Exponentiated
#'   internally to \eqn{\sigma^2 = \exp(\text{ln\_var\_value})}.
#' @param Var_Type Integer. Variance parameterisation: \code{0} = marginal
#'   (stationary) variance, where diagonal elements of \eqn{\Omega} are
#'   solved recursively via the accumulator \eqn{(I - B)^{-1}} to achieve a
#'   constant marginal variance \eqn{\sigma^2} at every node (slower);
#'   \code{1} = conditional (non-stationary) variance, where all diagonal
#'   elements of \eqn{\Omega} are set to \eqn{\sigma^2} directly (faster).
#'
#' @return A sparse \code{Matrix::sparseMatrix} precision matrix \eqn{Q} of
#'   dimension \eqn{(n_{\text{ages}} \times n_{\text{yrs}}) \times
#'   (n_{\text{ages}} \times n_{\text{yrs}})}, compatible with
#'   \code{RTMB::dgmrf}.
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom methods as
#' @keywords internal
Get_3d_precision <- function(n_ages, n_yrs, pcorr_age, pcorr_year, pcorr_cohort, ln_var_value, Var_Type){

    "c" <- RTMB::ADoverload("c")
    "[<-" <- RTMB::ADoverload("[<-")

    index = expand.grid(seq_len(n_ages), seq_len(n_yrs)) # create index combinations to loop through
    i = j = x = numeric(0) # initialize posiiton to fill in precision matrix
    var_value = exp(ln_var_value) # transform to normal space

    for(n in 1:nrow(index)){
      age = index[n,1] # get age index out of all index combinations
      year = index[n,2] # get year index out of all index combinations
      if(age > 1 ){
        i = c(i, n)
        j = c(j, which(index[,1] == (age-1) & index[,2] == year))
        x = c(x, pcorr_year) # year correaltion indexing
      }
      if(year > 1){
        i = c(i, n)
        j = c(j, which(index[,1]==age & index[,2]==(year-1)) )
        x = c(x, pcorr_age) # age correlation indexing
      }
      if( age>1 & year>1 ){
        i = c(i, n)
        j = c(j, which(index[,1]==(age-1) & index[,2] == (year-1)) )
        x = c(x, pcorr_cohort) # cohort correlation indexing
      }
    } # end n loop

    # create B path matrix
    B = matrix(0, nrow = n_ages * n_yrs, ncol = n_ages * n_yrs)
    B[cbind(i, j)] = x
    B = as(B, "sparseMatrix")

    # identity matrix
    I = as(diag(1, n_ages * n_yrs, n_ages * n_yrs), "sparseMatrix")

    # Solve Omega recursively for stationary variance (accumulator function)
    if(Var_Type == 0) {
      L = solve(I-B) # solve to get accumulator function for stationary variance
      d = rep(0, nrow(index))
      for(n in 1:nrow(index) ){
        if(n==1){
          d[n] = var_value
        }else{
          cumvar = sum(L[n,seq_len(n-1)] * d[seq_len(n-1)] * L[n,seq_len(n-1)])
          d[n] = (var_value-cumvar) / L[n,n]^2
        }
      } # end n loop
    } # end marginal variance (stationary variance)

    if(Var_Type == 1) d = var_value # conditional variance (non-stationary variance)

    # omega matrix
    Omega_inv = diag(1/d, n_ages * n_yrs, n_ages * n_yrs)
    Q = as((I-Matrix::t(B)) %*% Omega_inv %*% (I-B), "sparseMatrix") # solve for precision

    return(Q)
  }



