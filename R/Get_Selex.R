#' Calculate Selectivity
#'
#' Computes selectivity-at-bin using one of several parametric or
#' semi-parametric formulations. Supports constant, parameter-varying,
#' and semi-parametric (GMRF-based) time-varying selectivity structures.
#'
#' @param Selex_Model Integer specifying the parametric selectivity model:
#'   \describe{
#'     \item{0}{Logistic (b50, slope)}
#'     \item{1}{Gamma-shaped dome (bin-at-peak, delta)}
#'     \item{2}{Power function (monotonic decreasing)}
#'     \item{3}{Logistic (b50, b95 parameterization)}
#'     \item{4}{Double-normal dome with plateau and flexible tails}
#'   }
#'
#' @param TimeVary_Model Integer specifying the temporal deviation structure:
#'   \describe{
#'     \item{0}{No time variation}
#'     \item{1}{IID parameter deviations}
#'     \item{2}{Random walk parameter deviations}
#'     \item{3}{3D AR1-GMRF (marginal)}
#'     \item{4}{3D AR1-GMRF (conditional)}
#'     \item{5}{2D AR1-GMRF}
#'   }
#'
#' @param ln_Pars Numeric vector of log-scale selectivity parameters.
#'   Parameter interpretation and ordering depend on `Selex_Model`:
#'   \describe{
#'     \item{Model 0}{(ln b50, ln slope)}
#'     \item{Model 1}{(ln bmax, ln delta)}
#'     \item{Model 2}{(ln power)}
#'     \item{Model 3}{(ln b50, ln b95)}
#'     \item{Model 4}{Six parameters controlling peak, plateau width,
#'                    ascending/descending widths, and tail selectivity}
#'   }
#'
#' @param ln_seldevs Array of selectivity deviations. Expected dimension:
#'   `[n_regions, n_years, n_parameters_or_bins, n_sexes, 1]`.
#'
#'   For `TimeVary_Model` 1–2, deviations apply multiplicatively to the
#'   transformed parametric parameters.
#'
#'   For `TimeVary_Model` 3–5, deviations apply multiplicatively to the
#'   resulting selectivity curve at each bin (semi-parametric formulation).
#'
#' @param Region Integer index of region.
#' @param Year Integer index of year.
#' @param Bin Numeric vector of bins (e.g., age or length).
#' @param Sex Integer index of sex.
#'
#' @return
#' Numeric vector of selectivity values corresponding to `Bin`.
#'
#' @details
#' All parameters are supplied on the log scale and transformed internally
#' to enforce biological constraints (e.g., positivity or bounds in [0,1]).
#'
#' Parametric models (TimeVary_Model 0–2) modify model parameters directly.
#' Semi-parametric models (TimeVary_Model 3–5) apply log-scale deviations
#' multiplicatively to the fully constructed selectivity curve.
#'
#' Output selectivity is not explicitly normalized; downstream model
#' components are responsible for any scaling or standardization.
#'
#' @keywords internal
Get_Selex = function(Selex_Model,
                     TimeVary_Model,
                     ln_Pars,
                     ln_seldevs,
                     Region,
                     Year,
                     Bin,
                     Sex) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  selex = rep(0, length(Bin)) # Temporary container vector

  if(Selex_Model == 0) { # logistic selectivity (b50 and slope)
    # Extract out and exponentiate the parameters here
    b50 = exp(ln_Pars[1]); # b50
    k = exp(ln_Pars[2]); # slope

    if(TimeVary_Model %in% c(1:2)) {
      b50 = b50 * exp(ln_seldevs[Region, Year, 1, Sex, 1]) # b50 parameter varying
      k = k * exp(ln_seldevs[Region, Year, 2, Sex, 1]) # slope parameter varying
    } # end if iid or random walk

    selex = 1 / (1 + exp(-k * (Bin - b50))) # return parmetric form
  }

  if(Selex_Model == 1) { # gamma dome-shaped selectivity
    # Extract out and exponentiate the parameters here
    bmax = exp(ln_Pars[1]) # Bin at max selex
    delta = exp(ln_Pars[2]) # slope parameter

    if(TimeVary_Model %in% c(1:2)) {
      bmax = bmax * exp(ln_seldevs[Region, Year, 1, Sex, 1]) # bmax parameter varying
      delta = delta * exp(ln_seldevs[Region, Year, 2, Sex, 1]) # delta parameter varying
    } # end if iid or random walk

    # Now, calculate/derive power parameter + selex values
    p = 0.5 * (sqrt( bmax^2 + (4 * delta^2)) - bmax)
    selex = (Bin / bmax)^(bmax/p) * exp( (bmax - Bin) / p ) # return parametric form
  }

  if(Selex_Model == 2) { # power function selectivity
    # Extract out and exponentiate the parameters here
    power = exp(ln_Pars[1]); # power parameter

    if(TimeVary_Model %in% c(1:2)) {
      power = power * exp(ln_seldevs[Region, Year, 1, Sex, 1]) # power parameter varying
    } # end if iid or random walk

    selex = 1 / Bin^power # return parametric form
  }

  if(Selex_Model == 3) { # logistic selectivity (b50 and b95)

    # Extract out and exponentiate the parameters here
    b50 = exp(ln_Pars[1]); # b50
    b95 = exp(ln_Pars[2]); # b95

    if(TimeVary_Model %in% c(1:2)) {
      b50 = b50 * exp(ln_seldevs[Region, Year, 1, Sex, 1]) # b50 parameter varying
      b95 = b95 * exp(ln_seldevs[Region, Year, 2, Sex, 1]) # b95 parameter varying
    } # end if iid or random walk

    selex = 1 / (1+19^((b50-Bin)/b95)) # 19 b/c 0.95 / (1 - 0.95) return parametric form
  }

  if(Selex_Model == 4) {

    # define bin ranges for double normal here
    midbin <- Bin

    # Extract and transform parameters here
    p1trans <- min(Bin) + (max(Bin) - min(Bin)) * RTMB::plogis(ln_Pars[1]) # peak bin at plateau
    p2trans <- p1trans + 1 + (0.99 + max(Bin) - p1trans - 1)/(1 + exp(-1.0 * ln_Pars[2])) # width of plateau
    p3trans <- exp(ln_Pars[3]) # ascending width
    p4trans <- exp(ln_Pars[4]) # descending width
    p5trans <- 1/(1 + exp(-1.0 * ln_Pars[5])) # selectivity at first bin
    p6trans <- 1/(1 + exp(-1.0 * ln_Pars[6])) # selectivity at last bin

    if(TimeVary_Model %in% c(1:2)) {
      p1trans = p1trans * exp(ln_seldevs[Region, Year, 1, Sex, 1]) # p1 parameter varying
      p2trans = p2trans * exp(ln_seldevs[Region, Year, 2, Sex, 1]) # p2 parameter varying
      p3trans = p3trans * exp(ln_seldevs[Region, Year, 3, Sex, 1]) # p3 parameter varying
      p4trans = p4trans * exp(ln_seldevs[Region, Year, 4, Sex, 1]) # p4 parameter varying
      p5trans = p5trans * exp(ln_seldevs[Region, Year, 5, Sex, 1]) # p5 parameter varying
      p6trans = p6trans * exp(ln_seldevs[Region, Year, 6, Sex, 1]) # p6 parameter varying
    } # end if iid or random walk

    # construct selectivity function
    asc <- exp(-((midbin - p1trans)^2/p3trans))
    asc.scaled <- (p5trans + (1 - p5trans) * (asc - 0)/(1 - 0))
    desc <- exp(-((midbin - p2trans)^2/p4trans))
    stj <- exp(-((40 - p2trans)^2/p4trans))
    des.scaled <- (1 + (p6trans - 1) * (desc - 1) /(stj - 1))
    join1 <- 1/(1 + exp(-(20 * (midbin - p1trans)/(1 + abs(midbin - p1trans))))) # joiner functions
    join2 <- 1/(1 + exp(-(20 * (midbin - p2trans)/(1 + abs(midbin - p2trans))))) # joiner functions
    selex <- asc.scaled * (1 - join1) + join1 * ((1 - join2) + des.scaled * join2) # return parameteric form
    selex[1] <- p5trans # return parameteric form
  }

  # 3dgmrf model or 2dar1 (sel devs dimensioned as region, year, bin, sex)
  if(TimeVary_Model %in% c(3:5)) selex = selex * exp(ln_seldevs[Region,Year,,Sex, 1]) # varies semi-parametriclly

  return(selex)
} # end function

