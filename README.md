## SPoRC: Stochastic Population Model Over Regional Components <a href='https://github.com/chengmatt/SPoRC'><img src='man/figures/SPoRC_hex.png' align="right" style="height:139px;"/></a>

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/github/actions/workflow/status/chengmatt/SPoRC/R-CMD-check.yaml?label=R-CMD-check)](https://github.com/chengmatt/SPoRC/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://img.shields.io/github/actions/workflow/status/chengmatt/SPoRC/pkgdown.yaml?label=pkgdown)](https://github.com/chengmatt/SPoRC/actions/workflows/pkgdown.yaml)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.20543805-blue)](https://doi.org/10.5281/zenodo.20543805)
[![License: GPL-3](https://img.shields.io/badge/License-GPL3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![docs](https://img.shields.io/badge/docs-release-blue)](https://chengmatt.github.io/SPoRC/)
[![docs (dev)](https://img.shields.io/badge/docs-dev-orange)](https://chengmatt.github.io/SPoRC/dev/)
[![codecov](https://img.shields.io/codecov/c/github/chengmatt/SPoRC/master?label=codecov)](https://codecov.io/gh/chengmatt/SPoRC/tree/master)
<!-- badges: end -->
  
### Overview
`SPoRC` is a flexible modeling framework for spatially structured population dynamics. It accounts for stochasticity in vital rates and movement among geographically defined components. The framework supports:

- Integration of multiple data sources  
- Regional structuring  
- Age- and sex-specific processes  

Thus, `SPoRC` is suitable for both single-region and spatial stock assessment applications.

### Installation

`SPoRC` is implemented in `RTMB` and optionally relies on additional packages for plotting and diagnostics.

#### Prerequisites

Ensure the following packages are installed:

```
install.packages("devtools")       # Development tools
install.packages("TMB")            # Template Model Builder
install.packages("RTMB")           # R interface to TMB
TMB:::install.contrib("https://github.com/vtrijoulet/OSA_multivariate_dists/archive/main.zip") # Optional: multivariate OSA distributions
remotes::install_github("fishfollower/compResidual/compResidual") # Optional OSA residuals
```

#### Installing SPoRC
```
pak::pak("chengmatt/SPoRC")
```

For the dev branch, users can run the following:
```
pak::pak("chengmatt/SPoRC@dev-popn-seasons")
```

### Publications
1. Cheng, M.LH., Goethel, D.R., Cunningham, C.J., Hulson, P.F., Ianelli, J.N., Omori, K.L., 2026. The SPoRC Stock Assessment Package: A Generalized Next‐Generation Platform to Assess Spatial, Age and Sex‐Structured Populations. Fish and Fisheries faf.70082. https://doi.org/10.1111/faf.70082
2. Cheng, M.L.H., Thorson, J.T., Goethel, D.R., Cunningham, C.J., 2026. Parsimonious estimation of environment- and demography-dependent movement in spatially stratified stock assessment models. ICES Journal of Marine Science 83, fsag147. https://doi.org/10.1093/icesjms/fsag147
3. Cheng, M.L.H., Miller, T.J., Goethel, D.R., Cunningham, C.J., 2026. Ensuring consistency in spatially-stratified stock assessment models: An analytical solution for equilibrium population structure with connectivity for age-and size-structured models. Fisheries Research 300, 107800. https://doi.org/10.1016/j.fishres.2026.107800

### Code Examples
1. Alaska Sablefish Spatial Closed Loop Simulations: https://github.com/chengmatt/sablefish_cie_sims_2026
2. Alaska Sablefish 2026 CIE Review: https://github.com/dgoethel-noaa/2026_Sablefish_CIE
3. Alaska Sabelfish 2025 Stock Assessment: https://github.com/dgoethel-noaa/2025_Sablefish_SAFE
4. Continuous Time Markov Chain Movement: https://github.com/chengmatt/ctmc_movement
5. Northern Southeast Inside Waters Sablefish MSE: https://github.com/chengmatt/nsei_sablefish_mse

## Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.

