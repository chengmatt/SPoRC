# Gulf of Alaska walleye pollock data for the Rceattle case study

Inputs for the two models proposed for the 2026 Gulf of Alaska walleye
pollock assessment (models 26a and 26b), 1970 to 2026 at ages 1 to 10
with a plus group: catch and age compositions from one fishery, indices
and age compositions from four surveys, the covariate Shelikof Strait
catchability is measured through and the ten ecosystem series of model
26b, with the estimates `Rceattle` returns for the same inputs, used to
check the bridge.

## Usage

``` r
sgl_rg_goa_pollock_data
```

## Format

A list with components `inputs` (dimensions, weight, maturity and
natural mortality at age, the ageing error matrix and the two with young
ages pooled, the catch, index and composition observations with their
use flags and sample sizes, survey and spawn timing, the composition
constant, `Qcov`, the `covariates` data frame, the bottom trawl
catchability prior and the ages `Rceattle` rescales each selectivity
curve at), and `m26a` and `m26b` (`Rceattle`'s spawning biomass,
recruitment, mean recruitment and its deviations, annual log F,
selectivity parameters and the annual offsets of the fishery curve, log
catchability and its annual offsets, the Shelikof catchability effect,
standard deviation and correlation with the latent series `Qcov`
measures (`qcov_latent`), the Dirichlet-multinomial parameters,
predicted indices, marginal negative log likelihood and number of fixed
effects; `m26b` adds the dsem coefficients by arrow name, the covariate
means and the covariate cells). Composition sample sizes and
Dirichlet-multinomial parameters are scaled by one plus the number of
fitted bins times the composition constant, which is how `Rceattle`
counts them; every other estimate is in `Rceattle`'s own terms.

## Source

Monnahan, C., Adams, G., Lauffenberger, N. 2026. Gulf of Alaska pollock
model updates. Alaska Fisheries Science Center, September 2026. Inputs
from <https://github.com/grantdadams/Rceattle-models>, fit with
`Rceattle` 5.23.0 (Adams, G. D., Holsman, K. K., Barbeaux, S. J., Dorn,
M. W., Ianelli, J. N., Spies, I., Stewart, I. J., Punt, A. E. 2022. An
ensemble approach to understand predation mortality for groundfish in
the Gulf of Alaska. Fisheries Research 251: 106303).
