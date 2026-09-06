# Growth parameters in effect in one year

Applies the time-varying deviations of one stratum to its base
parameters. Under the log link a parameter in year \\y\\ is \\P
\exp(\delta_y)\\; under the logit link it is kept inside its bounds,
\$\$P_y = lo + (hi - lo)\\\mathrm{logit}^{-1}\\\left(\log\frac{P -
lo}{hi - P} + \delta_y\right)\$\$ A random walk's deviation array holds
the walk's position, so both structures read the same way here and
differ only in their penalty.

## Usage

``` r
get_growth_pars_year(ln_pars, ln_devs, tv_model, tv_link, bounds, y)
```

## Arguments

- ln_pars:

  Log-scale base parameters of the stratum, length `n_gpars` (five for
  the von Bertalanffy form, six with the Richards coefficient last).

- ln_devs:

  Matrix `[n_yrs x n_gpars]` of the stratum's deviations.

- tv_model:

  Integer vector `[n_gpars]`, 0 constant, 1 iid deviations, 2 random
  walk.

- tv_link:

  Integer, 0 log link, 1 logit link within `bounds`.

- bounds:

  Matrix `[n_gpars x 2]` of lower and upper bounds, read under the logit
  link.

- y:

  Year index.

## Value

Natural-scale parameter vector for the year, with a seventh element
`rho` of one when the base has five parameters.
