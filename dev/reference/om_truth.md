# One replicate's true value of a reported quantity or parameter

The operating model's own value where it holds the quantity under the
same name with a replicate dim last, so each replicate is compared with
what it ran on rather than with the fit. Where the fit's array runs
longer along one dim, a projected year say, the operating model's cells
fill the leading part and the fit's the rest; where the operating
model's runs longer, as recruitment deviations that stop before the last
year do, the fit's cells take its leading part. Anything else is the
fit's.

## Usage

``` r
om_truth(om_arr, fit_arr, sim, n_sims)
```

## Arguments

- om_arr:

  The operating model's array, or `NULL`.

- fit_arr:

  The fit's value, or under a joint self test the replicate's draw.

- sim:

  Replicate.

- n_sims:

  Number of replicates, which the operating model's last dim must be.

## Value

`fit_arr` with the operating model's values in it.
