# The recruitment bias ramp, year by year

The lognormal bias correction a recruitment deviation's penalty is
centered on, \\-b_t \sigma^2 / 2\\, as a factor \\b_t\\ per estimated
year. `do_rec_bias_ramp = 0` gives the full correction in every year
(\\b_t = 1\\); `1` ramps it up, holds it and ramps it down over the four
`bias_year` indices, scaled by `max_bias_ramp_fct`, so a ramp whose
years all sit at the last year is zero everywhere. Used by the objective
and by the setup checks that ask whether the correction touches a given
year.

## Usage

``` r
get_rec_bias_ramp(
  do_rec_bias_ramp,
  bias_year,
  n_est_rec_devs,
  max_bias_ramp_fct = 1
)
```

## Arguments

- do_rec_bias_ramp:

  Integer, `0` or `1`.

- bias_year:

  Integer vector of the four ramp years, as deviation indices.

- n_est_rec_devs:

  Number of estimated recruitment deviation years.

- max_bias_ramp_fct:

  Scale of the ramp at its plateau.

## Value

Numeric vector of length `n_est_rec_devs`.
