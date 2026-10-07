# One draw of an F deviation series from the process the estimation model penalizes

The reverse of
[`Get_Fdev_PE_loglik`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Fdev_PE_loglik.md)
for one region, season and fleet. Only estimated years are stepped
through, and a gap between two of them is bridged at the elapsed number
of years. A random walk's first year keeps the replicate's own value,
since the estimation model gives it a diffuse \\N(0, 5)\\ and leaves the
F level to the data, as the recruitment walk does in the operating
model; an AR1's starts at its stationary sd.

## Usage

``` r
draw_F_dev_series(devs, is_est, PE_model, sigma, rho, n_cond)
```

## Arguments

- devs:

  Numeric vector over years, the fit's deviations.

- is_est:

  Logical vector over years, where a deviation is estimated.

- PE_model:

  `1` iid, `2` random walk, `3` AR1.

- sigma:

  Process sd.

- rho:

  AR1 correlation on the natural scale.

- n_cond:

  Leading years that keep the fit's deviations.

## Value

`devs` with the estimated years after `n_cond` drawn.
