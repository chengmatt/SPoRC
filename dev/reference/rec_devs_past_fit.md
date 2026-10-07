# One replicate's recruitment deviations, the fit's first and drawn from its penalty after

Deviations over the first `n_cond_yrs` years are the fit's. After them,
every deviation the fit penalizes is drawn from the penalty
[`get_rec_devs_penalty`](https://chengmatt.github.io/SPoRC/dev/reference/get_rec_devs_penalty.md)
puts on it: independent deviations about minus half their variance times
the year's bias ramp, or about the fit's own mean under
`RecDevs_pen_center = 1`; a random walk from the previous calendar year;
or an AR1 about its centered mean. Each sd is divided by the square root
of the cell's `Wt_Rec`, since a weighted penalty is the same density at
that smaller sd.

## Usage

``` r
rec_devs_past_fit(data, pars, rep, n_cond_yrs, n_yrs, map_RecDevs = NULL)
```

## Arguments

- data, pars, rep:

  The fit's data list, parameter list and report, or one replicate's
  view of them.

- n_cond_yrs:

  Years whose deviations stay at the fit's.

- n_yrs:

  Years the operating model runs.

- map_RecDevs:

  The fit's map for `ln_RecDevs`, whose `NA` cells are fixed values
  rather than a process, or `NULL` when every cell is estimated.

## Value

Array `[n_pop, n_regions, n_yrs]` of log recruitment deviations.

## Details

A deviation the fit leaves unpenalized keeps the fit's value: one mapped
off, the first years under `dont_pen_recdev_first`, a walk's first year
under its diffuse start, a zero weight, or a region with no recruits.
Cells sharing a penalty level take one draw, and years past the fit's
deviations take none, as in the estimation model.
