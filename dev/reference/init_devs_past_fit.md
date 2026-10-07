# One replicate's initial age deviations, drawn from the penalty the fit puts on them

Each initial age deviation the fit penalizes is drawn from that penalty,
as
[`get_init_devs_penalty`](https://chengmatt.github.io/SPoRC/dev/reference/get_init_devs_penalty.md)
writes it: about minus half the early recruitment variance times the
bias ramp at the year the age was born (no correction under a walk), or
about the fit's own mean under `InitDevs_pen_center = 1`, at the early
sigma (its stationary value under an AR1) over the square root of the
cell's `Wt_Init_Rec`. A later sex tied to the first by
`Use_init_sex_pen` is drawn from both penalties together, given the
first sex's draw.

## Usage

``` r
init_devs_past_fit(data, pars, rep, map_InitDevs = NULL)
```

## Arguments

- data, pars, rep:

  The fit's data list, parameter list and report, or one replicate's
  view of them.

- map_InitDevs:

  The fit's map for `ln_InitDevs`, whose `NA` cells are fixed, or `NULL`
  when every cell is estimated.

## Value

Array `[n_pop, n_regions, n_ages - 1, n_sexes]` of initial age
deviations.

## Details

Deviations the fit leaves unpenalized keep its values, which is all of
them under `equil_init_age_strc` 0 or 4, and cells sharing a level take
one draw.
