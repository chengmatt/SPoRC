# Flag the selectivity parameters whose deviations are added

Flag the selectivity parameters whose deviations are added

## Usage

``` r
setup_sel_devs_form(x, cont_tv, sel_model, n_fleets, n_bins, what)
```

## Arguments

- x:

  Character vector or list with one element per fleet. Each element is
  `"multiplicative"` or `"additive"` for every parameter of that fleet,
  or one such string per selectivity parameter.

- cont_tv:

  Integer array `[n_regions, n_fleets]` of time variation codes.

- sel_model:

  Integer array `[n_regions, n_yrs, n_fleets]` of selectivity forms.

- n_fleets:

  Number of fleets.

- n_bins:

  Number of selectivity bins, the length of the deviation dim.

- what:

  Argument name for messages.

## Value

0/1 matrix `[n_bins, n_fleets]`, 1 where a parameter's deviation is
added.
