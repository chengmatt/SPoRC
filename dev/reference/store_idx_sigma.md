# Store the estimated part of an index sd

The operating model keeps the reported errors in `Obs*_SE` and draws at
their combination with the estimated part, the way the estimation
model's `build_idx_sd` forms its sd, so the reported errors it hands a
refit are the ones the fit read.

## Usage

``` r
store_idx_sigma(sim_list, data_name, form, ln_sigma, n_fleets)
```

## Arguments

- sim_list:

  Simulation list.

- data_name:

  `"FishIdx"`, `"FishIdx_pop"`, `"SrvIdx"` or `"SrvIdx_pop"`.

- form:

  Code 0 to 3, see
  [`combine_idx_sd`](https://chengmatt.github.io/SPoRC/dev/reference/combine_idx_sd.md).

- ln_sigma:

  Log of the estimated part, one per fleet, or `NULL` under form 0.

- n_fleets:

  Number of fleets.

## Value

`sim_list` with `sigma<data_name>_form` and `ln_sigma<data_name>`.
