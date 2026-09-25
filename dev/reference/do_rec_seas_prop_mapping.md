# Map recruitment seasonal apportionment parameters

Builds the factor map for `rec_seas_prop_pars` `[n_pop x (n_seas - 1)]`,
the logit-scale share of annual recruitment per season, under a softmax
with one reference season omitted. Called by
[`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md).
When `rec_lag = 0` and `spawn_seas > 1` the seasons before `spawn_seas`
are fixed at zero by a restricted softmax in the model, so the trailing
unused columns are forced to `NA` whatever `rec_seas_prop_spec` asks
for.

## Usage

``` r
do_rec_seas_prop_mapping(input_list, rec_seas_prop_spec)
```

## Arguments

- input_list:

  Named list with `$data`, `$par` and `$map`. Requires `$data$n_pop`,
  `$data$n_seas` and `$data$use_fixed_rec_seas_prop`.

- rec_seas_prop_spec:

  Seasonal apportionment structure. `"est_shared_pop"` estimates one set
  of `n_seas - 1` parameters for every population and is only valid when
  `n_seas > 1`. `"fix"` holds every parameter at its starting value.
  `NULL` estimates all `n_pop x (n_seas - 1)` independently. The first
  and last reset `use_fixed_rec_seas_prop` to `0` with a warning if it
  was `1`.

## Value

`input_list` with `$map$rec_seas_prop_pars` set to a factor vector of
length `n_pop * (n_seas - 1)`, or `NULL` when `n_seas = 1`. Both
`$par$rec_seas_prop_pars` and its map are `NULL` in that case.
`$data$use_fixed_rec_seas_prop` may change as a side effect.
