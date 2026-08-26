# Set the across-age correlation for one at-age stream

the ICES convention, per stream. Each stream is configured where its
data are configured, so the catch and discard streams are set in
[`Setup_Mod_Catch_and_F`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md)
and the index streams in their own setup functions.

## Usage

``` r
do_age_corr_setup(
  input_list,
  corr,
  stream,
  fleet_field,
  starting_values = list()
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par` and `$map`.

- corr:

  `"iid"` or `"1dar1"`.

- stream:

  Stream tag: `"catch"`, `"discard"`, `"fish_idx"` or `"srv_idx"`.

- fleet_field:

  `"n_fish_fleets"` or `"n_srv_fleets"`.

- starting_values:

  Named list from the caller's `...`.

## Value

`input_list` with the stream's correlation flag and its per-fleet
correlation parameter set.
