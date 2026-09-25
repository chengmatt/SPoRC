# Set up catchability deviations for a fishery or survey

Builds the deviation array, its sigma and correlation, and refuses the
combinations that cannot be identified. Called by
[`Setup_Mod_Fishsel_and_Q`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Fishsel_and_Q.md)
and
[`Setup_Mod_Srvsel_and_Q`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Srvsel_and_Q.md).

## Usage

``` r
setup_q_devs(
  input_list,
  q_model,
  sigma_q_spec,
  q_rho_spec,
  q_rw_init_sigma,
  q_type,
  prefix,
  fleet_field,
  use_field,
  fleet_label,
  starting_values
)
```

## Arguments

- input_list:

  List with `data`, `par` and `map`. The catchability blocks, the index
  use arrays and `q_type` must already be set.

- q_model:

  Character vector `[n_fleets]`, one of `"none"` (default), `"iid"`,
  `"rw"`, `"ar1"` or `"dsem"`.

- sigma_q_spec:

  Sharing string for the deviation sigma over region and fleet. Default
  `"est_all"`.

- q_rho_spec:

  Sharing string for the ar1 correlation. Default `"est_all"`.

- q_rw_init_sigma:

  Standard deviation of the first estimated year of a random walk. `NA`
  (default) starts the walk at zero under its own sigma, which keeps the
  block catchability as the level of the series.

- q_type:

  Character vector `[n_fleets]` of `"est"`, `"arith"` or `"geo"`, as
  already validated by the caller.

- prefix:

  `"fish"` or `"srv"`.

- fleet_field:

  Name of the fleet count in `data`.

- use_field:

  Name of the index use array in `data`.

- fleet_label:

  Label used in messages.

- starting_values:

  Named list of starting values.

## Value

`input_list` with the deviation parameters, their maps and
`<prefix>_q_model` in `data`.
