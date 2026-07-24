# Map fishery or survey catchability parameters

Constructs the factor map for `ln_fish_q` or `ln_srv_q`, controlling
whether catchability parameters are estimated independently per region
and time block or shared across regions. Cells with no index
observations (aggregated or population-specific) are automatically
mapped to `NA`. Consolidates what were previously two near-identical
functions, `do_fish_q_mapping` and `do_srv_q_mapping`.

## Usage

``` r
do_q_mapping(input_list, q_spec, prefix, fleet_field, fleet_label)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists.

- q_spec:

  Character vector of length `n_fish_fleets`/`n_srv_fleets`. Options:

  `"est_all"`

  :   Separate catchability per region × block × fleet.

  `"est_shared_r"`

  :   Single catchability shared across regions, unique per block ×
      fleet.

  `"fix"`

  :   All catchability parameters fixed (mapped to `NA`).

- prefix:

  Character, either `"fish"` or `"srv"`. Derives the data/parameter
  field names: `ln_<prefix>_q`, `<prefix>_q_blocks`, `Use<Fish/Srv>Idx`,
  `Use<Fish/Srv>Idx_pop`.

- fleet_field:

  Character. Name of the `$data` field giving the number of fleets.
  `"n_fish_fleets"` or `"n_srv_fleets"`.

- fleet_label:

  Character. Used only in the collected setup message, e.g.
  `"fishery fleet"` or `"survey fleet"`.

## Value

The input `input_list` with `$map$ln_<prefix>_q` set to a factor vector.
