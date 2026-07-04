# Map fishery catchability parameters

Constructs the factor map for `ln_fish_q`, controlling whether
catchability parameters are estimated independently per region and time
block or shared across regions. Cells with no fishery index observations
(`UseFishIdx == 0`) are automatically mapped to `NA`.

## Usage

``` r
do_fish_q_mapping(input_list, fish_q_spec)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists.

- fish_q_spec:

  Character vector of length `n_fish_fleets`. Options:

  `"est_all"`

  :   Separate catchability per region × block × fleet.

  `"est_shared_r"`

  :   Single catchability shared across regions, unique per block ×
      fleet.

  `"fix"`

  :   All catchability parameters fixed (mapped to `NA`).

## Value

The input `input_list` with `$map$ln_fish_q` set to a factor vector.
