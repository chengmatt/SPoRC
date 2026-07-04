# Map sigma_dmr (discard mortality process error SD) parameters

Constructs the `ln_sigma_dmr` factor map used by the TMB/RTMB objective
function to share or fix the logit-scale standard deviation of discard
mortality process error across regions, seasons, and fleets. All cells
within a shared group are assigned the same estimation index.

## Usage

``` r
do_sigma_dmr_mapping(input_list, sigma_dmr_spec)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists, as constructed
  by upstream setup functions.

- sigma_dmr_spec:

  Character string controlling the sharing and estimation structure for
  `ln_sigma_dmr`. One of:

  `"est_all"`

  :   Unique parameter per region × season × fleet combination.

  `"est_shared_r"`

  :   Shared across regions; unique per season × fleet.

  `"est_shared_seas"`

  :   Shared across seasons; unique per region × fleet.

  `"est_shared_f"`

  :   Shared across fleets; unique per region × season.

  `"est_shared_r_seas"`

  :   Shared across regions and seasons; unique per fleet.

  `"est_shared_r_f"`

  :   Shared across regions and fleets; unique per season.

  `"est_shared_seas_f"`

  :   Shared across seasons and fleets; unique per region.

  `"est_shared_r_seas_f"`

  :   Single parameter shared across all dimensions.

  `"fix"`

  :   All `ln_sigma_dmr` parameters fixed at starting values (mapped to
      `NA`).

## Value

The input `input_list` with `$map$ln_sigma_dmr` set to a factor vector
of length `prod(dim(par$ln_sigma_dmr))`. Each element is an integer
estimation index for shared or estimated configurations, or `NA` when
`sigma_dmr_spec = "fix"`.
