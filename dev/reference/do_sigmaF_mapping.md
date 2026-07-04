# Map sigmaF (fishing mortality process error SD) parameters

Constructs the `ln_sigmaF` factor map used by the TMB/RTMB objective
function to share or fix the log-scale standard deviation of fishing
mortality process error across regions, seasons, and fleets. All cells
within a shared group are assigned the same estimation index.

## Usage

``` r
do_sigmaF_mapping(input_list, sigmaF_spec)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists, as constructed
  by upstream setup functions.

- sigmaF_spec:

  Character string controlling the sharing and estimation structure for
  `ln_sigmaF`. One of:

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

  :   All `ln_sigmaF` parameters fixed at starting values (mapped to
      `NA`).

## Value

The input `input_list` with `$map$ln_sigmaF` set to a factor vector of
length `prod(dim(par$ln_sigmaF))`. Each element is an integer estimation
index for shared or estimated configurations, or `NA` when
`sigmaF_spec = "fix"`.
