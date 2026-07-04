# Map sigma_C (catch observation error SD) parameters

Constructs the `ln_sigmaC` factor map used by the TMB/RTMB objective
function to share or fix the log-scale standard deviation of catch
observation error across regions, years, seasons, and fleets. All cells
within a shared group are assigned the same estimation index.

## Usage

``` r
do_sigmaC_mapping(input_list, sigmaC_spec)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists, as constructed
  by upstream setup functions.

- sigmaC_spec:

  Character string controlling the sharing and estimation structure for
  `ln_sigmaC`. One of:

  `"est_all"`

  :   Unique parameter per region × year × season × fleet.

  `"est_shared_r"`

  :   Shared across regions.

  `"est_shared_y"`

  :   Shared across years.

  `"est_shared_seas"`

  :   Shared across seasons.

  `"est_shared_f"`

  :   Shared across fleets.

  `"est_shared_r_y"`

  :   Shared across regions and years.

  `"est_shared_r_seas"`

  :   Shared across regions and seasons.

  `"est_shared_r_f"`

  :   Shared across regions and fleets.

  `"est_shared_y_seas"`

  :   Shared across years and seasons.

  `"est_shared_y_f"`

  :   Shared across years and fleets.

  `"est_shared_seas_f"`

  :   Shared across seasons and fleets.

  `"est_shared_r_y_seas"`

  :   Shared across regions, years, and seasons.

  `"est_shared_r_y_f"`

  :   Shared across regions, years, and fleets.

  `"est_shared_r_seas_f"`

  :   Shared across regions, seasons, and fleets.

  `"est_shared_y_seas_f"`

  :   Shared across years, seasons, and fleets.

  `"est_shared_r_y_seas_f"`

  :   Single parameter shared across all dimensions.

  `"fix"`

  :   All `ln_sigmaC` parameters fixed at starting values (mapped to
      `NA`).

## Value

The input `input_list` with `$map$ln_sigmaC` set to a factor vector of
length `prod(dim(par$ln_sigmaC))`. Each element is an integer estimation
index, or `NA` when `sigmaC_spec = "fix"`.
