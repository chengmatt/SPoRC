# Map population-specific discard observation error SD parameters

Constructs the `ln_sigmaD_pop` factor map used by the TMB/RTMB objective
function to share or fix the log-scale standard deviation of
population-specific discard observation error across populations,
regions, years, seasons, and fleets. All cells within a shared group are
assigned the same estimation index.

## Usage

``` r
do_sigmaD_pop_mapping(input_list, sigmaD_pop_spec)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists, as constructed
  by upstream setup functions. Requires `$data$n_pop`,
  `$data$n_regions`, `$data$years`, `$data$n_seas`,
  `$data$n_fish_fleets`, and `$par$ln_sigmaD_pop` to be populated before
  calling.

- sigmaD_pop_spec:

  Character string controlling the sharing and estimation structure for
  `ln_sigmaD_pop`. One of:

  `"fix"`

  :   All parameters fixed at starting values (mapped to `NA`).

  `"est_all"`

  :   Unique parameter per population × region × year × season × fleet
      cell.

  `"est_shared_pop"`

  :   Shared across populations; unique per region × year × season ×
      fleet.

  `"est_shared_r"`

  :   Shared across regions; unique per population × year × season ×
      fleet.

  `"est_shared_y"`

  :   Shared across years; unique per population × region × season ×
      fleet.

  `"est_shared_seas"`

  :   Shared across seasons; unique per population × region × year ×
      fleet.

  `"est_shared_f"`

  :   Shared across fleets; unique per population × region × year ×
      season.

  `"est_shared_pop_r"`

  :   Shared across populations and regions.

  `"est_shared_pop_y"`

  :   Shared across populations and years.

  `"est_shared_pop_seas"`

  :   Shared across populations and seasons.

  `"est_shared_pop_f"`

  :   Shared across populations and fleets.

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

  `"est_shared_pop_r_y"`

  :   Shared across populations, regions, and years.

  `"est_shared_pop_r_seas"`

  :   Shared across populations, regions, and seasons.

  `"est_shared_pop_r_f"`

  :   Shared across populations, regions, and fleets.

  `"est_shared_pop_y_seas"`

  :   Shared across populations, years, and seasons.

  `"est_shared_pop_y_f"`

  :   Shared across populations, years, and fleets.

  `"est_shared_pop_seas_f"`

  :   Shared across populations, seasons, and fleets.

  `"est_shared_r_y_seas"`

  :   Shared across regions, years, and seasons.

  `"est_shared_r_y_f"`

  :   Shared across regions, years, and fleets.

  `"est_shared_r_seas_f"`

  :   Shared across regions, seasons, and fleets.

  `"est_shared_y_seas_f"`

  :   Shared across years, seasons, and fleets.

  `"est_shared_pop_r_y_seas"`

  :   Shared across populations, regions, years, and seasons.

  `"est_shared_pop_r_y_f"`

  :   Shared across populations, regions, years, and fleets.

  `"est_shared_pop_r_seas_f"`

  :   Shared across populations, regions, seasons, and fleets.

  `"est_shared_pop_y_seas_f"`

  :   Shared across populations, years, seasons, and fleets.

  `"est_shared_r_y_seas_f"`

  :   Shared across regions, years, seasons, and fleets.

  `"est_shared_pop_r_y_seas_f"`

  :   Single parameter shared across all dimensions.

## Details

The sharing specification encodes which dimensions are collapsed into a
single parameter via underscore-separated tokens. For example,
`"est_shared_pop_r"` shares across populations and regions (one
parameter per year × season × fleet combination), while
`"est_shared_r_f"` shares across regions and fleets (one parameter per
population × year × season combination).
