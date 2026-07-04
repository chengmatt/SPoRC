# Map discard mortality process error SD parameters

Constructs the `ln_sigma_dmr` factor map used in the TMB/RTMB objective
function to share or fix log-scale standard deviations of discard
mortality process error across regions, seasons, and fleets.

## Usage

``` r
do_dmr_mean_mapping(input_list, dmr_mean_spec)
```

## Arguments

- input_list:

  Named list containing `$data`, `$par`, and `$map`.

- dmr_mean_spec:

  Character string specifying sharing structure. Options include
  `"est_all"`, `"est_shared_r"`, `"est_shared_seas"`, `"est_shared_f"`,
  and combinations thereof, or `"fix"`.

## Value

Updated `input_list` with `$map$ln_sigma_dmr`.

## Details

The mapping assigns integer indices to parameter groups defined by
`dmr_mean_spec`. Cells within the same group share a single estimated
parameter.
