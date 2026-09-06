# Reject a bin restriction that leaves a data source nothing to fit

A composition fitted over a single bin has no information: the
normalized proportion in that bin is identically one whatever the model
says. Every family degenerates, and the routines around them degenerates
further. The logistic-normal families spend one bin as the additive
log-ratio reference and so have no free element left, which gives a
zero-length packed block, a zero-row label frame, and a zero-length
slice request in
[`eval_comp_osa`](https://chengmatt.github.io/SPoRC/dev/reference/eval_comp_osa.md).
The discrete families mark their one bin as the determined cell of the
multinomial, which leaves `get_osa` with nothing to keep and it fails
inside
[`RTMB::oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html).

## Usage

``` r
check_comp_bins_min(bins_arr, like_vals, what)
```

## Arguments

- bins_arr:

  `[n_obs_bins x n_fleets]` 0/1 array from
  [`parse_comp_bins`](https://chengmatt.github.io/SPoRC/dev/reference/parse_comp_bins.md).

- like_vals:

  Integer vector of likelihood codes, one per fleet. `999` marks a fleet
  that is not fitted. `NULL` checks every fleet.

- what:

  Character. Name of the bins argument, used in the error.

## Value

`bins_arr` invisibly.

## Details

Two bins is therefore the minimum, and it is checked at setup where the
message can name the argument and the fleet. Fleets whose likelihood is
`"none"` are skipped, since their bins are never read.
