# Helper that defines the names of the deviation series a dsem can link

A series is one time series of deviations, at fixed indices of every dim
but the year: recruitment on three regions is three series, one per
region, named `rec_Pop_1_Region_1` through `rec_Pop_1_Region_3`. An
arrow line names the series it acts on, and this writes those names from
array indices rather than by hand.

## Usage

``` r
dsem_series(input_list, process, ..., estimated = TRUE)
```

## Arguments

- input_list:

  List with `data`, `par` and `map`, after the setup function for that
  process has run.

- process:

  Deviations to name: `"rec"`, `"growth"`, `"growth_semipar"`, `"NAA"`,
  `"move"`, `"fish_q"` or `"srv_q"`.

- ...:

  Indices to keep, named by dim in any case: `pop`, `region`, `from`,
  `to`, `seas`, `age`, `sex`, `par` or `fleet`, whichever of those dims
  the process has. They are array indices and not labels, so `age = 3`
  is the third model age. A dim left out keeps every level.

- estimated:

  Whether to leave out series the map fixes, which an arrow has no
  parameter to link. `TRUE` by default.

## Value

Character vector of series names, in the order `Setup_Mod_DSEM` reads
them, and empty when the model has no such series, as movement on one
region.

## Details

Under CTMC movement a deviation sits on a region's preference rather
than on a pair of regions, so `to` has one level and `from` is the
region itself.

## See also

[`Setup_Mod_DSEM`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_DSEM.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  dsem_series(input_list, "rec")                          # every recruitment series
  s <- dsem_series(input_list, "move", from = 1, to = 2)   # one exchange, every age
  paste0(s, " <-> ", s, ", 0, sd_move")                    # arrow lines sharing one sd
} # }
```
