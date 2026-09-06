# Re-reconcile use flags against freshly simulated observations

[`drop_empty_fitted_blocks`](https://chengmatt.github.io/SPoRC/dev/reference/drop_empty_fitted_blocks.md)
runs once at setup, against the observations the model was built with. A
self test or closed loop replaces those observations replicate by
replicate while with the setup's use flags forward, so under a bin
restriction a simulated replicate can hold nothing in the fitted bins of
a block the flags still call used. This walks the marginal composition
data sources of a data list and reconciles them again.

## Usage

``` r
resync_fitted_blocks(data)
```

## Arguments

- data:

  A data list whose `Obs*` arrays have just been replaced.

## Value

`data`, with its use flags reconciled.

## Details

A no-op when no data source is restricted, which is the usual case.
