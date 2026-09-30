# Whether each linked series settles on one variance

The initial age deviations need a single spread, so a series is only
read through its marginal variance when that variance settles. A
moderated sd moves with its moderator and a lagged self path of one or
more never settles; both fall back to `ln_sigmaR`, the prior set for
exactly this.

## Usage

``` r
get_dsem_link_settles(dsem_model, link_col, sd_arrow)
```

## Arguments

- dsem_model:

  Output of
  [`read_dsem_arrows`](https://chengmatt.github.io/SPoRC/dev/reference/read_dsem_arrows.md).

- link_col:

  The grid column each linked series sits in.

- sd_arrow:

  Each series' own sd line, 0 when moderated.

## Value

Integer vector, 1 where the series settles.
