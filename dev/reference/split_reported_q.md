# Split a fit's reported catchability into its mean and its deviations

The reported catchability is \\\exp(\ln q\_{blk} + \epsilon)\\. The mean
is the level a process error walks from, the deviations are what the
conditioning years reproduce, and an analytically solved fleet reads no
`ln_q` so its reported value is already its mean.

## Usage

``` r
split_reported_q(rep_q, ln_q, q_blocks, q_type = NULL)
```

## Arguments

- rep_q:

  Reported catchability `[region, year, fleet]`.

- ln_q:

  Log block catchability `[region, block, fleet]`.

- q_blocks:

  Block index `[region, year, fleet]`.

- q_type:

  Integer vector `[fleet]`, zero where catchability is estimated and
  non-zero where it is solved from the index.

## Value

List of `q_mean` and `devs`, both `[region, year, fleet]`.
