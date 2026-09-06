# Drop blocks a bin restriction has emptied

A restriction can leave a region, year and season with no observations
at all in the bins being fitted, even though the full composition had
plenty. The fitting likelihood already skips such a block, since
normalizing it would divide by zero, but the one-step-ahead packer has
no way to know: it sees the use flag say "there is data here",
normalizes `(0 + addtocomp)` into a flat composition and fits that. The
packer and the evaluator cannot agree on an emptiness test between
themselves, because the evaluator never sees the observations, so the
two are reconciled here instead by clearing the use flag.

## Usage

``` r
drop_empty_fitted_blocks(obs, use, bins_arr, bin_dim, what)
```

## Arguments

- obs:

  Observation array for the data source.

- use:

  Use-flag array for the data source, whose dims are `obs` without its
  bin and sex dimensions.

- bins_arr:

  `[n_obs_bins x n_fleets]` 0/1 array.

- bin_dim:

  Integer. Which dimension of `obs` holds the bins. The sex dimension is
  taken to be the next one, and fleets the last.

- what:

  Character. Data source name, used in the message.

## Value

`use`, with emptied blocks cleared.

## Details

Clearing it changes nothing about the fit: the likelihood was already
contributing zero for those blocks. It only stops the residual routines
inventing an observation that was never there. Anything cleared is
reported, so a restriction that guts a data source is visible rather
than silent.

A block counts as empty when it holds no finite values at all, exactly
as when it sums to zero, since that is the test the likelihood's own
guard applies.
