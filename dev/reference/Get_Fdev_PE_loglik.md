# Compute Fishing Mortality Deviation Process Error Log-Likelihood (Negative Scale)

The negative log-likelihood of `ln_F_devs` under an iid, random walk or
AR1 process. Unlike
[`Get_PE_loglik`](https://chengmatt.github.io/SPoRC/dev/reference/Get_PE_loglik.md)
and
[`Get_move_PE_loglik`](https://chengmatt.github.io/SPoRC/dev/reference/Get_move_PE_loglik.md),
which return one positive scalar for the caller to negate, this returns
an already-negated array shaped like `ln_F_devs`, zero where catch is
not used, matching the `Fmort_nLL` reporting.

## Usage

``` r
Get_Fdev_PE_loglik(
  PE_model,
  ln_sigmaF,
  Fdev_rho,
  ln_F_devs,
  map_ln_F_devs,
  Fdev_pen_center = 0
)
```

## Arguments

- PE_model:

  Integer process error structure: `1` iid, `2` random walk with the
  first active year on a diffuse \\N(0, 5)\\, `3` AR1 with the first
  active year drawn from its stationary marginal.

- ln_sigmaF:

  Array `[n_regions x n_seas x n_fish_fleets]` of log-scale process
  error sd.

- Fdev_rho:

  Array `[n_regions x n_seas x n_fish_fleets]` of unconstrained AR1
  partial correlation, transformed to \\(-1, 1)\\ here and read under
  `PE_model == 3` only.

- ln_F_devs:

  Array `[n_regions x n_years x n_seas x n_fish_fleets]` of log-scale
  fishing mortality deviations.

- map_ln_F_devs:

  Array shaped like `ln_F_devs` mirroring `$map$ln_F_devs`: an
  estimation index where a deviation is estimated and `NA` where it is
  fixed. Only estimated deviations are penalized and they alone form the
  active sequence, so a fixed cell is skipped and widens the gap between
  the deviations either side of it.
  [`do_Fmort_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_Fmort_mapping.md)
  builds it from the catch use indicators, estimating a deviation
  wherever aggregated or population-specific catch is used, or where
  `ObsCatch` is `NA`, which is a missing value rather than a real
  closure; a recorded zero is a closure and is excluded.

- Fdev_pen_center:

  Integer. `1` centers on the deviations' own mean, `0` on zero.

## Value

Array shaped like `ln_F_devs` of the negative log-likelihood per cell.

## Details

The walk and the AR1 do not need catch-active years to be contiguous:
the transition between two active years is taken over the elapsed gap
between them, which is the same marginal as estimating deviations for
the closed years and integrating them out, and reduces to the
single-step transition when the gap is one year.
