# Compute Fishing Mortality Deviation Process Error Log-Likelihood (Negative Scale)

Calculates the negative log-likelihood contribution for fishing
mortality deviations (`ln_F_devs`) under an iid, random walk, or AR1
process error structure.

## Usage

``` r
Get_Fdev_PE_loglik(
  PE_model,
  ln_sigmaF,
  Fdev_rho,
  ln_F_devs,
  UseCatch,
  UseCatch_pop,
  missing_catch
)
```

## Arguments

- PE_model:

  Integer specifying the process error structure: `1` = IID (deviations
  drawn independently as \\N(0, \sigma^2)\\); `2` = random walk (first
  active year initialized with a diffuse \\N(0, 5)\\ prior); `3` = AR1
  (first active year drawn from its stationary marginal distribution
  \\N(0, \sigma^2 / (1 - \rho^2))\\).

- ln_sigmaF:

  Array `[n_regions x n_seas x n_fish_fleets]` of log-scale process
  error SD.

- Fdev_rho:

  Array `[n_regions x n_seas x n_fish_fleets]` of unconstrained AR1
  partial correlation (only used when `PE_model == 3`); transformed to
  \\(-1, 1)\\ via \\2 / (1 + e^{-2x}) - 1\\.

- ln_F_devs:

  Array `[n_regions x n_years x n_seas x n_fish_fleets]` of log-scale
  fishing mortality deviations.

- UseCatch, UseCatch_pop:

  Same binary catch-usage indicators as elsewhere; a cell is penalized
  if aggregated catch or any population-specific catch is used.

- missing_catch:

  Logical array `[n_regions x n_years x n_seas x n_fish_fleets]`, `TRUE`
  where the aggregate catch observation (`ObsCatch`) is missing (`NA`)
  rather than a true recorded zero. A cell with `UseCatch == 0` (and no
  population-specific catch used) is still penalized as an ordinary
  active year when `missing_catch` is `TRUE` there (fishing presumably
  continued, we simply lack a value to fit), as opposed to a true
  recorded zero, which is treated as a real closure and excluded from
  the gap count.

## Value

Array with the same dimensions as `ln_F_devs`: the negative
log-likelihood contribution per cell.

## Details

**Note:** Unlike
[`Get_sel_PE_loglik`](https://chengmatt.github.io/SPoRC/dev/reference/Get_sel_PE_loglik.md)
and
[`Get_move_PE_loglik`](https://chengmatt.github.io/SPoRC/dev/reference/Get_move_PE_loglik.md),
which return a single positive log-likelihood scalar to be negated by
the caller, this function returns an already-negated array with the same
dimensions as `ln_F_devs` (one value per region/year/season/fleet cell,
`0` where catch is not used), matching the existing `Fmort_nLL`
reporting convention.

Random walk and AR1 do not require catch-active years to be contiguous.
Instead, the transition between two active years is taken over the
elapsed gap \\d\\ between them – exactly the marginal transition you
would get from estimating deviations for the closed years in between and
integrating them out, without actually estimating them:

- Random walk:

  \\\delta_t \mid \delta_s \sim N(\delta_s, d\sigma^2)\\

- AR1:

  \\\delta_t \mid \delta_s \sim N(\rho^d \delta_s, \sigma^2
  \sum\_{i=0}^{d-1} \rho^{2i})\\, where the sum has closed form \\(1 -
  \rho^{2d}) / (1 - \rho^2)\\

Both reduce exactly to the standard single-step transition when \\d =
1\\.
