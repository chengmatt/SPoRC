# Map recruitment seasonal apportionment parameters

Internal helper called by
[`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md)
to construct the TMB/RTMB factor map for `rec_seas_prop_pars`, the
logit-scale parameters controlling the proportion of annual recruitment
assigned to each season. The array has dimensions
`[n_pop x (n_seas - 1)]`, using a sum-to-one soft-max parameterisation
with one reference season omitted.

## Usage

``` r
do_rec_seas_prop_mapping(input_list, rec_seas_prop_spec)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `$data$n_pop`, `$data$n_seas`, and `$data$use_fixed_rec_seas_prop` to
  be set by upstream setup functions.

- rec_seas_prop_spec:

  Character string specifying the seasonal apportionment structure, or
  `NULL` to estimate all proportions independently across populations
  and seasons. Options when non-`NULL`:

  `"est_shared_pop"`

  :   Estimate seasonal proportions but share them across populations,
      so a single set of `n_seas - 1` parameters applies to all
      populations. Only valid for seasonal models (`n_seas > 1`). Also
      resets `use_fixed_rec_seas_prop` to `0` if it was previously `1`.

  `"fix"`

  :   All `rec_seas_prop_pars` fixed at their starting values (mapped to
      `NA`).

  `NULL`

  :   Estimate all `n_pop x (n_seas - 1)` parameters independently. Also
      resets `use_fixed_rec_seas_prop` to `0` if it was previously `1`.

## Value

The input `input_list` with `$map$rec_seas_prop_pars` set to a factor
vector of length `n_pop * (n_seas - 1)`, or `NULL` when `n_seas = 1`.
`$data$use_fixed_rec_seas_prop` may be modified as a side effect when
estimation is requested alongside a previously fixed seasonal proportion
flag.

## Details

When `n_seas = 1`, the parameter is structurally irrelevant and both
`$par$rec_seas_prop_pars` and `$map$rec_seas_prop_pars` are set to
`NULL`. If estimation is requested but `use_fixed_rec_seas_prop = 1`, a
warning is issued and `$data$use_fixed_rec_seas_prop` is automatically
reset to `0`.

When `$data$rec_lag = 0` (age-0 recruitment) and `$data$spawn_seas > 1`,
seasons before `spawn_seas` are structurally fixed at zero by a
restricted softmax in the RTMB model function (recruits can't predate
the spawning event that produced them), so only the first
`n_seas - spawn_seas` columns of `rec_seas_prop_pars` are ever used as
free logits. This function forces the remaining, structurally-unused
trailing columns to `NA` regardless of `rec_seas_prop_spec`, so they
can't silently soak up estimation/gradient.
