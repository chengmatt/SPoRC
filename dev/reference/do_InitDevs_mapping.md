# Map initial age-structure deviation parameters

Internal helper called by
[`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md)
to construct the TMB/RTMB factor map for `ln_InitDevs`, the log-scale
deviations from the equilibrium initial age structure. The `ln_InitDevs`
array has dimensions `[n_pop x n_regions x (n_ages - 1)]`, where the age
dimension excludes the plus group by default (see `equil_init_age_strc`
below).

## Usage

``` r
do_InitDevs_mapping(input_list, InitDevs_spec, rec_dd, init_age_devs_shared)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `$data$equil_init_age_strc`, `$data$rec_region_prop_spec`,
  `$data$natal_region`, and `$data$rec_dd` to be set by upstream setup
  functions.

- InitDevs_spec:

  Character string specifying the sharing structure for `ln_InitDevs`,
  or `NULL` to estimate all deviations independently across all
  dimensions. Options when non-`NULL`:

  `"est_shared_pop_r"`

  :   A single set of age-specific deviations shared across all
      populations and regions. Each age class receives one estimated
      parameter regardless of how many populations or regions are
      modelled. Required when `rec_dd = "global"` and `n_regions > 1`.

  `"est_shared_r"`

  :   Separate age-specific deviations per population, shared across
      regions within each population. Regions within the same population
      are constrained to identical initial age structure. Also valid
      under global density dependence.

  `"fix"`

  :   All `ln_InitDevs` parameters fixed at zero (mapped to `NA`).
      Equivalent to assuming a fully deterministic initial age
      structure.

  `NULL`

  :   Estimate all deviations independently across populations, regions,
      and ages. Not permitted when `rec_region_prop_spec = 1` and
      `n_pop > 1`, as non-natal regions have no recruitment.

- rec_dd:

  Recruitment density-dependence structure inherited from
  [`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md).
  `"global"` restricts valid `InitDevs_spec` choices to `"est_shared_r"`
  or `"est_shared_pop_r"` when `n_regions > 1`.

- init_age_devs_shared:

  Integer vector of length `n_ages - 1` specifying an explicit
  parameter-sharing structure for `ln_InitDevs` along the age dimension.
  Each element gives the factor level assigned to that age position;
  positions sharing the same integer value are constrained to a single
  estimated parameter. Used in conjunction with
  `equil_init_age_strc = 3` (`"stoch_shared_ages"`), which activates
  user-defined age sharing while still estimating deviations
  independently across populations and regions. The sharing structure is
  also respected by `InitDevs_spec` options: `"est_shared_r"` applies
  the vector per population (with a population-level offset so pops
  remain independent), and `"est_shared_pop_r"` applies it globally (no
  offset, all pops and regions share the same parameters). A typical use
  case is replicating ADMB models where ages beyond the data plus group
  share the last estimated deviation, e.g. `c(1:42, rep(42, 9))` for a
  52-age model with 43 data ages, giving 42 free parameters. When `NULL`
  (default), age sharing follows the standard behaviour determined by
  `equil_init_age_strc` alone.

## Value

The input `input_list` with `$map$ln_InitDevs` set to a factor vector of
length `prod(dim(par$ln_InitDevs))`. Active parameters receive
sequential integer indices; plus-group slots (when
`equil_init_age_strc = 1`) and non-natal region slots (when
`rec_region_prop_spec = 1`) are `NA`. Starting values in
`$par$ln_InitDevs` are also reset to `0` for any fixed cells.

## Details

Mapping behaviour is governed by three interacting considerations:

1.  **Equilibrium initialisation** (`equil_init_age_strc`): if `0`, all
    deviations are fixed at zero (no stochastic initial structure). If
    `1`, plus-group deviations are fixed and the remaining ages are
    estimated or shared. If `2`, all ages including the plus group
    receive stochastic deviations.

2.  **Sharing specification** (`InitDevs_spec`): controls whether
    deviations are shared across regions and/or populations.

3.  **No-dispersal constraint**: when `rec_region_prop_spec = 1` and
    `n_pop > 1`, non-natal regions receive no recruitment and their
    initial age deviations are structurally zero; these are
    automatically fixed to `NA` regardless of `InitDevs_spec`, and the
    remaining indices are re-numbered sequentially.
