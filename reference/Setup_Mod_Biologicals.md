# Setup biological inputs for estimation model

Setup biological inputs for estimation model

## Usage

``` r
Setup_Mod_Biologicals(
  input_list,
  WAA,
  WAA_fish = NULL,
  WAA_srv = NULL,
  MatAA,
  addtocomp = 0.001,
  addtofishidx = 1e-04,
  addtosrvidx = 1e-04,
  addtotag = 1e-10,
  AgeingError = NULL,
  Use_M_prior = 0,
  M_prior = NA,
  fit_lengths = 0,
  SizeAgeTrans = NA,
  Selex_Type = "age",
  M_spec = "est_ln_M",
  M_ageblk_spec = "constant",
  M_regionblk_spec = "constant",
  M_yearblk_spec = "constant",
  M_sexblk_spec = "constant",
  Fixed_natmort = NULL,
  ...
)
```

## Arguments

- input_list:

  List containing data, parameter, and map lists for the model.

- WAA:

  Numeric array of weight-at-age (spawning), dimensioned
  `[n_regions, n_years, n_ages, n_sexes]`.

- WAA_fish:

  Numeric array of weight-at-age (fishery), dimensioned
  `[n_regions, n_years, n_ages, n_sexes, n_fish_fleets]`.

- WAA_srv:

  Numeric array of weight-at-age (survey), dimensioned
  `[n_regions, n_years, n_ages, n_sexes, n_srv_fleets]`.

- MatAA:

  Numeric array of maturity-at-age, dimensioned
  `[n_regions, n_years, n_ages, n_sexes]`.

- addtocomp:

  Numeric value for a constant to add to composition data. Default is
  1e-3. Not used if logistic normal likelihoods are utilized.

- addtofishidx:

  Numeric value for a constant to add to composition data. Default is
  1e-4.

- addtosrvidx:

  Numeric value for a constant to add to composition data. Default is
  1e-4.

- addtotag:

  Numeric value for a constant to add to composition data. Default is
  1e-10

- AgeingError:

  Numeric matrix or array representing the ageing error transition
  matrix. If a matrix (2D), dimensions should be
  `[number of modeled ages, number of observed composition ages]` and
  the ageing error is assumed to be constant over time. If an array
  (3D), dimensions should be
  `[number of years, number of modeled ages, number of observed composition ages]`
  allowing ageing error to vary by year. Defaults to an identity matrix
  (no ageing error) if not specified, assuming observed age bins exactly
  match modeled age bins.

  \*\*Note:\*\* If the observed age composition bins differ from the
  modeled age bins (e.g., observed ages 2–10 while modeled ages are
  1–10), the default identity matrix will cause a dimensional mismatch
  and misalignment. In such cases, users should provide a custom ageing
  error matrix mapping modeled to observed ages. For example, to drop
  the first modeled age bin, supply a matrix like `diag(1, 10)[, 2:10]`.
  This ensures proper alignment of age bins for likelihood calculations.

- Use_M_prior:

  Integer flag indicating whether to apply a natural mortality prior
  (`0` = no, `1` = yes).

- M_prior:

  Numeric vector of length two giving the mean (in normal space) and
  standard deviation of the natural mortality prior.

- fit_lengths:

  Integer flag indicating whether to fit length data (`0` = no, `1` =
  yes).

- SizeAgeTrans:

  Numeric array of size-at-age transition probabilities, dimensioned
  `[n_regions, n_years, n_lens, n_ages, n_sexes]`.

- Selex_Type:

  Character string specifying whether selectivity is age or
  length-based. Default is age-based

  - `"length"`: Length-based selectivity.

  - `"age"`: Age-based selectivity

- M_spec:

  Character string specifying natural mortality estimation approach.
  Defaults to `est_ln_M`, which estimates mortality to be invariant, if
  blocks are not specified. Options:

  - `"est_ln_M"`: Estimates natural mortality across the defined natural
    mortality blocks.

  - `"fix"`: Fix all natural mortality parameters using the provided
    array.

- M_ageblk_spec:

  Specification of age blocking for natural mortality estimation. Either
  a character string ("constant") or a list of index vectors, e.g.,
  `list(1:10, 11:30)`, which specifies 2 age blocks for M.

- M_regionblk_spec:

  Specification of regional blocking for natural mortality. Either a
  character string ("constant") or a list of index vectors, e.g.,
  `list(1:3, 4:5)`, which specifies 2 region blocks for M.

- M_yearblk_spec:

  Specification of year blocking for natural mortality. Either a
  character string ("constant") or a list of index vectors, e.g.,
  `list(1:10, 11:30)`, which specifies 2 year blocks for M.

- M_sexblk_spec:

  Specification of sex blocking for natural mortality. Either a
  character string ("constant") or a list of index vectors, e.g.,
  `list(1,2)`, which specifies sex-invariant M.

- Fixed_natmort:

  Numeric array of fixed natural mortality values, dimensioned
  `[n_regions, n_years, n_ages, n_sexes]`. Required if `M_spec = "fix"`.

- ...:

  Additional arguments for starting values such as `ln_M` and
  `M_offset.` These are ignored if `M_spec = fix`.

## See also

Other Model Setup:
[`Setup_Mod_Catch_and_F()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Catch_and_F.md),
[`Setup_Mod_Dim()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Dim.md),
[`Setup_Mod_FishIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_FishIdx_and_Comps.md),
[`Setup_Mod_Fishsel_and_Q()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Fishsel_and_Q.md),
[`Setup_Mod_Movement()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Movement.md),
[`Setup_Mod_Rec()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Rec.md),
[`Setup_Mod_SrvIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_SrvIdx_and_Comps.md),
[`Setup_Mod_Srvsel_and_Q()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Srvsel_and_Q.md),
[`Setup_Mod_Tagging()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Tagging.md),
[`Setup_Mod_Weighting()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Weighting.md)
