# Initial age deviation penalties

Population and region specific penalties on the initial age deviations
(`ln_InitDevs`), plus the tie holding each later sex's curve near the
first sex's. Called from
[`get_recruitment_penalty`](https://chengmatt.github.io/SPoRC/dev/reference/get_recruitment_penalty.md).

## Usage

``` r
get_init_devs_penalty(
  n_pop,
  n_regions,
  n_ages,
  rec_region_prop_spec,
  rec_region_prop,
  equil_init_age_strc,
  ln_InitDevs,
  init_age_devs_shared,
  ln_sigmaR,
  bias_ramp,
  InitDevs_pen_center = 0,
  init_devs_pen_use = NULL,
  Use_init_sex_pen = 0,
  ln_sigma_init_sex = 0,
  init_bias_ramp = NULL,
  map_ln_InitDevs = NULL
)
```

## Arguments

- n_pop, n_regions, n_ages:

  Dimension sizes.

- rec_region_prop_spec:

  Integer switch; when `1`, populations and regions with a fixed zero
  recruitment proportion are skipped.

- rec_region_prop:

  Array `[pop, region]` of recruitment regional apportionment.

- equil_init_age_strc:

  Integer switch naming which initial age deviations are penalized (`0`:
  none, `1`: all but the plus group, `2`: all, `3`: the shared subset).

- ln_InitDevs:

  Array `[pop, region, age, sex]` of initial age deviations. A 3-D
  `[pop, region, age]` array, the layout before the sex dim existed, is
  accepted and treated as one shared curve.

- init_age_devs_shared:

  Integer vector of shared initial age deviation indices, read when
  `equil_init_age_strc == 3`.

- ln_sigmaR:

  Array `[early/late, pop, region]` of log sigma. Initial ages read the
  early one.

- bias_ramp:

  Numeric vector `[year]` of bias ramp adjustment factors.

- InitDevs_pen_center:

  Integer. `1` centers on the deviations' own weighted mean, `0` on the
  bias-corrected mean.

- init_devs_pen_use:

  Array of 0/1 matching `ln_InitDevs`, naming which cells are penalized.
  Sexes sharing one parameter keep only the first sex's copy flagged so
  the shared parameter is not penalized twice; sex-specific deviations
  flag every sex. `NULL` penalizes only the first sex's slice, which is
  the pre-sex-dim behavior.

- Use_init_sex_pen:

  Integer (0/1). Whether each later sex's deviations are tied to the
  first sex's through a Gaussian on their difference at every penalized
  age. Only meaningful when the sexes have their own curves.

- ln_sigma_init_sex:

  Log standard deviation of that tie.

- init_bias_ramp:

  Numeric vector of length `n_ages - 1`, the bias ramp read at the year
  each initial age was born (deviation index `1 - age`). `NULL` reads
  the first model year's ramp value at every age.

- map_ln_InitDevs:

  Numeric array matching `ln_InitDevs` of map levels (`NA` where fixed).
  Cells sharing a level hold one parameter and split one penalty between
  them. `NULL` penalizes every cell in full.

## Value

List with `Init_Rec_nLL` (array `[pop, region, age, sex]`) and
`Init_Sex_nLL` (the same layout, the between-sex tie, zero for the first
sex and whenever the tie is off), each holding negative log-likelihood
penalties and zero where nothing is penalized.

## Details

Each penalized deviation is Gaussian on the log scale with the early
recruitment sigma, \\-\log \phi(d\_{a,s} \mid \mu, \sigma\_{R,1})\\,
where

- \\d\_{a,s}\\ is the deviation at age \\a\\ and sex \\s\\, log scale,
  estimated.

- \\\sigma\_{R,1} = \exp(\code{ln_sigmaR\[1,p,r\]})\\ is the early
  recruitment sigma, log scale.

- \\\mu\\ is the center, either the bias-corrected mean
  \\-\sigma\_{R,1}^2 b_a / 2\\ with \\b_a\\ the bias ramp read at the
  year age \\a\\ was born, or the deviations' own weighted mean pooled
  over ages and sexes.
