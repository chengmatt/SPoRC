# Setup tagging processes and parameters

Setup tagging processes and parameters

## Usage

``` r
Setup_Mod_Tagging(
  input_list,
  UseTagging = 0,
  tag_release_indicator = NULL,
  max_tag_liberty = 0,
  Tagged_Fish = NA,
  Obs_Tag_Recap = NA,
  Tag_LikeType = NA,
  mixing_period = 1,
  t_tagging = 0,
  tag_selex = NA,
  tag_natmort = NA,
  Use_TagRep_Prior = 0,
  TagRep_Prior = NULL,
  move_age_tag_pool = as.list(1:length(input_list$data$ages)),
  move_sex_tag_pool = as.list(1:input_list$data$n_sexes),
  Init_Tag_Mort_spec = NULL,
  Tag_Shed_spec = NULL,
  TagRep_spec = "fix",
  Tag_Reporting_blocks = NULL,
  ...
)
```

## Arguments

- input_list:

  List containing a data list, parameter list, and map list

- UseTagging:

  Numeric (0 or 1) indicating whether to use tagging data (1) or not (0)

- tag_release_indicator:

  Matrix \[n_tag_cohorts x 2\], where columns are release region and
  release year

- max_tag_liberty:

  Maximum number of years to track a tagged cohort

- Tagged_Fish:

  Array \[n_tag_cohorts x n_ages x n_sexes\] describing tagged fish
  releases

- Obs_Tag_Recap:

  Array \[max_tag_liberty x n_tag_cohorts x n_regions x n_ages x
  n_sexes\] observed tag recaptures

- Tag_LikeType:

  Character string specifying tag likelihood type. One of:

  - `"Poisson"`

  - `"NegBin"`

  - `"Multinomial_Release"`

  - `"Multinomial_Recapture"`

  - `"Dirichlet-Multinomial_Release"`

  - `"Dirichlet-Multinomial_Recapture"`

  Example: `Tag_LikeType = "NegBin"`

- mixing_period:

  Numeric indicating minimum years post-release to include in fitting

- t_tagging:

  Fractional year when tagging occurs (e.g., 0.5 for mid-year)

- tag_selex:

  Character string specifying tag recovery selectivity. One of:

  - `"Uniform_DomFleet"`

  - `"SexAgg_DomFleet"`

  - `"SexSp_DomFleet"`

  - `"Uniform_AllFleet"`

  - `"SexAgg_AllFleet"`

  - `"SexSp_AllFleet"`

  Example: `tag_selex = "SexSp_AllFleet"`

- tag_natmort:

  Character string specifying tag natural mortality parameterization.
  One of:

  - `"AgeAgg_SexAgg"`

  - `"AgeSp_SexAgg"`

  - `"AgeAgg_SexSp"`

  - `"AgeSp_SexSp"`

  Example: `tag_natmort = "AgeSp_SexSp"`

- Use_TagRep_Prior:

  Numeric (0 or 1) whether to use tag reporting rate prior

- TagRep_Prior:

  Data frame containing prior specifications for tag reporting
  parameters. Must include columns: `region` (region index), `block`
  (time block index), `mu` (Numeric mean for tag reporting prior (normal
  space); `NA` if symmetric beta is used), `sd` (Numeric standard
  deviation for tag reporting prior (normal space)), and `type` (0 ==
  symmetric beta, 1 == regular beta). Each row specifies a beta prior
  for one tag reporting parameter. Only parameters with rows in this
  data frame will have priors applied.

- move_age_tag_pool:

  List or character specifying pooling of tagging data by age groups.
  Default does not pool ages. Examples:

  - `list(1:5, 6:11, 12:20)` pools these age groups together

  - `"all"` pools all ages together (internally converted to
    `list(1:n_ages)`)

  - `as.list(1:n_ages)` fits each sex separately

- move_sex_tag_pool:

  List or character specifying pooling of tagging data by sex groups.
  Default do not pool sexes. Examples:

  - `list(1:2)` pools sexes together

  - `"all"` pools all sexes together (internally converted to
    `list(1:n_sexes)`)

  - `list(1, 2)` fits each sex separately

- Init_Tag_Mort_spec:

  Character string `"fix"` or `"est"` specifying if initial tag
  mortality is fixed or estimated

- Tag_Shed_spec:

  Character string `"fix"` or `"est"` specifying if chronic tag shedding
  is fixed or estimated

- TagRep_spec:

  Character string specifying tag reporting rate estimation scheme:

  - `"est_all"` estimates rates for all blocks and regions independently

  - `"est_shared_r"` estimates rates shared across regions but varying
    by block

  - `"fix"` fixes all reporting rates (no estimation)

- Tag_Reporting_blocks:

  Character vector specifying blocks of years and regions for tag
  reporting rates. Default is a single block for all regions. Format
  examples:

  - `"Block_1_Year_1-15_Region_1"`

  - `"Block_2_Year_16-terminal_Region_2"`

  - `"none_Region_3"` (means no block, constant for that region)

- ...:

  Additional starting values for tagging parameters such as
  `ln_Init_Tag_Mort`, `ln_Tag_Shed`, `ln_tag_theta`,
  `Tag_Reporting_Pars`

## See also

Other Model Setup:
[`Setup_Mod_Biologicals()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Biologicals.md),
[`Setup_Mod_Catch_and_F()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Catch_and_F.md),
[`Setup_Mod_Dim()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Dim.md),
[`Setup_Mod_FishIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_FishIdx_and_Comps.md),
[`Setup_Mod_Fishsel_and_Q()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Fishsel_and_Q.md),
[`Setup_Mod_Movement()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Movement.md),
[`Setup_Mod_Rec()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Rec.md),
[`Setup_Mod_SrvIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_SrvIdx_and_Comps.md),
[`Setup_Mod_Srvsel_and_Q()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Srvsel_and_Q.md),
[`Setup_Mod_Weighting()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Weighting.md)
