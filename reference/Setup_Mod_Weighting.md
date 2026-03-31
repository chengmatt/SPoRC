# Set up SPoRC model weighting

Set up SPoRC model weighting

## Usage

``` r
Setup_Mod_Weighting(
  input_list,
  Wt_Catch = 1,
  Wt_FishIdx = 1,
  Wt_SrvIdx = 1,
  Wt_Rec = 1,
  Wt_F = 1,
  Wt_Tagging = 1,
  Wt_FishAgeComps,
  Wt_SrvAgeComps,
  Wt_FishLenComps,
  Wt_SrvLenComps
)
```

## Arguments

- input_list:

  List containing data, parameter, and map lists.

- Wt_Catch:

  Either a numeric scalar (lambda) applied to the overall catch dataset
  or an array of lambdas (i.e., weights can change by year and fleet)
  dimensioned by n_regions, n_years, n_fish_fleets.

- Wt_FishIdx:

  Either a numeric scalar (lambda) applied to the overall fishery index
  dataset or an array of lambdas (i.e., weights can change by year and
  fleet) dimensioned by n_regions, n_years, n_fish_fleets.

- Wt_SrvIdx:

  Either a numeric scalar (lambda) applied to the overall survey index
  dataset or an array of lambdas (i.e., weights can change by year and
  fleet) dimensioned by n_regions, n_years, n_srv_fleets.

- Wt_Rec:

  Numeric weight (lambda) applied to the recruitment penalty.

- Wt_F:

  Numeric weight (lambda) applied to the fishing mortality penalty.

- Wt_Tagging:

  Numeric weight (lambda) applied to tagging data.

- Wt_FishAgeComps:

  Numeric weight (lambda) applied to fishery age composition data.

- Wt_SrvAgeComps:

  Numeric weight (lambda) applied to survey age composition data.

- Wt_FishLenComps:

  Numeric weight (lambda) applied to fishery length composition data.

- Wt_SrvLenComps:

  Numeric weight (lambda) applied to survey length composition data.

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
[`Setup_Mod_Tagging()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Tagging.md)
