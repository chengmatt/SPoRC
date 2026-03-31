# Restructure composition values, used within a variety of functions to either do Francis reweighting or get observed and expected composition values

Restructure composition values, used within a variety of functions to
either do Francis reweighting or get observed and expected composition
values

## Usage

``` r
Restrc_Comps(Exp, Obs, Comp_Type, age_or_len, AgeingError, comp_agg_type)
```

## Arguments

- Exp:

  Expected values (catch at age or survey index at age) indexed for a
  given year and fleet (structured as a matrix by age and sex)

- Obs:

  Observed values (catch at age or survey index at age) indexed for a
  given year and fleet (structured as a matrix by age and sex)

- Comp_Type:

  Composition Parameterization Type (== 0, aggregated comps by sex, ==
  1, split comps by sex (no implicit sex ratio information), == 2, joint
  comps across sexes (implicit sex ratio information)

- age_or_len:

  Age or length comps (== 0, Age, == 1, Length)

- AgeingError:

  Ageing Error matrix

## Value

Returns a list of array observed and expected values for a given year
and fleet

## Examples

``` r
if (FALSE) { # \dontrun{
comps <- Restrc_Comps(Exp, Obs, Comp_Type, age_or_len, AgeingError)
comps$Exp; comps$Obs
} # }
```
