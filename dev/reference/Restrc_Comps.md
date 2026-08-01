# Restructure Composition Values

Restructures observed and expected composition values (catch-at-age or
survey index-at-age) for use in Francis reweighting or other
composition-based analyses. The function can handle aggregated, split,
or joint sex composition parameterizations and can apply ageing error if
provided.

## Usage

``` r
Restrc_Comps(Exp, Obs, Comp_Type, age_or_len, AgeingError)
```

## Arguments

- Exp:

  Array; expected composition values indexed by
  `[region, pop, year, age/length bins, sex, fleet]`.

- Obs:

  Array; observed composition values indexed similarly to `Exp`.

- Comp_Type:

  Integer; composition parameterization type:

  - 0 = aggregated across sexes

  - 1 = split by sex (no implicit sex ratio information)

  - 2 = joint across sexes (implicit sex ratio information)

- age_or_len:

  Integer; 0 for age compositions, 1 for length compositions.

- AgeingError:

  Matrix; ageing error transition matrix (applied to age compositions
  only).

## Value

A list with elements:

- `Exp`: array of expected composition values in observed bins

- `Obs`: array of observed composition values in observed bins
