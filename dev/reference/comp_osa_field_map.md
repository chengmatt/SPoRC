# Map a composition data source to its internal-OSA field names

Translates a composition data source identifier into the exact field
names used in the model `data` list and the RTMB-tracked OSA vector
name, following the naming convention used throughout `SPoRC_rtmb.R`
(e.g. `ObsFishAgeComps`, `ISS_FishAgeComps`, `FishAgeComps_Type`,
`ObsFishAgeComps_osa_discrete`).

## Usage

``` r
comp_osa_field_map(comp_source, pop = FALSE, discard = FALSE)
```

## Arguments

- comp_source:

  One of `"FishAge"`, `"FishLen"`, `"SrvAge"`, `"SrvLen"`.

- pop:

  Logical; population-specific composition source.

- discard:

  Logical; discard composition source (only valid for Fish\* sources).

## Value

A named list of field names: `Obs`, `ISS`, `Wt`, `Use`, `Type`,
`LikeType`, `n_fleets_field`.
