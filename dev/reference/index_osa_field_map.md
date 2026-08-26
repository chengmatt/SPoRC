# Map an index-type data source to its internal-OSA field names

Translates an index-type data source identifier into the exact field
names used in the model `data` list, following the naming convention
used throughout `SPoRC_rtmb.R` (e.g. `ObsFishIdx`, `UseFishIdx`).

## Usage

``` r
index_osa_field_map(index_source, pop = FALSE)
```

## Arguments

- index_source:

  One of `"Catch"`, `"Discard"`, `"FishIdx"`, `"SrvIdx"`, or their
  at-age forms `"CatchAA"`, `"DiscardAA"`, `"FishIdxAA"`, `"SrvIdxAA"`.
  At-age sources return an extra `age` column.

- pop:

  Logical; population-specific index source.

## Value

A named list of field names: `Obs`, `Use`.
