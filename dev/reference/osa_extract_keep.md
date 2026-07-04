# Extract the keep indicator from an OSA observation

Returns the per-bin `keep` indicator carried by an OSA-tagged
observation slice. During ordinary fitting (plain numeric input) all
bins are treated as kept.

## Usage

``` r
osa_extract_keep(xobs, n)
```

## Arguments

- xobs:

  Either an object of class `"osa"` or a plain numeric vector.

- n:

  Integer length used to construct the default all-ones indicator when
  `xobs` is not an `"osa"` object.

## Value

A numeric/AD vector of keep indicators (length `n`).
