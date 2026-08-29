# Every argument the setup stages accept

Assembled from [`formals()`](https://rdrr.io/r/base/formals.html) and
the package's own Rd, so it covers the API as it currently stands rather
than as it stood when someone last wrote it down. Regenerating the
vignette regenerates this.

## Usage

``` r
option_reference(stages = setup_stage_order(), guide = NULL)
```

## Arguments

- stages:

  Function names to document. Defaults to the eleven setup stages.

- guide:

  Path to the hand-written options guide, checked so the reference can
  report which settings it also discusses. `NULL` skips the check.

## Value

A data frame with one row per argument: `stage`, `argument`, `default`,
`description`, and `in_guide`.

## Examples

``` r
if (FALSE) { # \dontrun{
ref <- option_reference()
subset(ref, stage == "Tagging")
} # }
```
