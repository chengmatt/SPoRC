# Every argument the setup stages accept

Assembled from [`formals()`](https://rdrr.io/r/base/formals.html) and
the package's own Rd, so it covers every argument as the package
currently defines them, not as someone last wrote them down.
Regenerating the vignette regenerates this.

## Usage

``` r
option_reference(
  stages = c("Setup_Mod_Dim", "Setup_Mod_Rec", "Setup_Mod_Biologicals",
    "Setup_Mod_Movement", "Setup_Mod_Tagging", "Setup_Mod_Catch_and_F",
    "Setup_Mod_FishIdx_and_Comps", "Setup_Mod_SrvIdx_and_Comps",
    "Setup_Mod_Fishsel_and_Q", "Setup_Mod_Srvsel_and_Q", "Setup_Mod_DSEM",
    "Setup_Mod_Weighting", "Setup_Mod_Projection"),
  guide = NULL
)
```

## Arguments

- stages:

  Function names to document. Defaults to every setup stage, in the
  order a model is built.

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
