# Keep A Projection's Catch Within some percentage of Last Year's

Builds a constraint for `F_constraint` in
[`Do_Population_Projection`](https://chengmatt.github.io/SPoRC/dev/reference/Do_Population_Projection.md).
Where the retained catch the rule's F would take in a year is more than
`down` below or `up` above last year's catch, F is lowered or raised
until the catch sits at the nearer edge of that band; elsewhere F is
left as the rule set it. The comparison is made on the plain numbers
([`ad_value`](https://chengmatt.github.io/SPoRC/dev/reference/ad_value.md)),
so the constraint also runs inside the model, which keeps the side of
the comparison that was true when it was built.
[`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md)
builds it at the estimates.

## Usage

``` r
catch_stability(down = 0.2, up = 0.2, scope = "region")
```

## Arguments

- down:

  Fraction between 0 and 1. The most the catch may fall below last
  year's. Default 0.2, a catch no lower than 80 percent of last year's.

- up:

  Non-negative fraction. The most the catch may rise above last year's.
  Default 0.2. `Inf` leaves the rise unbounded, so only the floor
  applies.

- scope:

  `"region"` (default) holds each region's own catch to the band of its
  own last year's catch. `"global"` holds the catch summed over regions
  to the band of last year's total, and shares the allowed total across
  regions in proportion to the catch the rule's F would have taken.

## Value

A function `(F, state)` for `F_constraint`.

## Details

Last year's catch is read from the projection, so in projection year 2
it is the terminal year's replayed catch. A region (or, under
`"global"`, a projection) that took no catch last year has no catch to
hold to, and is left at the rule's F. Applied after `catch_cap` in a
list, a floor can raise the catch back above the cap, so put the cap
last when nothing is to exceed it.

## See also

Other Reference Points and Projections:
[`Do_Population_Projection()`](https://chengmatt.github.io/SPoRC/dev/reference/Do_Population_Projection.md),
[`Get_Reference_Point_Uncertainty()`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Point_Uncertainty.md),
[`Get_Reference_Points()`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Points.md),
[`Setup_Mod_Projection()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Projection.md),
[`catch_cap()`](https://chengmatt.github.io/SPoRC/dev/reference/catch_cap.md),
[`get_key_quants()`](https://chengmatt.github.io/SPoRC/dev/reference/get_key_quants.md),
[`hcr_threshold()`](https://chengmatt.github.io/SPoRC/dev/reference/hcr_threshold.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# advice may change by at most 15 percent a year, and never exceed 25 kt
Do_Population_Projection(..., F_constraint = list(catch_stability(down = 0.15, up = 0.15),
                                                  catch_cap(25)))
} # }
```
