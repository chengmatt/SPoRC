# Cap The Catch A Projection Takes

Builds a constraint for `F_constraint` in
[`Do_Population_Projection`](https://chengmatt.github.io/SPoRC/dev/reference/Do_Population_Projection.md).
Where the retained catch the rule's F would take in a year is above the
cap, F is lowered until the catch equals it; elsewhere F is left as the
rule set it. The comparison is made on the plain numbers
([`ad_value`](https://chengmatt.github.io/SPoRC/dev/reference/ad_value.md)),
so the constraint also runs inside the model, which keeps the side of
the comparison that was true when it was built.
[`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md)
builds it at the estimates.

## Usage

``` r
catch_cap(cap, scope = "region")
```

## Arguments

- cap:

  Numeric. The cap on retained catch in biomass, in the units of
  `proj_Catch`. Under `scope = "region"`, one value for every region and
  year, a vector `[n_regions]`, or a matrix `[n_regions, n_proj_yrs]`
  indexed by the year the catch is taken, as `catch_input` is. Under
  `scope = "global"`, one value or a vector `[n_proj_yrs]`.

- scope:

  `"region"` (default) caps each region's own catch. `"global"` caps the
  catch summed over regions and shares the capped total across regions
  in proportion to the catch the rule's F would have taken, so the
  rule's regional pattern is kept. Pairs with
  `fmort_opt = "HCR_global"`, but runs under any rule.

## Value

A function `(F, state)` for `F_constraint`.

## See also

Other Reference Points and Projections:
[`Do_Population_Projection()`](https://chengmatt.github.io/SPoRC/dev/reference/Do_Population_Projection.md),
[`Get_Reference_Point_Uncertainty()`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Point_Uncertainty.md),
[`Get_Reference_Points()`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Points.md),
[`Setup_Mod_Projection()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Projection.md),
[`catch_stability()`](https://chengmatt.github.io/SPoRC/dev/reference/catch_stability.md),
[`get_key_quants()`](https://chengmatt.github.io/SPoRC/dev/reference/get_key_quants.md),
[`hcr_threshold()`](https://chengmatt.github.io/SPoRC/dev/reference/hcr_threshold.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# a 25 kt cap on every region and year
Do_Population_Projection(..., F_constraint = catch_cap(25))
# a cap on the total catch, compared against one global stock status
Do_Population_Projection(..., fmort_opt = "HCR_global",
                         F_constraint = catch_cap(60, scope = "global"))
} # }
```
