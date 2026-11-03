# Threshold Harvest Control Rule

The target F above the trigger biomass `brp`, no fishing below `alpha`
times it, and a straight ramp between. Status is compared on its plain
numbers through
[`ad_value`](https://chengmatt.github.io/SPoRC/dev/reference/ad_value.md),
so the rule can run inside the model despite its `if`. The model keeps
the side of the comparison that was true when it was built, and
[`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md)
builds it at the estimates, so that side is the estimated stock's.

## Usage

``` r
hcr_threshold(x, frp, brp, alpha = 0.05)
```

## Arguments

- x:

  Spawning biomass the rule reads.

- frp:

  Fishing mortality reference point.

- brp:

  Biomass reference point.

- alpha:

  Share of `brp` below which F is zero. Default 0.05.

## Value

The F to apply. Inside the model it still depends on `x` and `frp`, so
`sdreport` accounts for their uncertainty.

## See also

Other Reference Points and Projections:
[`Do_Population_Projection()`](https://chengmatt.github.io/SPoRC/dev/reference/Do_Population_Projection.md),
[`Get_Reference_Point_Uncertainty()`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Point_Uncertainty.md),
[`Get_Reference_Points()`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Points.md),
[`Setup_Mod_Projection()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Projection.md),
[`catch_cap()`](https://chengmatt.github.io/SPoRC/dev/reference/catch_cap.md),
[`catch_stability()`](https://chengmatt.github.io/SPoRC/dev/reference/catch_stability.md),
[`get_key_quants()`](https://chengmatt.github.io/SPoRC/dev/reference/get_key_quants.md)

## Examples

``` r
hcr_threshold(x = 80, frp = 0.1, brp = 100) # on the ramp: 0.1 * (0.8 - 0.05) / 0.95
#> [1] 0.07894737
```
