# A mean and sd beta moved so its mode sits at a given value

The concentration `mu * (1 - mu) / sd^2 - 1` is kept and the mean and sd
are read back from the moved shape parameters, as
[`get_tagrep_prior`](https://chengmatt.github.io/SPoRC/dev/reference/get_tagrep_prior.md)
forms them. A beta whose concentration is at or below two has no
interior mode and is returned as given.

## Usage

``` r
beta_mode_at(mu, sd, x)
```

## Arguments

- mu, sd:

  The prior's mean and sd on the natural scale.

- x:

  Value in (0, 1) the mode is put at.

## Value

The moved `c(mu, sd)`.
