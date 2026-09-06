# Transform at-age observations onto the scale their likelihood is written on

A lognormal data source is fit on the log scale and a normal data source
on the natural scale, and the choice is per fleet, so the transformation
is applied cell by cell before
[`OBS`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html) registration.
Registration must happen against the name `getAll` supplied, so the
caller does it: a vector registered under a local name does not link to
the data element, and the objective then diverges from the reported
likelihood.

## Usage

``` r
prep_at_age_obs(obs, use, like_type, const = 0)
```

## Arguments

- obs:

  Observation array, shaped like `use`.

- use:

  Integer array flagging which cells are fit.

- like_type:

  Integer per fleet. `0` lognormal, `1` normal.

- const:

  Small constant added inside the log of a lognormal cell.

## Value

A numeric vector, one element per flagged cell in
[`which()`](https://rdrr.io/r/base/which.html) order, on the scale its
fleet's likelihood uses.
