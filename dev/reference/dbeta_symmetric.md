# Evaluate a symmetric beta log-density

Computes a log-density that penalizes a parameter value toward the
midpoint of `[p_lb, p_ub]` using a symmetric beta-like kernel. The
penalty strengthens as `p_prsd` increases and diffuses as `p_prsd`
decreases. Used in SPoRC as a prior for tag reporting rates to
discourage values near the boundaries of the unit interval.

## Usage

``` r
dbeta_symmetric(p_val, p_ub, p_lb, p_prsd, log = TRUE)
```

## Arguments

- p_val:

  Numeric. Parameter value to evaluate; must lie in `(p_lb, p_ub)`.

- p_ub:

  Numeric. Upper bound of the parameter support.

- p_lb:

  Numeric. Lower bound of the parameter support.

- p_prsd:

  Numeric. Pseudo-standard-deviation controlling prior concentration.
  Larger values produce a stronger penalty toward the midpoint
  \\(p\_{ub} + p\_{lb}) / 2\\; smaller values produce a more diffuse
  prior.

- log:

  Logical. If `TRUE` (default), returns the log-density; otherwise
  returns the density on the probability scale.

## Value

Numeric. Log-density (or density if `log = FALSE`).
