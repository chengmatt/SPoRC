# Centering penalty on a set of selectivity fixed-effect parameters

Penalizes the squared log of the mean exponentiated value of a named set
of selectivity parameters, \\w \\ (\log \overline{e^{\theta}})^2\\. A
non-parametric selectivity curve is only identified up to a scalar once
catchability or fishing mortality is free to absorb it, and this pins
that scalar by pushing the set's average selectivity toward one, which
is a softer constraint than fixing a bin outright.

## Usage

``` r
get_selex_fixed_penalty(selex_penalty, fixed_sel_pars)
```

## Arguments

- selex_penalty:

  Data frame with columns `region`, `fleet`, `block`, `sex`, `par`, and
  `wt`, one row per penalized set. `par` is a list column of integer
  vectors naming the parameter indices in the set; `wt` is the weight.

- fixed_sel_pars:

  Array `[region, par, block, sex, fleet]` of selectivity fixed effects.

## Value

Numeric scalar negative log-likelihood contribution, summed across all
rows of `selex_penalty`.

## Details

Each row of the table names one set, so the penalty applies to the group
jointly rather than to each parameter separately. Because the expression
averages on the natural scale, it is meant for parameter sets held on
the log scale; a set stored on the logit scale (the non-parametric form,
or the asymptote of the asymptotic logistic forms) would not average to
anything interpretable as selectivity.
