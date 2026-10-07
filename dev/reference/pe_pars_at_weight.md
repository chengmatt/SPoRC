# Process error parameters a weighted penalty draws at

A penalty multiplied by `wt` is, up to a constant, the density of the
same process with its variance divided by `wt`: each log sd is lowered
by `log(wt) / 2` under iid, a random walk or the separable AR1, and the
3D GMRF's log variance by `log(wt)`. The correlations are unchanged.

## Usage

``` r
pe_pars_at_weight(pe_pars, PE_model, wt)
```

## Arguments

- pe_pars:

  Array `[1, 1, slot, sex]` for one region and fleet.

- PE_model:

  Integer process code, as `Get_PE_loglik` reads it.

- wt:

  Penalty weight.

## Value

`pe_pars` at that weight.
