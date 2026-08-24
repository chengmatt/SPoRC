# Selection-weighted weight at age

The mean weight of the fish a length-selective gear takes at each age,
\\\sum_l P(l \mid a) s(l) w(l) / \sum_l P(l \mid a) s(l)\\. This is the
weight a catch in biomass is made of when selectivity acts on length.

## Usage

``` r
get_selected_waa(key, sel_l, w_len)
```

## Arguments

- key:

  Matrix `[n_lens x n_ages]`, \\P(l \mid a)\\.

- sel_l:

  Selectivity at length, length `n_lens`.

- w_len:

  Weight at the bin midpoints, length `n_lens`.

## Value

Vector of length `n_ages`.
