# Refuse a 2d logistic normal composition that is not joint by sex

The 2d logistic normal correlates bins within and across sexes, so it
reads one composition stacked over sexes. Aggregated or split by sex,
the operating model would cut the bins in half as if they were two
sexes, or fail inside the draw. The estimation model refuses the same
pairing.

## Usage

``` r
check_sim_2d_comp(comp_like, comp_type, like_name, type_name)
```

## Arguments

- comp_like:

  Likelihood code per fleet, 4 and 7 the 2d forms.

- comp_type:

  Composition type per year and fleet, or per fleet.

- like_name, type_name:

  Argument names, for the message.

## Value

`NULL`, invisibly. Called for the error it raises.
