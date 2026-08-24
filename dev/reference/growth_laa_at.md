# Mean length and spread of every age at one point in a year

The size each integer age has reached at elapsed time `e` into year `y`.
Under curve growth (`cohort = 0`) every age is read from the year's
curve at its real age, with the plus group grown from its adjusted size.
Under cohort growth the ages that are propagated start from the
start-of-year state `L_beg` and grow by the year's increment; the ages
still in the linear phase take the length at `A1` their own cohort was
born with; and the first integer age past `A1` sits on the current
year's curve, which is where the propagation picks it up. The CV at age
is that of the year's own curve under curve growth and is held at the
first year's sizes under cohort growth.

## Usage

``` r
growth_laa_at(
  e,
  gp,
  ages,
  growth_A1,
  growth_A2,
  growth_L0,
  growth_cv_type,
  growth_sd_type,
  cohort,
  L_beg,
  L1_birth,
  cv_ref,
  a_prop,
  len_devs = NULL,
  growth_L2_asymptote = 0
)
```

## Arguments

- len_devs:

  Optional vector of log deviations on mean length at age, one per age,
  or `NULL` for a purely parametric curve.
