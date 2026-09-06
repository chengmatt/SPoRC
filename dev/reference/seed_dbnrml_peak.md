# Sensible starting values for the double normal's peak

The double normal holds the bin at which its plateau begins on the bin
scale, so a parameter left at zero puts the peak at bin zero. Where the
first bin is itself zero the ascending limb has no extent, its rescaling
divides by zero, and the curve comes back as `NaN`. A starting value in
the middle of the bin range gives a curve.

## Usage

``` r
seed_dbnrml_peak(pars, sel_model_arr, bin_vec, sex_offset = NULL)
```

## Arguments

- pars:

  Array `[region, par, block, sex, fleet]` of fixed-effect selectivity
  parameters.

- sel_model_arr:

  Integer array `[region, year, fleet]` of functional forms.

- bin_vec:

  Numeric vector of the bins selectivity is evaluated over.

- sex_offset:

  Character vector `[n_fleets]` of sex-offset specifications, used to
  tell a peak from an offset.

## Value

`pars` with the peak slot of every double normal fleet set to the middle
of the bin range, for the sexes that have a peak there.

## Details

Only a sex whose slot holds a peak is seeded. Under a par sex offset the
slots of every sex beyond the first hold offsets on the first sex's
parameters, where zero is the right starting value and the middle of the
bin range would shift that sex's curve by half the bin range.
