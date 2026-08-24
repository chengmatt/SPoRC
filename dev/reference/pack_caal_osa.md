# Pack observed CAAL data into a single flat OBS vector (OSA)

Produces the flat tracked OBS vector required by
[`RTMB::oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)
for conditional age-at-length data. Mirrors
[`pack_comp_osa`](https://chengmatt.github.io/SPoRC/dev/reference/pack_comp_osa.md)
with a length-bin loop inserted inside the season loop, so each length
bin becomes its own tracked group of age bins. Only the discrete
families are handled, since those are the only ones the CAAL likelihood
supports.

## Usage

``` r
pack_caal_osa(
  ObsArr,
  ISSArr,
  WtArr,
  UseArr,
  TypeMat,
  LikeTypeVec,
  n_yrs,
  n_seas,
  n_lens,
  n_fleets,
  n_sexes,
  addtocomp,
  return_labels = FALSE
)
```

## Arguments

- ObsArr:

  Observed CAAL array \\\[region \times year \times season \times len
  \times age \times sex \times fleet\]\\.

- ISSArr:

  Input sample sizes \\\[region \times year \times season \times len
  \times sex \times fleet\]\\.

- WtArr:

  Multinomial weights, same shape as `ISSArr`.

- UseArr:

  Use flags \\\[region \times year \times season \times len \times
  fleet\]\\.

- TypeMat:

  Composition type matrix \\\[year \times fleet\]\\.

- LikeTypeVec:

  Likelihood type per fleet.

- n_yrs:

  Number of model years.

- n_seas:

  Number of seasons per year.

- n_lens:

  Number of length bins.

- n_fleets:

  Number of fleets.

- n_sexes:

  Number of sexes.

- addtocomp:

  Small constant added to proportions before normalization.

- return_labels:

  Logical; if `TRUE`, also builds a per-element label data.frame
  identifying the origin of every entry, in the same order, for post-hoc
  relabeling of
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
  residuals.

## Value

If `return_labels = FALSE` (default), the flat OBS vector, or `NULL`
when no fleet carries CAAL data. If `return_labels = TRUE`, a list with
`vec` and `labels`.

## Details

Counts are formed the same way as for the marginal compositions:

- Multinomial (0): counts = round(prop x ISS x Wt)

- Dirichlet-multinomial (1): counts = round(prop x ISS)

The order is year, fleet, season, length bin, then region-fastest within
the group, and
[`eval_caal_osa`](https://chengmatt.github.io/SPoRC/dev/reference/eval_caal_osa.md)
walks the vector in exactly that order.
