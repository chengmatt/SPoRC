# Pack one composition data source's observations for OSA residuals

Builds the two flat observation vectors
[`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)
needs, one holding the fleets on a discrete likelihood and one the
fleets on a continuous one. Either is `NULL` when no fleet uses that
family.

## Usage

``` r
pack_comp_source_osa(
  ObsArr,
  ISSArr,
  WtArr,
  UseArr,
  TypeMat,
  LikeTypeVec,
  comp_bins_spec,
  n_yrs,
  n_seas,
  n_fleets,
  n_sexes,
  n_pop = 1,
  pop = FALSE,
  addtocomp = 0,
  do_internal_comp_osa = FALSE
)
```

## Arguments

- ObsArr:

  Observed compositions.

- ISSArr:

  Input sample size.

- WtArr:

  Likelihood weight applied to a multinomial.

- UseArr:

  Integer array flagging which cells are fit.

- TypeMat:

  Integer codes naming how each year and fleet's compositions are split,
  over year by fleet.

- LikeTypeVec:

  Integer likelihood per fleet.

- comp_bins_spec:

  Array of 0/1 flags over bin by fleet naming the bins this data source
  is fit over. All ones means the whole range.

- pop:

  Logical. `TRUE` for the population-specific data source, whose arrays
  have a leading population dimension and are never summed over
  populations.

- addtocomp:

  Small constant added to a composition.

- do_internal_comp_osa:

  Logical. `TRUE` hands the data source to
  [`eval_comp_source_osa`](https://chengmatt.github.io/SPoRC/dev/reference/eval_comp_source_osa.md),
  which reads the vectors `pack_comp_source_osa` built and the call site
  registered.

## Value

A list with `discrete` and `continuous`.

## Details

The result is registered with
[`OBS`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html) at the call
site rather than here. `OBS` takes an observation's name from the
variable it is called on, and
[`get_osa`](https://chengmatt.github.io/SPoRC/dev/reference/get_osa.md)
looks each data source up under that name, so registering inside this
function would give every data source the same one.
