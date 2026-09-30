# Gates on a linked recruitment series

One definition of each, read by the objective and by the operating
model, so the two cannot decide either differently. The arrows take over
an iid recruitment penalty: `dsem_rec_margvar_on` says that penalty's
marginal variance is the one in force, which is what the initial age
deviations are held at, and `dsem_rec_correction_on` says the lognormal
correction is applied on top. They part company where a model asks for
no correction but still takes its initial age spread from the arrows.

## Usage

``` r
dsem_rec_margvar_on(rec_linked, RecDevs_model)

dsem_rec_correction_on(
  rec_linked,
  RecDevs_model,
  RecDevs_pen_center,
  bias_correct_pe
)
```

## Arguments

- rec_linked:

  Whether any series is linked to `ln_RecDevs`.

- RecDevs_model:

  The deviations' own form, 1 for iid. A declared dsem stores 1, since
  the arrows stand in for that penalty.

- RecDevs_pen_center:

  1 when the penalty centers on the deviations' own mean, which leaves
  nothing to correct.

- bias_correct_pe:

  0 for `"none"`, which does the same.

## Value

`TRUE` when that gate is open.

## Details

The bias ramp has no part in either. A linked deviation is a random
effect, which takes the whole correction or none, so `Setup_Mod_Rec`
refuses a ramp alongside a declared dsem and `bias_correct_pe` decides
alone.
