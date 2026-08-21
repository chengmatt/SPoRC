# Set up sex offsets on selectivity for one selectivity stream

Parses the per-fleet sex-offset specification, stores the model flags,
and creates the curve scale-offset parameters with their factor map.
Under a `"par"` offset the sexes beyond the first store additive offsets
on the first sex's transformed-scale fixed-effect parameters, which
needs no new parameters, only the flag. Under a `"scale"` offset each
sex beyond the first carries a constant log-scale offset on its whole
realized curve, estimated per region and block.

## Usage

``` r
setup_sel_sex_offset(
  input_list,
  sex_offset,
  prefix,
  n_fleets,
  fleet_label,
  sel_model_arr,
  cont_tv_mat,
  max_blks,
  sel_blocks = NULL,
  fixed_spec = NULL,
  starting_values = list()
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map`.

- sex_offset:

  Character vector `[n_fleets]`: `"none"`, `"par"`, `"scale"`,
  `"par_scale"`, `"apical"`, or `"par_apical"`.

- prefix:

  One of `"fish"`, `"ret"`, or `"srv"`.

- n_fleets:

  Integer. Number of fleets in this stream.

- fleet_label:

  Character used in messages.

- sel_model_arr:

  Integer array `[region, year, fleet]` of functional forms, used to
  refuse a scale offset on forms whose standardization would cancel it
  (non-parametric forms 5 and 9).

- cont_tv_mat:

  Integer matrix `[region, fleet]` of continuous time-varying
  structures, used to refuse a scale offset under the semi-parametric
  structures (3-5) for the same reason.

- max_blks:

  Integer. Maximum number of selectivity blocks, sizing the scale
  parameter array.

- sel_blocks:

  Integer array `[region, year, fleet]` of selectivity block indices, so
  a scale offset is only estimated for the blocks a fleet actually has
  (the array is padded to `max_blks`).

- fixed_spec:

  Character vector `[n_fleets]` of fixed-parameter sharing
  specifications, or `NULL`. A `"par"` offset reads the later sexes'
  slots as offsets on the first sex's, so a specification that shares
  those slots across sexes (`"est_shared_s"`, `"est_shared_r_s"`) would
  silently double the first sex's parameters and is refused.

- starting_values:

  Named list of user-supplied starting values.

## Value

`input_list` with the flags, scale parameters, and map added.
