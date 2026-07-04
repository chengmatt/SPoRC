# Map fishery selectivity deviation parameters

Constructs the factor map for `ln_fishsel_devs`, the annual deviations
in continuous time-varying fishery selectivity. For iid and random-walk
forms, the deviation dimension corresponds to selectivity parameters (up
to 6 for double-normal); for semi-parametric forms (3D GMRF, 2D AR1), it
corresponds to age or length bins. Cells with no time-variation
(`cont_tv_fish_sel == 0`) or no catch data are mapped to `NA`.

## Usage

``` r
do_fishsel_devs_mapping(
  input_list,
  fish_sel_devs_spec,
  fishsel_devs_shared_bins,
  bins
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists.

- fish_sel_devs_spec:

  Character vector of length `n_fish_fleets`. Options:

  `"est_all"`

  :   Separate deviation series per region × sex × bin.

  `"est_shared_r"`

  :   Shared across regions.

  `"est_shared_s"`

  :   Shared across sexes.

  `"est_shared_r_s"`

  :   Shared across regions and sexes.

  `"est_shared_b"`

  :   Shared across bin groups defined by `fishsel_devs_shared_bins`.

  `"est_shared_r_b"`

  :   Shared across regions and bin groups.

  `"est_shared_b_s"`

  :   Shared across bin groups and sexes.

  `"est_shared_r_b_s"`

  :   Shared across regions, bin groups, and sexes.

  `"est_shared_f_x"`

  :   Copy deviation map from fleet `x`.

  `"fix"` or `"none"`

  :   All deviations fixed at zero (mapped to `NA`).

  Bin-sharing options (`"est_shared_b"`, etc.) are only valid with
  semi-parametric time-varying forms (`cont_tv_fish_sel` is in
  `c(3,4,5)`).

- fishsel_devs_shared_bins:

  List of integer vectors defining bin groupings for bin-sharing
  options. Each element groups bins that share a single deviation
  series, e.g., `list(1:5, 6:10, 11:30)`. Only used when
  `fish_sel_devs_spec` includes `"est_shared_b"`.

- bins:

  Number of selectivity bins

## Value

The input `input_list` with `$map$ln_fishsel_devs` set to a factor
vector and `$data$map_ln_fishsel_devs` set to the corresponding integer
array (for use in the RTMB template).

## Details

Bin-sharing (`"est_shared_b"` and related options) groups bins into
blocks defined by `fishsel_devs_shared_bins`, reducing the number of
estimated deviation series. Fleet sharing (`"est_shared_f_x"`) is
handled in a second pass. The resulting integer map is also stored as
`$data$map_ln_fishsel_devs` for use in the RTMB objective function.
