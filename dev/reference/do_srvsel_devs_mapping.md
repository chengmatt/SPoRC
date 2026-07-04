# Map survey selectivity deviation parameters

Internal helper called by
[`Setup_Mod_Srvsel_and_Q`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Srvsel_and_Q.md)
to construct the TMB/RTMB factor map for `ln_srvsel_devs`
`[n_regions × (n_years + n_proj_yrs_devs) × n_bins × n_sexes × n_srv_fleets]`,
where `n_bins` is `n_ages` for age-based selectivity or `n_lens` for
length-based selectivity. For iid and random walk forms, the bin
dimension indexes selectivity parameters (1, 2, or 6 depending on
functional form); for semi-parametric forms (3D GMRF, 2D AR1), it
indexes age or length bins directly.

## Usage

``` r
do_srvsel_devs_mapping(
  input_list,
  srv_sel_devs_spec,
  srvsel_devs_shared_bins,
  bins
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists. Requires
  `$data$n_srv_fleets`, `$data$n_regions`, `$data$n_sexes`,
  `$data$cont_tv_srv_sel`, `$data$srv_sel_model`, `$data$UseSrvIdx`,
  `$data$UseSrvIdx_pop`, `$data$n_proj_yrs_devs`, and
  `$data$srv_selex_type`.

- srv_sel_devs_spec:

  Character vector `[n_srv_fleets]` specifying the sharing structure for
  selectivity deviations. One of:

  `"est_all"`

  :   Separate deviation time series per region, sex, bin, and fleet.

  `"est_shared_r"`

  :   Shared across regions; separate by sex and bin.

  `"est_shared_s"`

  :   Shared across sexes; separate by region and bin.

  `"est_shared_r_s"`

  :   Shared across regions and sexes.

  `"est_shared_b"`

  :   Shared across bin groups defined by `srvsel_devs_shared_bins`.
      Semi-parametric forms only.

  `"est_shared_r_b"`

  :   Shared across regions and bin groups. Semi-parametric forms only.

  `"est_shared_b_s"`

  :   Shared across bin groups and sexes. Semi-parametric forms only.

  `"est_shared_r_b_s"`

  :   Shared across regions, bin groups, and sexes. Semi-parametric
      forms only.

  `"est_shared_f_x"`

  :   Copy the map from fleet `x`. Fleet `x` must not itself use
      `"est_shared_f_x"`.

  `"fix"`/`"none"`

  :   All deviations fixed at zero (mapped to `NA`).

- srvsel_devs_shared_bins:

  List of integer vectors defining bin groups for age/length sharing.
  Each element specifies the bin indices within one group (e.g.,
  `list(1:5, 6:10, 11:30)`). Required when `srv_sel_devs_spec` includes
  `"est_shared_b"` or its variants; ignored otherwise.

- bins:

  Number of seletivity bins

## Value

The input `input_list` with `$map$ln_srvsel_devs` set to a factor vector
and `$data$map_ln_srvsel_devs` set to the equivalent integer array
`[n_regions × (n_years + n_proj_yrs_devs) × n_bins × n_sexes × n_srv_fleets]`.
Active parameters receive sequential integer indices; inactive
parameters are `NA`.

## Details

Bin-sharing options (`"est_shared_b"` and variants) are only valid for
semi-parametric forms and require `srvsel_devs_shared_bins` to define
which bin groups share a common deviation. Fleet-sharing
(`"est_shared_f_x"`) is handled in a second pass. Parameters are
automatically fixed for fleets with no continuous time-variation or no
active survey index data. The integer-valued map is also stored in
`$data$map_ln_srvsel_devs` for use within the TMB objective function.
