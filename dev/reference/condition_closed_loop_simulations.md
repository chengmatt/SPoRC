# Construct and Condition Closed-Loop Simulation Inputs

Builds the simulation object a closed-loop projection runs on.
Historical years are conditioned on the fitted model's report, and
projection years are extended by the `*_fill` rules or by inputs
supplied through `...`, which replace any internally generated
component.

## Usage

``` r
condition_closed_loop_simulations(
  closed_loop_yrs,
  n_sims,
  data,
  parameters,
  mapping,
  sd_rep,
  rep,
  random = random,
  FishIdx_SE_fill = "mean",
  SrvIdx_SE_fill = "mean",
  FishIdx_SE_pop_fill = "mean",
  SrvIdx_SE_pop_fill = "mean",
  ISS_FishAgeComps_fill = "mean",
  ISS_FishLenComps_fill = "mean",
  ISS_FishAgeComps_discard_fill = "mean",
  ISS_FishLenComps_discard_fill = "mean",
  ISS_SrvAgeComps_fill = "mean",
  ISS_SrvLenComps_fill = "mean",
  ISS_FishAgeComps_pop_fill = "mean",
  ISS_FishLenComps_pop_fill = "mean",
  ISS_FishAgeComps_discard_pop_fill = "mean",
  ISS_FishLenComps_discard_pop_fill = "mean",
  ISS_SrvAgeComps_pop_fill = "mean",
  ISS_SrvLenComps_pop_fill = "mean",
  ...
)
```

## Arguments

- closed_loop_yrs:

  Integer projection years beyond the fitted data period.

- n_sims:

  Integer number of stochastic replicates.

- data:

  List. Data object used to fit the assessment model.

- parameters:

  List. Parameter vector from the fitted model.

- mapping:

  List. Parameter mapping used during estimation.

- sd_rep:

  List. Standard deviation report from the fitted model.

- rep:

  List. Report object from the fitted model.

- random:

  Character vector of estimated random effects.

- FishIdx_SE_fill, SrvIdx_SE_fill, FishIdx_SE_pop_fill,
  SrvIdx_SE_pop_fill:

  How the pooled and population-specific index standard errors are
  extended into the projection years. The population ones default to
  `"mean"`.

- ISS_FishAgeComps_fill, ISS_FishLenComps_fill, ISS_SrvAgeComps_fill,
  ISS_SrvLenComps_fill, ISS_FishAgeComps_discard_fill,
  ISS_FishLenComps_discard_fill:

  How the pooled input sample sizes are extended into the projection
  years.

- ISS_FishAgeComps_pop_fill, ISS_FishLenComps_pop_fill,
  ISS_SrvAgeComps_pop_fill, ISS_SrvLenComps_pop_fill,
  ISS_FishAgeComps_discard_pop_fill, ISS_FishLenComps_discard_pop_fill:

  The population-specific counterparts, each defaulting to `"mean"`.

  Every `*_fill` argument takes `"zeros"`, `"last"` (repeat the final
  observed year), `"mean"`, `"F_pattern"` (scale the sample sizes with
  the simulated fishing mortality pattern, pooled fishery input sample
  sizes only), or a constant. An array passed instead is taken as the
  fully specified input and the fill rule is ignored.

- ...:

  Optional named simulation inputs overriding what is generated
  internally. Any argument of
  [`Setup_Sim_Fishing`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Fishing.md),
  [`Setup_Sim_Survey`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Survey.md),
  [`Setup_Sim_Biologicals`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Biologicals.md),
  [`Setup_Sim_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Rec.md)
  or
  [`Setup_Sim_Tagging`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Tagging.md)
  may be given; the common ones are `Fmort_input`, `fish_sel_input`,
  `fish_q_input`, `srv_sel_input`, `srv_q_input`, `WAA_input`,
  `MatAA_input`, `natmort_input`, `R0_input`, `rinit_input`, `h_input`,
  `Rec_input`, `conv_tag_fish_reporting_input` and `Movement`.
  Dimensions must match the model structure,
  `length(data$years) + closed_loop_yrs` years and `n_sims` simulations.

## Value

A `sim_list` holding the model dimensions and simulation containers, the
biological inputs, the pooled and population-specific fishing and survey
processes, recruitment, tagging and movement, each replicated across
`n_sims`.

## Details

The conditioning years are the fitted model's own, reconstructed from
its report objects; the projection years are simulated under closed-loop
management and extended by the fill rules or the supplied inputs. By
default the biological inputs, selectivity and catchability are extended
at the final estimated year's values, fishing mortality starts at zero
in the projection years, recruitment is simulated forward unless
`Rec_input` specifies it, and the population-specific data sources fall
back to uninformative defaults when their `Use*_pop` flags hold no ones.
Feedback begins in the first projection year.

## See also

Other Closed Loop Simulations:
[`catch_to_F_multifleet()`](https://chengmatt.github.io/SPoRC/dev/reference/catch_to_F_multifleet.md),
[`catch_to_F_singlefleet()`](https://chengmatt.github.io/SPoRC/dev/reference/catch_to_F_singlefleet.md),
[`get_closed_loop_reference_points()`](https://chengmatt.github.io/SPoRC/dev/reference/get_closed_loop_reference_points.md)
