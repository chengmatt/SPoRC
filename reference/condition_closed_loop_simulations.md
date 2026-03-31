# Set up simulation list for closed-loop projections

This function creates and initializes a simulation list for closed-loop
projections of population dynamics. All components of the simulation
list must match the expected names used internally by the setup
functions. Users can provide \*\*custom definitions\*\* for any
component by passing them through \`...\` using the correct name (e.g.,
\`WAA_input\`, \`fish_sel_input\`).

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
  ISS_FishAgeComps_fill = "mean",
  ISS_FishLenComps_fill = "mean",
  ISS_SrvAgeComps_fill = "mean",
  ISS_SrvLenComps_fill = "mean",
  ...
)
```

## Arguments

- closed_loop_yrs:

  Integer. Number of years to project in the closed loop.

- n_sims:

  Integer. Number of simulation replicates.

- data:

  List. Observed data and configuration for the population.

- parameters:

  List. Parameter values for the model.

- mapping:

  List. Mapping of parameters for optimization.

- sd_rep:

  List. Standard deviation reports from fitted model.

- rep:

  List. Report from fitted model.

- random:

  Character vector of random effects estimated

- FishIdx_SE_fill:

  Character or numeric. Fill method for fishery index standard errors
  when extending to simulation years. Options are: - \`"zeros"\`: fill
  with zeros - \`"last"\`: repeat last non-NA slice - \`"mean"\`: fill
  with the mean of the observed series - Numeric: constant scalar or
  array value

- SrvIdx_SE_fill:

  Character or numeric. Fill method for survey index standard errors.
  Same options as \`FishIdx_SE_fill\`.

- ISS_FishAgeComps_fill:

  Character or numeric. Fill method for fishery age composition input
  sample sizes. Options are: - \`"zeros"\`, \`"last"\`, \`"mean"\` (as
  above) - \`"F_pattern"\`: fill based on fishing mortality pattern in
  the closed-loop simulation - Numeric: constant scalar or array value

- ISS_FishLenComps_fill:

  Character or numeric. Fill method for fishery length composition input
  sample sizes. Same options as \`ISS_FishAgeComps_fill\`.

- ISS_SrvAgeComps_fill:

  Character or numeric. Fill method for survey age composition input
  sample sizes. Options are \`"zeros"\`, \`"last"\`, \`"mean"\`, or a
  numeric constant.

- ISS_SrvLenComps_fill:

  Character or numeric. Fill method for survey length composition input
  sample sizes. Same options as \`ISS_SrvAgeComps_fill\`.

- ...:

  Optional named arguments for custom inputs. Each name must correspond
  to a component expected by the simulation setup functions, and be
  dimensioned appropriately: \`Setup_Sim_Fishing()\`,
  \`Setup_Sim_Survey()\`, \`Setup_Sim_Biologicals()\`,
  \`Setup_Sim_Rec()\`, and \`Setup_Sim_Tagging()\`. Examples include:

  \- \*\*Fishing\*\*: \`fish_sel_input\`, \`ln_sigmaC\`,
  \`Fmort_input\`, \`fish_q_input\`, etc. - \*\*Survey\*\*:
  \`srv_sel_input\`, \`srv_q_input\`, \`ObsSrvIdx_SE\`, etc. -
  \*\*Biologicals\*\*: \`WAA_input\`, \`MatAA_input\`,
  \`natmort_input\`, \`AgeingError_input\`, \`SizeAgeTrans_input\` -
  \*\*Recruitment inputs\*\*: - \`R0_input\`, \`h_input\`,
  \`sexratio_input\`, \`ln_InitDevs_input\`, \`Rec_input\` -
  \`Rec_input\`: - If shorter than the number of projection years, new
  recruitment deviates will be simulated based on \`recruitment_opt\`,
  \`R0_input\`, and \`h_input\` (supports changing regimes across
  years). - If you want fixed recruitment for all projection years,
  provide a \`Rec_input\` array that spans all years. - \*\*Tagging\*\*:
  \`Tag_Reporting_input\`, \`ln_Init_Tag_Mort\`, \`ln_Tag_Shed\`,
  \`tag_selex\`, \`tag_natmort\`, \`n_tags\`, \`n_tags_rel_input\` -
  \*\*Movement\*\*: \`Movement\` (must match the expected dimensions and
  be named exactly \`Movement\`)

  The values must have the correct dimensions expected by each
  component. If a component is not provided, default behavior will
  extend the last year (or zeros for fishing mortality, which can filled
  in subsequently)

## See also

Other Closed Loop Simulations:
[`catch_to_F_multifleet()`](https://chengmatt.github.io/SPoRC/reference/catch_to_F_multifleet.md),
[`catch_to_F_singlefleet()`](https://chengmatt.github.io/SPoRC/reference/catch_to_F_singlefleet.md),
[`get_closed_loop_reference_points()`](https://chengmatt.github.io/SPoRC/reference/get_closed_loop_reference_points.md)
