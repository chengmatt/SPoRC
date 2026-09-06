# Construct and Condition Closed-Loop Simulation Inputs

Initializes and conditions a simulation object used for closed-loop
population projections. The function reconstructs fitted model
quantities, extends time-varying processes into projection years, and
prepares all simulation components required by the operating model.

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

  Integer. Number of projection years added beyond the fitted data
  period.

- n_sims:

  Integer. Number of stochastic simulation replicates.

- data:

  List. Data object used to fit the assessment model.

- parameters:

  List. Parameter vector from the fitted model.

- mapping:

  List. Parameter mapping object used during estimation.

- sd_rep:

  List. Standard deviation report from the fitted model.

- rep:

  List. Model report object produced by the fitted model.

- random:

  Character vector of estimated random effects.

- FishIdx_SE_fill:

  Character or numeric specifying how pooled fishery index standard
  errors are extended into projection years.

- SrvIdx_SE_fill:

  Character or numeric specifying how pooled survey index standard
  errors are extended into projection years.

- FishIdx_SE_pop_fill:

  Character or numeric specifying how population-specific fishery index
  standard errors are extended into projection years. Default `"mean"`.

- SrvIdx_SE_pop_fill:

  Character or numeric specifying how population-specific survey index
  standard errors are extended into projection years. Default `"mean"`.

- ISS_FishAgeComps_fill:

  Character or numeric specifying how pooled fishery age-composition
  input sample sizes are extended into projection years.

- ISS_FishLenComps_fill:

  Same behavior as `ISS_FishAgeComps_fill` for pooled fishery length
  compositions.

- ISS_FishAgeComps_discard_fill:

  Same behavior as `ISS_FishAgeComps_fill` for pooled fishery length
  compositions.

- ISS_FishLenComps_discard_fill:

  Same behavior as `ISS_FishAgeComps_fill` for pooled fishery length
  compositions.

- ISS_SrvAgeComps_fill:

  Character or numeric specifying how pooled survey age-composition
  input sample sizes are extended into projection years.

- ISS_SrvLenComps_fill:

  Same behavior as `ISS_SrvAgeComps_fill` for pooled survey length
  compositions.

- ISS_FishAgeComps_pop_fill:

  Character or numeric specifying how population-specific fishery
  age-composition input sample sizes are extended into projection years.
  Default `"mean"`.

- ISS_FishLenComps_pop_fill:

  Same behavior as `ISS_FishAgeComps_pop_fill` for population-specific
  fishery length compositions. Default `"mean"`.

- ISS_FishAgeComps_discard_pop_fill:

  Same behavior as `ISS_FishAgeComps_pop_fill` for population-specific
  fishery length compositions. Default `"mean"`.

- ISS_FishLenComps_discard_pop_fill:

  Same behavior as `ISS_FishLenComps_pop_fill` for population-specific
  fishery length compositions. Default `"mean"`.

  Extension rules for all `*_fill` arguments may be:

  `"zeros"`

  :   Fill projection years with zeros.

  `"last"`

  :   Repeat the final observed year.

  `"mean"`

  :   Use the mean of the historical series.

  `"F_pattern"`

  :   Scale fishery composition sample sizes proportionally to the
      simulated fishing mortality pattern (pooled fishery ISS only).

  numeric

  :   Use a supplied constant or array.

  If a fill argument receives an array directly, it is interpreted as
  the fully specified input and stored accordingly; the fill rule is
  ignored.

- ISS_SrvAgeComps_pop_fill:

  Character or numeric specifying how population-specific survey
  age-composition input sample sizes are extended into projection years.
  Default `"mean"`.

- ISS_SrvLenComps_pop_fill:

  Same behavior as `ISS_SrvAgeComps_pop_fill` for population-specific
  survey length compositions. Default `"mean"`.

- ...:

  Optional named simulation inputs that override internally generated
  components. Any argument expected by
  [`Setup_Sim_Fishing`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Fishing.md),
  [`Setup_Sim_Survey`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Survey.md),
  [`Setup_Sim_Biologicals`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Biologicals.md),
  [`Setup_Sim_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Rec.md),
  or
  [`Setup_Sim_Tagging`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Tagging.md)
  may be provided.

  Common overrides include:

  **Fishing processes**

  - `Fmort_input`

  - `fish_sel_input`

  - `fish_q_input`

  **Survey processes**

  - `srv_sel_input`

  - `srv_q_input`

  **Biological processes**

  - `WAA_input`

  - `MatAA_input`

  - `natmort_input`

  **Recruitment processes**

  - `R0_input`

  - `rinit_input`

  - `h_input`

  - `Rec_input`

  **Tagging processes**

  - `conv_tag_fish_reporting_input`

  **Movement processes**

  - `Movement`

  Supplied objects must have dimensions consistent with model structure,
  the total number of years (`length(data$years) + closed_loop_yrs`),
  and the number of simulations (`n_sims`).

## Value

A fully initialized `sim_list` object containing:

- model dimensions and simulation containers

- biological process inputs

- pooled and population-specific fishing and survey processes

- recruitment dynamics

- tagging processes

- spatial movement structures

- replicated arrays across `n_sims`

## Details

Simulation inputs are generated for biological dynamics, fishing,
surveys, recruitment, tagging, and movement processes. Historical years
are conditioned on outputs from the fitted assessment model, while
projection years are initialized using extension rules or user-supplied
overrides.

Users may replace any internally generated component by supplying a
correctly named object via `...`.

Simulation years consist of two periods:

- **Conditioning years**: historical years corresponding to the fitted
  assessment model.

- **Projection years**: future years simulated under closed-loop
  management.

During conditioning years, model processes are reconstructed from the
fitted model report objects. For projection years, quantities are
extended using the specified fill rules or user-supplied inputs.

By default:

- Biological inputs, selectivity, and catchability are extended using
  values from the final estimated year.

- Fishing mortality is initialized to zero in projection years.

- Recruitment is simulated forward when not fully specified by
  `Rec_input`.

- Population-specific data sources (`ObsFishIdx_pop_SE`,
  `ObsSrvIdx_pop_SE`, and all `*_pop` ISS arrays) fall back to
  uninformative defaults when the corresponding `Use*_pop` flags contain
  no ones.

Closed-loop feedback begins in the first projection year and allows
management actions (e.g., fishing mortality adjustments) to update
dynamically during the simulation.

## See also

Other Closed Loop Simulations:
[`catch_to_F_multifleet()`](https://chengmatt.github.io/SPoRC/dev/reference/catch_to_F_multifleet.md),
[`catch_to_F_singlefleet()`](https://chengmatt.github.io/SPoRC/dev/reference/catch_to_F_singlefleet.md),
[`get_closed_loop_reference_points()`](https://chengmatt.github.io/SPoRC/dev/reference/get_closed_loop_reference_points.md)
