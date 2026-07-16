# Defining Priors

In some cases, prior information might be available for use to inform
model estimation. Priors can either serve as informative priors that
have a strong influence on model estimation (e.g., natural mortality,
catchability), but can also serve as regularizing priors that help
stabilize model estimation. These regularizing priors can be helpful
particularly for stabilizing tag reporting, movement, and selectivity
parameters because some parameters are not well informed by data (e.g.,
slope of selectivity when it approaches infinity) such that the gradient
is zero with respect to the likelihood surface.

In `SPoRC`, priors can be placed on several parameters. These include:

1.  Steepness,
2.  Natural Mortality,
3.  Catchability,
4.  Selectivity,
5.  Movement,
6.  Tag Reporting Rates,
7.  Recruitment Regional, Seasonal Proportions, and R0,
8.  Stray rates

In the following, we will demonstrate how priors can be specified for
each of these parameters. For further mathematical details on the
formulation of these priors, refer to the vignette describing model
equations.

Let us first load in any required packages, data we may use, and define
model dimensions.

``` r

library(SPoRC) 
#> Loading required package: RTMB
data("sgl_rg_sable_data")

input_list <- Setup_Mod_Dim(years = 1:length(sgl_rg_sable_data$years),
                            ages = 1:length(sgl_rg_sable_data$ages),
                            lens = seq(41, 99, 2),
                            n_regions = 1,
                            n_sexes = sgl_rg_sable_data$n_sexes,
                            n_fish_fleets = sgl_rg_sable_data$n_fish_fleets,
                            n_srv_fleets = sgl_rg_sable_data$n_srv_fleets,
                            n_seas = 1,
                            n_pop = sgl_rg_sable_data$n_pop,
                            verbose = TRUE)
#> Number of Years: 65
#> Number of Seasons: 1
#> Duration of season 1: 1
#> Number of Projection Years for Dev Pars: 0
#> Number of Regions: 1
#> Number of Populations: 1
#> Number of Age Bins: 30
#> Number of Length Bins: 30
#> Number of Sexes: 2
#> Number of Fishery Fleets: 2
#> Number of Survey Fleets: 3
```

## Steepness

To specify priors on the steepness of a Beverton-Holt stock-recruitment
relationship, we use the `Setup_Mod_Rec` function. First, we define a
dataframe that includes the population (`pop`), region, mean (`mu`), and
standard deviation (`sd`) of the steepness prior. These values are
specified in normal space but are internally transformed to a beta
distribution bounded between 0.2 and 1. To activate steepness priors,
set `Use_h_prior = 1` and provide the dataframe to the `h_prior`
argument:

``` r

steepness_prior <- data.frame(
  pop    = 1,   # population 1
  region = 1,   # region 1
  mu     = 0.6, # prior mean
  sd     = 0.2  # prior sd
)

input_list <- Setup_Mod_Rec(
  input_list = input_list,
  rec_model = "bh_rec",
  Use_h_prior = 1,
  h_prior = steepness_prior,
  h_spec = "est_shared_pop_r"
)
#> Recruitment is specified as: bh_rec
#> Recruitment Density Dependence is specified as: global
#> Recruitment and SSB lag is specified as: 1
#> Recruitment regional proportion priors are: Not Used
#> Recruitment seasonal proportion priors are: Not Used
#> Recruitment seasonal proportions is: fix
#> Stray rates for population 1 specified with 1 block(s).
#> Stray rate prior is: Not Used
#> Steepness priors are: Used
#> Sex Ratios specified with 1 block for population 1 and region 1
#> Recruitment Bias Ramp is: Off
#> Initial Age Structure is: Movement and Matrix Geometric Series
#> Recruitment deviations for every year are estimated
#> Spawning season occurs in season 1
#> Recruitment Variability is specified as: est_all
#> Initial Age Deviations is stochastic for all ages, but the plus group follows equilibrium calculations.
#> Recruitment Deviations is estimated for all dimensions
#> Steepness is specified as: est_shared_pop_r
#> Sex ratio is specified as: fix
#> Stray rates fixed (n_pop == 1, straying not applicable).
```

In a spatial model where region-specific steepness priors are required,
the prior dataframe can include multiple entries. Each row corresponds
to a population–region combination along with its associated mean and
standard deviation:

``` r

steepness_prior <- data.frame(
  pop    = 1,
  region = c(1, 2, 3),
  mu     = c(0.6, 0.8, 0.3),
  sd     = c(0.2, 0.2, 0.2)
)

steepness_prior
#>   pop region  mu  sd
#> 1   1      1 0.6 0.2
#> 2   1      2 0.8 0.2
#> 3   1      3 0.3 0.2
```

## Natural Mortality

To define natural mortality priors, we use the `Setup_Mod_Biologicals`
function. Setting `Use_M_prior = 1` activates the prior, and a dataframe
of values is provided via the `M_prior` argument. Each row corresponds
to a specific combination of population, region, year, age, and sex
blocks, along with their associated mean and standard deviation.

``` r

M_prior <- data.frame( 
  popblk    = 1,
  regionblk = 1,
  yearblk   = 1,
  ageblk    = 1,
  sexblk    = c(1, 2),
  mu        = 0.085,
  sd        = 0.1
)

input_list <- Setup_Mod_Biologicals(
  input_list = input_list,
  WAA          = sgl_rg_sable_data$WAA,
  MatAA        = sgl_rg_sable_data$MatAA,
  AgeingError  = as.matrix(sgl_rg_sable_data$age_error),
  SizeAgeTrans = sgl_rg_sable_data$SizeAgeTrans,
  Use_M_prior  = 1,
  M_prior      = M_prior,
  M_spec       = "est_ln_M",
  M_sexblk_spec = list(1, 2)
)
#> Length Composition data are: Not Used
#> Natural Mortality priors are: Used
#> WAA_fish was specified at NULL. Using the spawning WAA for WAA_fish
#> WAA_srv was specified at NULL. Using the spawning WAA for WAA_srv
#> Ageing Error is specified to be time-invariant
#> Natural Mortality specified as: est_ln_M
#> Natural Mortality Population Blocks is specified as: 1
#> Natural Mortality Region Blocks is specified as: 1
#> Natural Mortality Year Blocks is specified as: 1
#> Natural Mortality Age Blocks is specified as: 1
#> Natural Mortality Sex Blocks is specified as: 2
```

## Catchability

Defining catchability priors for the fishery and survey follows the same
process. For demonstration purposes, the following section will
illustrate how catchability priors can be defined for the survey using
the `Setup_Mod_Srvsel_and_Q` function. A dataframe of catchability
priors is first defined that determines the region, fleet, and block to
which the prior is applied. The dataframe also includes columns defining
the mean (`mu`) and standard deviation (`sd`) of the catchability prior.
These priors are provided in normal space, but are internally
transformed to a lognormal density. Setting `Use_srv_q_prior = 1`
activates the prior.

``` r

srv_q_prior <- data.frame(
  region = 1,
  fleet  = c(1, 2, 3),
  block  = 1,
  mu     = c(6, 0.5, 4),
  sd     = c(0.1, 0.35, 0.1)
)

srv_q_prior
#>   region fleet block  mu   sd
#> 1      1     1     1 6.0 0.10
#> 2      1     2     1 0.5 0.35
#> 3      1     3     1 4.0 0.10

input_list <- Setup_Mod_Srvsel_and_Q(
  input_list = input_list,
  cont_tv_srv_sel        = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3"),
  srv_sel_blocks         = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3"),
  srv_sel_model          = c("logist1_Fleet_1", "exponential_Fleet_2", "logist1_Fleet_3"),
  srv_q_blocks           = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3"),
  srv_fixed_sel_pars_spec = c("est_all", "est_all", "est_all"),
  srv_q_spec             = c("est_all", "est_all", "est_all"),
  Use_srv_q_prior        = 1,
  srv_q_prior            = srv_q_prior
)
#> Survey Catchability priors are: Used
#> Survey Selectivity priors are: Not Used
#> Survey Selectivity is aged-based.
#> Continuous survey time-varying selectivity specified as: none for survey fleet 1
#> Continuous survey time-varying selectivity specified as: none for survey fleet 2
#> Continuous survey time-varying selectivity specified as: none for survey fleet 3
#> Survey Selectivity Time Blocks for survey 1 is specified at: 1
#> Survey Selectivity Time Blocks for survey 2 is specified at: 1
#> Survey Selectivity Time Blocks for survey 3 is specified at: 1
#> Survey selectivity functional form specified as:logist1 for survey fleet 1
#> Survey selectivity functional form specified as:exponential for survey fleet 2
#> Survey selectivity functional form specified as:logist1 for survey fleet 3
#> Survey Catchability Time Blocks for survey 1 is specified at: 1
#> Survey Catchability Time Blocks for survey 2 is specified at: 1
#> Survey Catchability Time Blocks for survey 3 is specified at: 1
#> srv_fixed_sel_pars_spec is specified as: est_all for survey fleet 1
#> srv_fixed_sel_pars_spec is specified as: est_all for survey fleet 2
#> srv_fixed_sel_pars_spec is specified as: est_all for survey fleet 3
#> srv_q_spec is specified as: est_all for survey fleet 1
#> srv_q_spec is specified as: est_all for survey fleet 2
#> srv_q_spec is specified as: est_all for survey fleet 3
```

## Selectivity

Defining selectivity priors for the fishery and survey follows the same
process. Here, we will focus on defining selectivity priors for the
survey for brevity, using the `Setup_Mod_Srvsel_and_Q` function. We
first define a dataframe with the following columns: `region`, `par`,
`block`, `sex`, `fleet`, `mu`, and `sd`, which defines the region,
parameter, time block, sex, and fleet combination to apply the prior to.
Given the complexity of the number of selectivity parameters that can be
estimated across model partitions (regions, sexes, blocks, and fleets),
users will need to be cautious about how the selectivity prior dataframe
is defined to ensure they are correctly applied.

In the following, we will apply priors to the logistic selectivity
parameters for both sexes in fleet 1, and priors to the power parameter
describing exponential selectivity in fleet 2.

``` r

srv_selex_prior <- data.frame(
  region = 1,
  par    = c(1, 2, 1, 2, 1, 1),
  sex    = c(1, 1, 2, 2, 1, 2),
  fleet  = c(1, 1, 1, 1, 2, 2),
  block  = 1,
  mu     = rep(0, 6),
  sd     = rep(5, 6)
)

srv_selex_prior
#>   region par sex fleet block mu sd
#> 1      1   1   1     1     1  0  5
#> 2      1   2   1     1     1  0  5
#> 3      1   1   2     1     1  0  5
#> 4      1   2   2     1     1  0  5
#> 5      1   1   1     2     1  0  5
#> 6      1   1   2     2     1  0  5
```

In the example provided above, the first row of the dataframe indicates
that a prior is applied to region 1, parameter 1, sex 1, fleet 1, and
block 1, which corresponds to the $`b_{50}`$ parameter for logistic
selectivity specified for the first survey fleet. The last row indicates
that a prior is applied to region 1, parameter 1, sex 2, fleet 2, and
block 1, which corresponds to the male $`\phi`$ power parameter for
exponential selectivity specified for the second survey fleet.

Next, we can activate selectivity priors by setting
`Use_srv_selex_prior = 1` and providing the associated selectivity prior
dataframe to `srv_selex_prior`.

``` r

input_list <- Setup_Mod_Srvsel_and_Q(
  input_list = input_list,
  cont_tv_srv_sel         = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3"),
  srv_sel_blocks          = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3"),
  srv_sel_model           = c("logist1_Fleet_1", "exponential_Fleet_2", "logist1_Fleet_3"),
  srv_q_blocks            = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3"),
  srv_fixed_sel_pars_spec = c("est_all", "est_all", "est_all"),
  srv_q_spec              = c("est_all", "est_all", "est_all"),
  Use_srv_selex_prior     = 1,
  srv_selex_prior         = srv_selex_prior
)
#> Survey Catchability priors are: Not Used
#> Survey Selectivity priors are: Used
#> Survey Selectivity is aged-based.
#> Continuous survey time-varying selectivity specified as: none for survey fleet 1
#> Continuous survey time-varying selectivity specified as: none for survey fleet 2
#> Continuous survey time-varying selectivity specified as: none for survey fleet 3
#> Survey Selectivity Time Blocks for survey 1 is specified at: 1
#> Survey Selectivity Time Blocks for survey 2 is specified at: 1
#> Survey Selectivity Time Blocks for survey 3 is specified at: 1
#> Survey selectivity functional form specified as:logist1 for survey fleet 1
#> Survey selectivity functional form specified as:exponential for survey fleet 2
#> Survey selectivity functional form specified as:logist1 for survey fleet 3
#> Survey Catchability Time Blocks for survey 1 is specified at: 1
#> Survey Catchability Time Blocks for survey 2 is specified at: 1
#> Survey Catchability Time Blocks for survey 3 is specified at: 1
#> srv_fixed_sel_pars_spec is specified as: est_all for survey fleet 1
#> srv_fixed_sel_pars_spec is specified as: est_all for survey fleet 2
#> srv_fixed_sel_pars_spec is specified as: est_all for survey fleet 3
#> srv_q_spec is specified as: est_all for survey fleet 1
#> srv_q_spec is specified as: est_all for survey fleet 2
#> srv_q_spec is specified as: est_all for survey fleet 3
```

## Movement

Defining priors for movement, tag reporting rates, and recruitment
proportions (see sections below) is best illustrated with a spatial
model. Here, we will load in data from the multi-region sablefish case
study.

``` r

data(mlt_rg_sable_data)

input_list <- Setup_Mod_Dim(years = 1:length(mlt_rg_sable_data$years),
                            ages = 1:length(mlt_rg_sable_data$ages),
                            lens = mlt_rg_sable_data$lens,
                            n_regions = mlt_rg_sable_data$n_regions,
                            n_sexes = mlt_rg_sable_data$n_sexes,
                            n_fish_fleets = mlt_rg_sable_data$n_fish_fleets,
                            n_srv_fleets = mlt_rg_sable_data$n_srv_fleets,
                            n_seas = 1,
                            n_pop = mlt_rg_sable_data$n_pop,
                            verbose = TRUE)
#> Number of Years: 62
#> Number of Seasons: 1
#> Duration of season 1: 1
#> Number of Projection Years for Dev Pars: 0
#> Number of Regions: 5
#> Number of Populations: 1
#> Number of Age Bins: 30
#> Number of Length Bins: 30
#> Number of Sexes: 2
#> Number of Fishery Fleets: 2
#> Number of Survey Fleets: 2
```

Because movement parameters are correlated and must sum to 1, `SPoRC`
uses a Dirichlet prior to model movement probabilities. Movement priors
are enabled by setting `Use_Movement_Prior = 1`. The `Movement_prior`
argument must be provided as a dataframe where each row corresponds to a
combination of population, origin region, year, age, season, and sex
indices, and includes a list-column `alpha` specifying the Dirichlet
concentration parameters that define the prior over movement
probabilities to all destination regions.

``` r

prior <- expand.grid(
  pop         = 1,
  region_from = 1:input_list$data$n_regions,
  year        = 1,
  age         = 1,
  seas        = 1,
  sex         = 1,
  alpha       = I(list(rep(2, input_list$data$n_regions)))
)

input_list <- Setup_Mod_Movement(
  input_list           = input_list,
  Movement_ageblk_spec = "constant",
  Movement_yearblk_spec = "constant",
  Movement_sexblk_spec = "constant",
  Use_Movement_Prior   = 1,
  Movement_prior       = prior
)
#> Movement is: Estimated
#> Movement priors are: Used
#> Recruits are: Not Moving
#> Continuous movement specification is: none
#> Continuous movement process error specification is: none
#> Movement type is: Unstructured Markov
#> Movement fixed effect blocks are population-invariant
#> Movement fixed effect blocks are season-invariant
#> Movement fixed effect blocks are sex-invariant
#> Movement fixed effect blocks are time-invariant
#> Movement fixed effect blocks are age-invariant
```

## Tag Reporting Rates

Defining priors for tag reporting rates is done by constructing a
dataframe and passing it to the `Setup_Mod_Tagging` function. The
dataframe must contain the following columns: `region`, `block`,
`fleet`, `mu`, `sd`, and `type`. These columns define the region and
time block the prior applies to, the mean and standard deviation of the
prior distribution, and the type of beta prior to use. Two prior types
are supported: a symmetric beta prior (`type = 0`), which penalizes
extreme values near 0 and 1, and a standard beta prior (`type = 1`),
which allows users to specify both the mean and standard deviation
directly.

In the example below, symmetric beta priors are applied to two time
blocks in region 1. Because the symmetric form does not require a mean,
the `mu` column is set to `NA`:

``` r

tag_prior <- data.frame(
  region = 1,
  block  = c(1, 2),
  fleet  = 1,
  mu     = NA,
  sd     = 5,
  type   = 0
)
tag_prior
#>   region block fleet mu sd type
#> 1      1     1     1 NA  5    0
#> 2      1     2     1 NA  5    0
```

The first row of the dataframe indicates that a symmetric beta prior is
applied to tag reporting rates in region 1 and block 1, with a standard
deviation of 5. The second row applies the same prior structure to block
2. Because `type = 0`, the prior is symmetric around 0.5 and discourages
extreme reporting rate estimates. If `type = 1` were used instead, the
`mu` column would need to be specified with the prior mean.

Tag reporting rate priors are activated by setting
`use_conv_tag_fishrep_prior = 1` and providing the dataframe to the
`conv_tag_fishrep_prior` argument.

``` r

input_list <- Setup_Mod_Tagging(
  input_list = input_list,
  use_conv_fish_tagging      = c(1, 0),
  conv_tag_max_liberty       = 15,

  # Data inputs
  conv_tag_release_indicator = mlt_rg_sable_data$conv_tag_release_indicator,
  conv_tagged_fish           = mlt_rg_sable_data$conv_tagged_fish,
  obs_conv_tag_fish_recap    = mlt_rg_sable_data$obs_conv_tag_fish_recap,

  # Model options
  conv_fish_tag_like           = "Poisson",
  use_conv_tag_fishrep_prior   = 1,
  conv_tag_fishrep_prior       = tag_prior,
  init_conv_tag_mort_spec      = "fix",
  conv_tag_shed_spec           = "fix",
  conv_tagrep_spec             = "est_shared_r_f",
  conv_tag_fish_reporting_blocks = c(
    apply(expand.grid(1:input_list$data$n_regions, 1:input_list$data$n_fish_fleets), 1, function(x)
      paste0("Block_1_Year_1-35_Region_", x[1], "_Fleet_", x[2])),
    apply(expand.grid(1:input_list$data$n_regions, 1:input_list$data$n_fish_fleets), 1, function(x)
      paste0("Block_2_Year_36-terminal_Region_", x[1], "_Fleet_", x[2]))
  ),

  # Starting / fixed values
  ln_init_conv_tag_mort = log(0.1),
  ln_conv_tag_shed      = log(0.02)
)
#> Tagging priors are used
#> Conventional Tag Likelihood specified as: Poisson
#> Conventional Tagging data are fit to 1 population groups
#> Conventional Tagging data are fit to 30 age groups
#> Conventional Tagging data are fit to 2 sex groups
#> Conventional Tag Reporting estimated with 2 blocks for region 1 and fleet 1
#> Conventional Tag Reporting estimated with 2 blocks for region 2 and fleet 1
#> Conventional Tag Reporting estimated with 2 blocks for region 3 and fleet 1
#> Conventional Tag Reporting estimated with 2 blocks for region 4 and fleet 1
#> Conventional Tag Reporting estimated with 2 blocks for region 5 and fleet 1
#> Conventional Tag Reporting estimated with 2 blocks for region 1 and fleet 2
#> Conventional Tag Reporting estimated with 2 blocks for region 2 and fleet 2
#> Conventional Tag Reporting estimated with 2 blocks for region 3 and fleet 2
#> Conventional Tag Reporting estimated with 2 blocks for region 4 and fleet 2
#> Conventional Tag Reporting estimated with 2 blocks for region 5 and fleet 2
#> Conventional Initial Tag Mortality is specified as: fix
#> Conventional Chronic Tag Shedding is specified as: fix
#> Conventional Tag Reporting is specified as: est_shared_r_f
```

## Recruitment Proportions and R0

### Recruitment Regional Apportionment

In the context of a spatial model, Dirichlet priors can also be placed
on recruitment proportions. This is done through the `Setup_Mod_Rec`
function by setting `use_rec_region_prop_prior = 1` and providing a
dataframe to `rec_region_prop_prior`. Each row of the dataframe
corresponds to a population and contains a list-column `alpha`
specifying the Dirichlet concentration parameters for movement among
destination regions. If a uniform Dirichlet prior is desired, all
concentration parameters can be set to the same value. Providing
different values for each region will apply a weighted Dirichlet prior.

``` r

input_list <- Setup_Mod_Rec(
  input_list = input_list,
  do_rec_bias_ramp = 0,
  rec_model = "mean_rec",
  use_rec_region_prop_prior = 1,
  rec_region_prop_prior = data.frame(
    pop   = 1,
    alpha = I(list(rep(5, mlt_rg_sable_data$n_regions)))
  )
)
#> Recruitment is specified as: mean_rec
#> Recruitment Density Dependence is specified as: global
#> Recruitment regional proportion priors are: Used
#> Recruitment seasonal proportion priors are: Not Used
#> Recruitment seasonal proportions is: fix
#> Stray rates for population 1 specified with 1 block(s).
#> Stray rate prior is: Not Used
#> Sex Ratios specified with 1 block for population 1 and region 1
#> Sex Ratios specified with 1 block for population 1 and region 2
#> Sex Ratios specified with 1 block for population 1 and region 3
#> Sex Ratios specified with 1 block for population 1 and region 4
#> Sex Ratios specified with 1 block for population 1 and region 5
#> Recruitment Bias Ramp is: Off
#> Initial Age Structure is: Movement and Matrix Geometric Series
#> Recruitment deviations for every year are estimated
#> Spawning season occurs in season 1
#> Recruitment Variability is specified as: est_all
#> Initial Age Deviations is stochastic for all ages, but the plus group follows equilibrium calculations.
#> Recruitment Deviations is estimated for all dimensions
#> Sex ratio is specified as: fix
#> Stray rates fixed (n_pop == 1, straying not applicable).
```

### Recruitment Seasonal Apportionment

Likewise, priors on recruitment seasonal apportionment can be applied in
an analogous manner. This is done by setting
`use_rec_seas_prop_prior = 1` and providing a dataframe to
`rec_seas_prop_prior`. Note that seasonal recruitment proportion priors
are only relevant when `n_seas > 1` and `use_fixed_rec_seas_prop = 0`
(i.e., seasonal proportions are being estimated rather than fixed). Each
row of the dataframe corresponds to a population and contains a
list-column `alpha` specifying the Dirichlet concentration parameters
for the seasonal apportionment of recruitment.

### R0

Lognormal priors on $`R_0`$ can be applied per population via
`use_r0_prior` and `r0_prior`. Each row of `r0_prior` corresponds to a
population, with `mu` specifying the prior mean on the natural scale and
`sd` the standard deviation on the log scale.

``` r

input_list <- Setup_Mod_Rec(
  input_list = input_list,
  use_r0_prior = 1,
  do_rec_bias_ramp = 0,
  rec_model = "mean_rec",
  r0_prior = data.frame(
    pop = 1,
    mu  = 1e6,
    sd  = 1
  )
)
#> Recruitment is specified as: mean_rec
#> Recruitment Density Dependence is specified as: global
#> Recruitment regional proportion priors are: Not Used
#> Recruitment seasonal proportion priors are: Not Used
#> Recruitment seasonal proportions is: fix
#> Stray rates for population 1 specified with 1 block(s).
#> Stray rate prior is: Not Used
#> Sex Ratios specified with 1 block for population 1 and region 1
#> Sex Ratios specified with 1 block for population 1 and region 2
#> Sex Ratios specified with 1 block for population 1 and region 3
#> Sex Ratios specified with 1 block for population 1 and region 4
#> Sex Ratios specified with 1 block for population 1 and region 5
#> Recruitment Bias Ramp is: Off
#> Initial Age Structure is: Movement and Matrix Geometric Series
#> Recruitment deviations for every year are estimated
#> Spawning season occurs in season 1
#> Recruitment Variability is specified as: est_all
#> Initial Age Deviations is stochastic for all ages, but the plus group follows equilibrium calculations.
#> Recruitment Deviations is estimated for all dimensions
#> Sex ratio is specified as: fix
#> Stray rates fixed (n_pop == 1, straying not applicable).
```

### Stray Rates

When stray rates are estimated (`stray_rate_spec` is not `"fix"`), a
beta prior should almost always be used to keep estimation stable.
Estimation of stray rates are only applicable when `n_pop > 1`. The
prior is specified via a dataframe passed to `stray_rate_prior`, with
columns `pop`, `block`, `mu`, and `sd`. The mean (`mu`) and standard
deviation (`sd`) are supplied in probability space and are internally
converted to beta shape parameters via method-of-moments. Note that `sd`
must satisfy $`\sigma < \sqrt{\mu(1-\mu)}`$ to ensure a valid (unimodal)
beta density.

In the example below, stray rates are estimated independently per
population using a single time block, with population 1 given a prior
centred at 0.1 and population 2 centred at 0.5:

``` r


input_list$data$n_pop <- 2 # modify to two populations for demo purposes

stray_prior <- data.frame(
  pop   = c(1, 2),
  block = c(1, 1),
  mu    = c(0.1, 0.5),
  sd    = c(0.1, 0.1)
)

input_list <- Setup_Mod_Rec(
  input_list       = input_list,
  do_rec_bias_ramp = 0,
  rec_model        = "bh_rec",
  rec_dd           = "local",
  stray_rate_blocks = c(
    "Block_1_Year_1-terminal_Pop_1",
    "Block_1_Year_1-terminal_Pop_2"
  ),
  stray_rate_spec      = "est_all",
  use_fixed_stray_rate = 0,
  use_stray_rate_prior = 1,
  stray_rate_prior     = stray_prior
)
```

Tighter values of `sd` produce a stronger regularizing effect and are
recommended when data contain little information about cross-population
spawning contributions. If `sd` is set too large relative to
`mu * (1 - mu)`, the beta density becomes U-shaped (mass near 0 and 1),
which can cause numerical instability during optimization; `SPoRC`
guards against this internally by squishing the logistic transform away
from the boundaries, but specifying a sensible `sd` is still good
practice.

### Time-varying stray rates

Stray rates can also be allowed to vary over time by defining multiple
blocks. The `stray_rate_blocks` argument follows the same blocking
syntax used elsewhere in `SPoRC`. For example, to allow stray rates to
differ between an early and late period for each population:

``` r

stray_prior_tv <- data.frame(
  pop   = c(1, 1, 2, 2),
  block = c(1, 2, 1, 2),
  mu    = c(0.1, 0.2, 0.4, 0.5),
  sd    = c(0.05, 0.05, 0.1, 0.1)
)

input_list <- Setup_Mod_Rec(
  input_list       = input_list,
  do_rec_bias_ramp = 0,
  rec_model        = "bh_rec",
  rec_dd           = "local",
  stray_rate_blocks = c(
    "Block_1_Year_1-20_Pop_1",
    "Block_2_Year_21-terminal_Pop_1",
    "Block_1_Year_1-20_Pop_2",
    "Block_2_Year_21-terminal_Pop_2"
  ),
  stray_rate_spec      = "est_all",
  use_fixed_stray_rate = 0,
  use_stray_rate_prior = 1,
  stray_rate_prior     = stray_prior_tv
)
```

Each unique block index in `stray_rate_blocks` requires a corresponding
row in `stray_rate_prior`, one per population × block combination.
