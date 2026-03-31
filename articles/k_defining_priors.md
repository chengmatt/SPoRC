# Defining Priors

In some cases, prior information might be available for use to inform
model estimation. Priors can either serve as informative priors that
have a strong influence on model estimation (e.g., natural mortality,
catchability), but can also serve as regularizing priors that help
stabilize model estimation. These regularizing priors can be helpful
particularly for stabalizing tag reporting, movement, and selectivity
parameters because some parameters are not well informed by data (e.g.,
slope of selectivity when it approaches infinity) such that the gradient
is zero, with respect to the likelihood surface.

In `SPoRC`, priors can be placed on several parameters. These include:

1.  Steepness,
2.  Natural Mortality,
3.  Catchability,
4.  Selectivity,
5.  Movement,
6.  Tag Reporting Rates, and
7.  Recruitment Proportions

In the following, we will demonstrate how priors can be specified for
each of these parameters. For further mathematical details on the
formulation of these priors, refer to the vignette describing model
equations.

Let us first load in any required packages, data we may use, and define
model dimensions.

``` r
# Load in packages
library(SPoRC) 
#> Loading required package: RTMB
data("sgl_rg_sable_data") # load in data

# Define model dimensions
input_list <- Setup_Mod_Dim(years = 1:length(sgl_rg_sable_data$years), # vector of years 
                            # (corresponds to year 1960 - 2024)
                            ages = 1:length(sgl_rg_sable_data$ages), # vector of ages
                            lens = seq(41,99,2), # number of lengths
                            n_regions = 1, # number of regions
                            n_sexes = sgl_rg_sable_data$n_sexes, # number of sexes == 1,
                            # female, == 2 male
                            n_fish_fleets = sgl_rg_sable_data$n_fish_fleets, # number of fishery
                            # fleet == 1, fixed gear, == 2 trawl gear
                            n_srv_fleets = sgl_rg_sable_data$n_srv_fleets, # number of survey fleets
                            verbose = TRUE
                            )
#> Number of Years: 65
#> Number of Projection Years for Dev Pars: 0
#> Number of Regions: 1
#> Number of Age Bins: 30
#> Number of Length Bins: 30
#> Number of Sexes: 2
#> Number of Fishery Fleets: 2
#> Number of Survey Fleets: 3
```

## Steepness

To specify priors on the steepness of a Beverton-Holt stock-recruitment
relationship, we use the `Setup_Mod_Rec` function with default settings
for simplicity. First, we define a dataframe that includes the region
where the prior applies, the mean (`mu`), and standard deviation (`sd`)
of the steepness prior. These values are specified in normal space but
are internally transformed to a beta distribution bounded between 0.2
and 1. To activate steepness priors, set `Use_h_prior = 1` and provide
the dataframe to the `h_prior` argument:

``` r
# set up dataframe specifying the steepness prior and the dimensions it is applied to
steepness_prior <- data.frame(
  region = 1, # apply steepness prior to region 1
  mu = 0.6, # mean of steepness
  sd = 0.2 # sd of steepness
)

# setup recruitment options
input_list <- Setup_Mod_Rec(
  input_list = input_list,        # input data list from above
  # Model options
  rec_model = "bh_rec",          # recruitment model type
  Use_h_prior = 1,               # whether steepness priors are used
  h_prior = steepness_prior      # Steepness prior dataframe
)
#> Recruitment is specified as: bh_rec
#> Recruitment and SSB lag is specified as: 1
#> Recruitment proportion priors are: Not Used
#> Steepness priors are: Used
#> Sex Ratios estimated with 1 block for region 1
#> Recruitment Bias Ramp is: Off
#> Initial Age Structure is: Movement and Matrix Geometric Series
#> Recruitment deviations for every year are estimated
#> Recruitment Variability is estimated for both early and late periods
#> Initial Age Deviations is stochastic for all ages, but the plus group follows equilibrium calculations.
#> Recruitment Deviations is estimated for all dimensions
#> Steepness is estimated for all dimensions
#> Sex ratio is specified as: fix
```

In a spatial model where region-specific steepness priors are required,
the prior dataframe can include multiple entries. In particular, each
row corresponds to a region-specific steepness parameter, along with its
associated mean (`mu`) and standard deviation (`sd`). For instance:

``` r
steepness_prior <- data.frame(
  region = c(1,2,3), # apply steepness prior to region 1, 2 and, 3
  mu = c(0.6, 0.8, 0.3), # mean of steepness prior corresponding to the regions defined
  sd = c(0.2, 0.2, 0.2) # sd of steepness corresponding to the regions defined
)

steepness_prior
#>   region  mu  sd
#> 1      1 0.6 0.2
#> 2      2 0.8 0.2
#> 3      3 0.3 0.2
```

## Natural Mortality

To define natural mortality priors, we use the `Setup_Mod_Biologicals`
function to specify a natural mortality prior. To do so, we will set
`Use_M_prior = 1` and provide a dataframe of values corresponding the
region, year, age, and sex blocks specified for natural mortality, along
with their associated mean and standard deviation.

``` r
# Define natural mortality prior
M_prior <- data.frame( 
  regionblk = 1, # region block
  yearblk = 1, # year block
  ageblk = c(1), # age block
  sexblk = c(1,2), # sex block
  mu = 0.085, # mean
  sd = 0.1 #sd
)

input_list <- Setup_Mod_Biologicals(input_list = input_list,
                                    # Data inputs
                                    WAA = sgl_rg_sable_data$WAA, # weight-at-age
                                    MatAA = sgl_rg_sable_data$MatAA, # maturity at age
                                    AgeingError = as.matrix(sgl_rg_sable_data$age_error), # ageing error
                                    SizeAgeTrans = sgl_rg_sable_data$SizeAgeTrans, # size age transition matrix
                                    # Model options
                                    Use_M_prior = 1, # use natural mortality prior
                                    M_prior = M_prior, # mean and sd for M prior
                                    M_spec = "est_ln_M", # estmating M
                                    M_sexblk_spec = list(1,2) # specifying M as sex-specific
)
#> Length Composition data are: Not Used
#> Natural Mortality priors are: Used
#> Selectivity is aged-based.
#> WAA_fish was specified at NULL. Using the spawning WAA for WAA_fish
#> WAA_srv was specified at NULL. Using the spawning WAA for WAA_srv
#> Ageing Error is specified to be time-invariant
#> Natural Mortality specified as: est_ln_M
#> Natural Mortality Region Blocks is specified as: 1
#> Natural Mortality Year Blocks is specified as: 1
#> Natural Mortality Age Blocks is specified as: 1
#> Natural Mortality Sex Blocks is specified as: 2
```

## Catchability

Defining catchability priors for the fishery and survey follows the same
process. For demonstration purposes, the following section will
illustrate how catchability priors can be defined for the survey using
the `Setup_Mod_Srvsel_and_Q` function. Similar to the steepness section,
a dataframe of catchability priors is first defined that determines the
region, fleet, and block the catchability prior is applied to. The
dataframe also includes columns defining the the mean (`mu`), and
standard deviation (`sd`) of the catchability prior. Again, these priors
are provided in normal space, but are internally transformed to be a
lognormal density. The `Use_srv_q_prior = 1` argument is utilized to
activate the prior, and the dataframe is provided to the `srv_q_prior`
argument.

``` r
# set up dataframe specifying the survey catchability prior and the dimensions it is applied to
srv_q_prior <- data.frame(
  region = 1, # apply catchability prior to region 1
  fleet = c(1, 2, 3), # apply catchability prior to region 1, fleets 1 - 3
  block = 1, # apply catchability prior to region 1, fleets 1 - 3, block 1
  mu = c(6, 0.5, 4), # corresponding mean catchability prior values
  sd = c(0.1, 0.35, 0.1) # corresponding standard deviation prior values
)

srv_q_prior
#>   region fleet block  mu   sd
#> 1      1     1     1 6.0 0.10
#> 2      1     2     1 0.5 0.35
#> 3      1     3     1 4.0 0.10

# Setup survey selectivity and catchability
input_list <- Setup_Mod_Srvsel_and_Q(input_list = input_list,

                                     # Model options
                                     # survey selectivity, whether continuous time-varying
                                     cont_tv_srv_sel = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3"),
                                     # survey selectivity blocks
                                     srv_sel_blocks = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3"), 
                                     # survey selectivity form
                                     srv_sel_model = c("logist1_Fleet_1", "exponential_Fleet_2", "logist1_Fleet_3"),
                                     # survey catchability blocks
                                     srv_q_blocks = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3"),
                                     # whether to estiamte all fixed effects for survey selectivity
                                     srv_fixed_sel_pars_spec = c("est_all", "est_all", "est_all"),
                                     # whether to estiamte all fixed effects for survey catchability
                                     srv_q_spec = c("est_all", "est_all", "est_all"),
                                     # Using catchability priors
                                     Use_srv_q_prior = 1,
                                     # Input survey catchability prior dataframe
                                     srv_q_prior = srv_q_prior
)
#> Survey Catchability priors are: Used
#> Survey Selectivity priors are: Not Used
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
`block`, `sex`, `fleet`, `mu`, `sd`, which defines the region,
parameter, time block, sex, and fleet combination to apply the prior to.
Given the complexity of the number of selectivity parameters that can be
estimated across model partitions (regions, sexes, blocks, and fleets),
users will need to be cautious how the selectivity prior dataframe is
defined to ensure they are correctly applied. In the following, we will
apply priors to the logistic selectivity parameters for both sexes in
fleet 1, and priors to the power parameter describing exponential
selectivity in fleet 2.

``` r
srv_selex_prior <- data.frame(
  region = 1, # only a single region
  par = c(1,2,1,2,1,1), # parameter to apply prior to
  sex = c(1,1,2,2,1,2), # sexes to apply prior to
  fleet = c(1,1,1,1,2,2), # fleets to apply prior to
  block = 1, # only a single time block
  mu = rep(0, 6), # mean specified at 0
  sd = rep(5, 6) # sd specified to be wide for all parameters
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
that a prior is pplied to region 1, parameter 1, sex 1, fleet 1, and
block 1, which corresponds to the $b_{50}$ parameter for logistic
selectivity specified for the first survey fleet. The last row of the
datframe indicates that a prior is applied to region 1, parameter 1, sex
2, fleet 2, and block 1, which corresponds to the male $\phi$ power
paramter for exponential selectivity specified for the second survey
fleet.

Next, we can activate selectivity priors by setting
`Use_srv_selex_prior = 1` and providing the associated selectivity prior
dataframe to `srv_selex_prior`.

``` r
# Setup survey selectivity and catchability
input_list <- Setup_Mod_Srvsel_and_Q(input_list = input_list,

                                     # Model options
                                     # survey selectivity, whether continuous time-varying
                                     cont_tv_srv_sel = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3"),
                                     # survey selectivity blocks
                                     srv_sel_blocks = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3"), 
                                     # survey selectivity form
                                     srv_sel_model = c("logist1_Fleet_1", "exponential_Fleet_2", "logist1_Fleet_3"),
                                     # survey catchability blocks
                                     srv_q_blocks = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3"),
                                     # whether to estiamte all fixed effects for survey selectivity
                                     srv_fixed_sel_pars_spec = c("est_all", "est_all", "est_all"),
                                     # whether to estiamte all fixed effects for survey catchability
                                     srv_q_spec = c("est_all", "est_all", "est_all"),
                                     # Using selex priors
                                     Use_srv_selex_prior = 1,
                                     # Input survey selex priors
                                     srv_selex_prior = srv_selex_prior
)
#> Survey Catchability priors are: Not Used
#> Survey Selectivity priors are: Used
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
proportions (see section below) is best illustrated with setting up a
spatial model. Here, we will load in data from the multi-region
sablefish case study to demonstrate defining priors for these processes.

``` r
data(mlt_rg_sable_data) # load in data

# Initialize model dimensions and data list
input_list <- Setup_Mod_Dim(years = 1:length(mlt_rg_sable_data$years),
                            # vector of years (1 - 62)
                            ages = 1:length(mlt_rg_sable_data$ages),
                            # vector of ages (1 - 30)
                            lens = mlt_rg_sable_data$lens,
                            # number of lengths (41 - 99)
                            n_regions = mlt_rg_sable_data$n_regions,
                            # number of regions (5)
                            n_sexes = mlt_rg_sable_data$n_sexes,
                            # number of sexes (2)
                            n_fish_fleets = mlt_rg_sable_data$n_fish_fleets,
                            # number of fishery fleet (2)
                            n_srv_fleets = mlt_rg_sable_data$n_srv_fleets,
                            # number of survey fleets (2)
                            verbose = TRUE
)
#> Number of Years: 62
#> Number of Projection Years for Dev Pars: 0
#> Number of Regions: 5
#> Number of Age Bins: 30
#> Number of Length Bins: 30
#> Number of Sexes: 2
#> Number of Fishery Fleets: 2
#> Number of Survey Fleets: 2
```

Because movement parameters are correlated and must sum to 1, `SPoRC`
uses a Dirichlet prior to model movement probabilities. To illustrate,
we set up the movement specification using the `Setup_Mod_Movement`
function. Movement priors are enabled by setting
`Use_Movement_Prior = 1`. The `Movement_prior` argument must be provided
as a dataframe where each row corresponds to an origin region (and age
or sex index) and includes a list-column alpha specifying the Dirichlet
shape parameters that define the prior over movement probabilities to
all destination regions.

``` r
prior <- expand.grid(
    region_from = 1:input_list$data$n_regions, # regions
    age = 1, # age blocks
    sex = 1, # sex
    alpha = I(list(rep(2, input_list$data$n_regions))) # prior alpha to each row
  )

# setting up movement parameterization
input_list <- Setup_Mod_Movement(input_list = input_list,
                                 # Model options
                                 Movement_ageblk_spec = "constant",
                                 # estimating movement in 3 age blocks
                                 # (ages 1-6, ages 7-15, ages 16-30)
                                 Movement_yearblk_spec = "constant", # time-invariant movement
                                 Movement_sexblk_spec = "constant", # sex-invariant movement
                                 Use_Movement_Prior = 1, # priors used for movement
                                 Movement_prior = prior  # vague prior to penalize movement away from the extremes
                                 )
#> Movement is: Estimated
#> Movement priors are: Used
#> Recruits are: Not Moving
#> Continuous movement specification is: none
#> Continuous movement process error specification is: none
#> Movement type is: Unstructured Markov
#> Movement fixed effect blocks are sex-invariant
#> Movement fixed effect blocks are time-invariant
#> Movement fixed effect blocks are age-invariant
```

## Tag Reporting Rates

Defining priors for tag reporting rates is done by constructing a
dataframe and passing it to the `Setup_Mod_Tagging` function. The
dataframe must contain the following columns: `region`, `block`, `mu`,
`sd`, and `type.` These columns define the region and time block the
prior applies to, the mean and standard deviation of the prior
distribution, and the type of beta prior to use. Two prior types are
supported: a symmetric beta prior (`type = 0`), which penalizes extreme
values near 0 and 1, and a standard beta prior (`type = 1`), which
allows users to specify both the mean and standard deviation directly.

In the example below, symmetric beta priors are applied to two time
blocks in region 1. Because the symmetric form does not require a mean,
the `mu` column is set to NA:

``` r
# setup tagging priors
tag_prior <- data.frame(
  region = 1,         # apply to region 1
  block  = c(1, 2),   # two time blocks
  mu     = NA,        # not needed for symmetric beta
  sd     = 5,         # standard deviation
  type   = 0          # symmetric beta prior
)
tag_prior
#>   region block mu sd type
#> 1      1     1 NA  5    0
#> 2      1     2 NA  5    0
```

In this example, the first row of the dataframe indicates that a
symmetric beta prior is applied to tag reporting rates in region 1 and
block 1, with a standard deviation of 5. The second row applies the same
prior structure to block 2. Because `type = 0`, the prior is symmetric
around 0.5 and discourages extreme reporting rate estimates. If
`type = 1` were used instead, the mu column would need to be specified
with the prior mean. Tag reporting rate priors can then be applied by
setting `Use_TagRep_Prior = 1` and inputting the dataframe into the
`TagRep_Prior` argument.

``` r
# setting up tagging parameterization
input_list <- Setup_Mod_Tagging(input_list = input_list,
                                UseTagging = 1, # using tagging data
                                max_tag_liberty = 15, # maximum number of years to track a cohort

                                # Data Inputs
                                tag_release_indicator = mlt_rg_sable_data$tag_release_indicator,
                                Tagged_Fish = mlt_rg_sable_data$Tagged_Fish, # Released fish
                                Obs_Tag_Recap = mlt_rg_sable_data$Obs_Tag_Recap,

                                # Model options
                                Use_TagRep_Prior = 1, # tag reporting rate priors are used
                                TagRep_Prior = tag_prior,
                                # sex-specific data when fitting
                                Init_Tag_Mort_spec = "fix", # fixing initial tag mortality
                                Tag_Shed_spec = "fix", # fixing chronic shedding
                                TagRep_spec = "est_shared_r", # tag reporting rates are
                                # not region specific
                                # Time blocks for tag reporting rates
                                Tag_Reporting_blocks = c(
                                  paste("Block_1_Year_1-35_Region_",
                                        c(1:input_list$data$n_regions), sep = ''),
                                  paste("Block_2_Year_36-terminal_Region_",
                                        c(1:input_list$data$n_regions), sep = '')
                                ),

                                # Specify starting values or fixing values
                                ln_Init_Tag_Mort = log(0.1), # fixing initial tag mortality
                                ln_Tag_Shed = log(0.02)  # fixing tag shedding
)
#> Tagging priors are used
#> Tagging data are fit to 30 age groups
#> Tagging data are fit to 2 sex groups
#> Tag Reporting estimated with 2 block for region 1
#> Tag Reporting estimated with 2 block for region 2
#> Tag Reporting estimated with 2 block for region 3
#> Tag Reporting estimated with 2 block for region 4
#> Tag Reporting estimated with 2 block for region 5
#> Initial Tag Mortality is specified as: fix
#> Chronic Tag Shedding is specified as: fix
```

## Recruitment Proportions

In the context of a spatial model, Dirichlet priors can also be placed
on recruitment proportions. This process is relatively straightforward
and uses the `Setup_Mod_Rec` functions to do so. In particular, users
will need to define whether the prior is being used by setting
`Use_Rec_prop_Prior = 1`. Following that, a scalar or vector of
(`n_regions`) can be defined for `Rec_prop_prior`. In the case where
`Rec_prop_prior` is a scalar, a uniform Dirichlet prior will be applied.
If a vector of different values for each region are provided, then the
Dirichlet prior will be applied according to those concentration
parameters.

``` r
# Setup recruitment stuff (using defaults for other stuff)
input_list <- Setup_Mod_Rec(input_list = input_list, # input data list from above
                            do_rec_bias_ramp = 0,

                            # Model options
                            rec_model = "mean_rec", # recruitment model
                            Use_Rec_prop_Prior = 1, # using recruitment proportion prior
                            Rec_prop_prior = 5 # prior value for recruitment proportions
)
#> Recruitment is specified as: mean_rec
#> Recruitment proportion priors are: Used
#> Sex Ratios estimated with 1 block for region 1
#> Sex Ratios estimated with 1 block for region 2
#> Sex Ratios estimated with 1 block for region 3
#> Sex Ratios estimated with 1 block for region 4
#> Sex Ratios estimated with 1 block for region 5
#> Recruitment Bias Ramp is: Off
#> Initial Age Structure is: Movement and Matrix Geometric Series
#> Recruitment deviations for every year are estimated
#> Recruitment Variability is estimated for both early and late periods
#> Initial Age Deviations is stochastic for all ages, but the plus group follows equilibrium calculations.
#> Recruitment Deviations is estimated for all dimensions
#> Sex ratio is specified as: fix
```
